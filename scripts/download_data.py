#!/usr/bin/env python3
"""Download the raw input data listed in metadata/sources.tsv.

Uses only the Python standard library, so it runs before the project's
Pixi environment exists.

For each row of metadata/sources.tsv the script:
  1. skips the file if it is already present and matches its SHA-256 in
     metadata/checksums.sha256;
  2. otherwise removes any stale copy, downloads the file into a temporary
     folder inside data/, and checks it against the md5 published by ENCODE;
  3. unpacks it if the `unpack` column says `tar`, then moves it into place;
  4. makes the file read-only.

Finally it verifies every file listed in metadata/checksums.sha256.

Usage (from anywhere; paths are resolved from this script's location):
  python scripts/download_data.py                    # fetch what is missing or wrong
  python scripts/download_data.py --verify-only      # check only, download nothing
  python scripts/download_data.py --write-checksums  # (re)create checksums.sha256
"""

import argparse
import csv
import hashlib
import shutil
import stat
import sys
import tarfile
import tempfile
import time
import urllib.request
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SOURCES = PROJECT_ROOT / "metadata" / "sources.tsv"
CHECKSUMS = PROJECT_ROOT / "metadata" / "checksums.sha256"
DATA_DIR = PROJECT_ROOT / "data"

TIMEOUT_S = 60      # seconds without data before a download attempt fails
ATTEMPTS = 3        # tries per file
CHUNK = 1024 * 1024  # bytes read per step


def file_hash(path, algorithm):
    """Return the hex digest of a file, read in chunks."""
    h = hashlib.new(algorithm)
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()


def read_sources():
    """Read metadata/sources.tsv into a list of dicts."""
    with open(SOURCES, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    required = {"name", "url", "md5", "unpack", "dest"}
    missing = required - set(rows[0].keys() if rows else [])
    if missing:
        sys.exit(f"ERROR: {SOURCES} lacks columns: {', '.join(sorted(missing))}")
    return rows


def read_checksums():
    """Read metadata/checksums.sha256 into {relative path: sha256}."""
    if not CHECKSUMS.exists():
        return {}
    sums = {}
    for line in CHECKSUMS.read_text().splitlines():
        if line.strip():
            digest, rel = line.split(maxsplit=1)
            sums[rel.strip()] = digest
    return sums


def make_read_only(path):
    mode = path.stat().st_mode
    path.chmod(mode & ~(stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH))


def remove(path):
    """Delete a file or folder even if it was made read-only."""
    if path.is_dir():
        for p in path.rglob("*"):
            p.chmod(p.stat().st_mode | stat.S_IWUSR)
        path.chmod(path.stat().st_mode | stat.S_IWUSR)
        shutil.rmtree(path)
    elif path.exists():
        path.chmod(path.stat().st_mode | stat.S_IWUSR)
        path.unlink()


def kept_files(row):
    """Files a source row puts under data/raw (one file, or a tar's contents)."""
    dest = PROJECT_ROOT / row["dest"]
    if row["unpack"] == "tar":
        return sorted(p for p in dest.rglob("*") if p.is_file()) if dest.exists() else []
    return [dest] if dest.exists() else []


def is_current(row, sums):
    """True if every file from this source is present and matches its SHA-256."""
    files = kept_files(row)
    if not files:
        return False
    for f in files:
        rel = f.relative_to(PROJECT_ROOT).as_posix()
        if rel not in sums or file_hash(f, "sha256") != sums[rel]:
            return False
    return True


def download(url, target):
    """Download url to target, retrying on failure."""
    for attempt in range(1, ATTEMPTS + 1):
        try:
            with urllib.request.urlopen(url, timeout=TIMEOUT_S) as resp, open(target, "wb") as out:
                total = int(resp.headers.get("Content-Length") or 0)
                done = 0
                for block in iter(lambda: resp.read(CHUNK), b""):
                    out.write(block)
                    done += len(block)
                    if total:
                        print(f"\r    {done / 1e6:8.1f} / {total / 1e6:.1f} MB", end="", flush=True)
            print()
            return
        except OSError as err:  # includes URLError and timeouts
            print(f"\n    attempt {attempt}/{ATTEMPTS} failed: {err}")
            if attempt == ATTEMPTS:
                raise
            time.sleep(5 * attempt)


def fetch(row):
    """Download one source, check its md5, unpack if needed and move it into place."""
    dest = PROJECT_ROOT / row["dest"]
    for stale in kept_files(row) + ([dest] if dest.is_dir() else []):
        remove(stale)
    dest.parent.mkdir(parents=True, exist_ok=True)

    # Temporary files stay on the same disk as the data, so the final move is instant.
    with tempfile.TemporaryDirectory(dir=DATA_DIR, prefix=".download-") as tmp:
        tmp_file = Path(tmp) / Path(row["url"]).name
        print(f"  downloading {row['name']} from {row['url']}")
        download(row["url"], tmp_file)

        got = file_hash(tmp_file, "md5")
        if got != row["md5"]:
            sys.exit(f"ERROR: md5 mismatch for {row['name']}: expected {row['md5']}, got {got}")

        if row["unpack"] == "tar":
            unpacked = Path(tmp) / "unpacked"
            with tarfile.open(tmp_file) as tar:
                try:
                    tar.extractall(unpacked, filter="data")  # Python 3.12+: refuse unsafe paths
                except TypeError:
                    tar.extractall(unpacked)
            shutil.move(str(unpacked), dest)
        elif row["unpack"] == "none":
            shutil.move(str(tmp_file), dest)
        else:
            sys.exit(f"ERROR: unknown unpack value '{row['unpack']}' for {row['name']}")

    for f in kept_files(row):
        make_read_only(f)


def write_checksums(rows):
    """Write SHA-256 checksums of every kept raw file to metadata/checksums.sha256."""
    lines = []
    for row in rows:
        for f in kept_files(row):
            lines.append(f"{file_hash(f, 'sha256')}  {f.relative_to(PROJECT_ROOT).as_posix()}")
    CHECKSUMS.write_text("\n".join(lines) + "\n")
    print(f"Wrote {len(lines)} checksums to {CHECKSUMS.relative_to(PROJECT_ROOT)}")


def verify(rows, sums):
    """Check every listed file: present, SHA-256 matches, read-only. Return True if all pass."""
    problems = []
    expected = {f"{row['dest']}" for row in rows if row["unpack"] == "none"}
    unlisted = sorted(expected - set(sums))
    if unlisted:
        problems.append(f"not in {CHECKSUMS.name}: {', '.join(unlisted)}")
    for rel, digest in sums.items():
        path = PROJECT_ROOT / rel
        if not path.exists():
            problems.append(f"missing: {rel}")
        elif file_hash(path, "sha256") != digest:
            problems.append(f"sha256 mismatch: {rel}")
        elif path.stat().st_mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH):
            problems.append(f"writable (should be read-only): {rel}")
    if problems:
        print("FAILED:")
        for p in problems:
            print(f"  {p}")
        return False
    print(f"OK: all {len(sums)} raw files present, sha256 verified, read-only")
    return True


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--verify-only", action="store_true",
                       help="check the raw files without downloading anything")
    group.add_argument("--write-checksums", action="store_true",
                       help="write metadata/checksums.sha256 from the files now in data/raw")
    args = parser.parse_args()

    rows = read_sources()

    if args.write_checksums:
        write_checksums(rows)
        return

    sums = read_checksums()
    if not args.verify_only:
        if not sums:
            print(f"Note: {CHECKSUMS.name} not found; files are checked by md5 only.")
        for row in rows:
            if sums and is_current(row, sums):
                print(f"  ok: {row['name']}")
                continue
            dest = PROJECT_ROOT / row["dest"]
            if not sums and dest.exists():
                print(f"  present (no sha256 to check yet): {row['name']}")
                continue
            fetch(row)
        if not sums:
            print(f"Run with --write-checksums to create {CHECKSUMS.name}, then commit it.")
            return

    if not verify(rows, read_checksums()):
        sys.exit(1)


if __name__ == "__main__":
    main()

