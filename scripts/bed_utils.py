"""Small helpers shared by the analysis scripts: read and write BED files, run bedtools.

bedtools must be on PATH (it comes from the Pixi `pipeline` environment).
"""

import io
import subprocess
import tempfile
from contextlib import contextmanager
from pathlib import Path

import pandas as pd

NARROWPEAK = ["chrom", "start", "end", "name", "score", "strand",
              "signal", "neglog10_p", "neglog10_q", "summit"]


def read_bed(path, names=None):
    """Read a BED-like file (plain or .gz), skipping '#', 'track' and 'browser' lines."""
    df = pd.read_csv(path, sep="\t", header=None, comment="#", dtype={0: str})
    df = df[~df[0].str.startswith(("track", "browser"))].reset_index(drop=True)
    if names is not None:
        df = df.iloc[:, :len(names)]
        df.columns = names
    return df


def write_bed(df, path):
    """Write a DataFrame as a headerless, tab-separated BED file, creating parent folders."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", header=False, index=False)


def write_table(df, path):
    """Write a DataFrame as a TSV with a header row, creating parent folders."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def bedtools(*args):
    """Run a bedtools subcommand and return its standard output as text."""
    result = subprocess.run(["bedtools", *map(str, args)],
                            check=True, capture_output=True, text=True)
    return result.stdout


def bedtools_df(*args):
    """Run a bedtools subcommand and return its output as a DataFrame (empty if no output)."""
    out = bedtools(*args)
    if not out.strip():
        return pd.DataFrame()
    return pd.read_csv(io.StringIO(out), sep="\t", header=None, dtype={0: str})


@contextmanager
def temp_bed(df):
    """Write a DataFrame to a temporary BED file and yield its path; delete it afterwards."""
    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as fh:
        df.to_csv(fh, sep="\t", header=False, index=False)
        path = fh.name
    try:
        yield path
    finally:
        Path(path).unlink(missing_ok=True)


def sort_bed(df):
    """Sort a BED DataFrame by chromosome name, then start, then end."""
    return df.sort_values(["chrom", "start", "end"], kind="mergesort").reset_index(drop=True)


def windows(df, width):
    """Return a copy of df with start/end extended by `width` bp on each side (start clipped at 0)."""
    out = df.copy()
    out["start"] = (out["start"] - width).clip(lower=0)
    out["end"] = out["end"] + width
    return out
