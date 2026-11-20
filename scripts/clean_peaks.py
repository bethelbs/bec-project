#!/usr/bin/env python3
"""Clean one ENCODE narrowPeak file: keep chosen chromosomes, drop blacklist overlaps,
apply a q-value threshold, remove duplicate intervals and sort.

Output keeps all ten narrowPeak columns.

Example:
  python scripts/clean_peaks.py \
      --peaks data/raw/chip/bmec_h3k27ac_peaks.bed.gz \
      --blacklist data/raw/reference/grch38_blacklist.bed.gz \
      --chromosomes chr1 chr2 ... --min-neglog10-q 0 \
      --out data/processed/peaks/bmec_h3k27ac_peaks.bed
"""

import argparse

from bed_utils import NARROWPEAK, bedtools_df, read_bed, sort_bed, temp_bed, write_bed


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--peaks", required=True)
    parser.add_argument("--blacklist", required=True)
    parser.add_argument("--chromosomes", nargs="+", required=True)
    parser.add_argument("--min-neglog10-q", type=float, default=0.0)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    peaks = read_bed(args.peaks, names=NARROWPEAK)
    n_raw = len(peaks)

    peaks = peaks[peaks["chrom"].isin(args.chromosomes)]
    n_chrom = len(peaks)

    peaks = peaks[peaks["neglog10_q"] >= args.min_neglog10_q]
    n_q = len(peaks)

    # bedtools intersect -v keeps peaks with no overlap at all with the blacklist.
    with temp_bed(peaks) as path:
        kept = bedtools_df("intersect", "-v", "-a", path, "-b", args.blacklist)
    if kept.empty:
        raise SystemExit(f"ERROR: no peaks left in {args.peaks} after filtering")
    kept.columns = NARROWPEAK
    n_black = len(kept)

    kept = sort_bed(kept.drop_duplicates(["chrom", "start", "end"]))

    write_bed(kept, args.out)
    print(f"{args.peaks}: {n_raw} raw, {n_chrom} on kept chromosomes, {n_q} pass q, "
          f"{n_black} outside blacklist, {len(kept)} unique -> {args.out}")


if __name__ == "__main__":
    main()
