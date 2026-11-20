#!/usr/bin/env python3
"""Call active promoters and enhancers from cleaned histone-mark peaks.

Promoters: H3K4me3 peaks that lie within --tss-window bp of an annotated TSS.
  Each promoter lists the gene(s) whose TSS it covers, and whether it is also
  marked by H3K27ac (active).

Enhancers: H3K27ac peaks farther than --tss-window bp from every TSS
  (optionally also not overlapping any H3K4me3 peak). Each enhancer records
  whether it overlaps the H3K4me1 proxy; with --require-h3k4me1 only those are kept.

Outputs (BED-like, with a header row):
  promoters: chrom start end promoter_id genes n_genes h3k27ac signal
  enhancers: chrom start end enhancer_id h3k4me1_proxy signal

Example:
  python scripts/call_elements.py \
      --h3k27ac data/processed/peaks/bmec_h3k27ac_peaks.bed \
      --h3k4me3 data/processed/peaks/bmec_h3k4me3_peaks.bed \
      --h3k4me1 data/processed/peaks/huvec_h3k4me1_peaks.bed \
      --tss data/processed/annotation/tss.bed --tss-window 2500 \
      --promoters-out data/processed/elements/promoters.tsv \
      --enhancers-out data/processed/elements/enhancers.tsv
"""

import argparse

from bed_utils import (NARROWPEAK, bedtools_df, read_bed, sort_bed, temp_bed,
                       windows, write_table)

TSS_COLUMNS = ["chrom", "start", "end", "gene_name", "gene_id", "strand"]


def overlaps(a, b_path):
    """Boolean Series: does each row of a overlap at least one interval in b_path?"""
    with temp_bed(a[["chrom", "start", "end"]].assign(row=a.index)) as path:
        hits = bedtools_df("intersect", "-u", "-a", path, "-b", b_path)
    hit_rows = set(hits[3]) if not hits.empty else set()
    return a.index.to_series().isin(hit_rows)


def call_promoters(h3k4me3, tss_windows, h3k27ac_path):
    """H3K4me3 peaks overlapping a TSS window, with the genes they cover."""
    with temp_bed(h3k4me3[["chrom", "start", "end"]].assign(row=h3k4me3.index)) as a, \
         temp_bed(tss_windows) as b:
        hits = bedtools_df("intersect", "-wa", "-wb", "-a", a, "-b", b)
    # columns: 0-2 peak, 3 row id, 4-9 TSS window (gene_name is column 7)
    genes = (hits.groupby(3)[7]
             .agg(lambda g: ",".join(sorted(set(g))))
             .rename("genes"))
    prom = h3k4me3.loc[genes.index, ["chrom", "start", "end", "signal"]].join(genes)
    prom["n_genes"] = prom["genes"].str.count(",") + 1
    prom["h3k27ac"] = overlaps(prom, h3k27ac_path)
    prom = sort_bed(prom)
    prom["promoter_id"] = [f"P{i:06d}" for i in range(1, len(prom) + 1)]
    return prom[["chrom", "start", "end", "promoter_id", "genes", "n_genes", "h3k27ac", "signal"]]


def call_enhancers(h3k27ac, tss_windows_path, h3k4me3_path, h3k4me1_path,
                   exclude_h3k4me3, require_h3k4me1):
    """Distal H3K27ac peaks, flagged for overlap with the H3K4me1 proxy."""
    enh = h3k27ac[~overlaps(h3k27ac, tss_windows_path)]
    n_distal = len(enh)
    if exclude_h3k4me3:
        enh = enh[~overlaps(enh, h3k4me3_path)]
    n_no_k4me3 = len(enh)
    enh = enh.assign(h3k4me1_proxy=overlaps(enh, h3k4me1_path))
    if require_h3k4me1:
        enh = enh[enh["h3k4me1_proxy"]]
    enh = sort_bed(enh[["chrom", "start", "end", "h3k4me1_proxy", "signal"]])
    enh["enhancer_id"] = [f"E{i:06d}" for i in range(1, len(enh) + 1)]
    counts = {"distal_h3k27ac": n_distal, "without_h3k4me3": n_no_k4me3, "kept": len(enh)}
    return enh[["chrom", "start", "end", "enhancer_id", "h3k4me1_proxy", "signal"]], counts


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--h3k27ac", required=True)
    parser.add_argument("--h3k4me3", required=True)
    parser.add_argument("--h3k4me1", required=True, help="H3K4me1 peaks (proxy)")
    parser.add_argument("--tss", required=True)
    parser.add_argument("--tss-window", type=int, required=True)
    parser.add_argument("--exclude-h3k4me3", action="store_true")
    parser.add_argument("--require-h3k4me1", action="store_true")
    parser.add_argument("--promoters-out", required=True)
    parser.add_argument("--enhancers-out", required=True)
    args = parser.parse_args()

    h3k27ac = read_bed(args.h3k27ac, names=NARROWPEAK)
    h3k4me3 = read_bed(args.h3k4me3, names=NARROWPEAK)
    tss_windows = windows(read_bed(args.tss, names=TSS_COLUMNS), args.tss_window)

    promoters = call_promoters(h3k4me3, tss_windows, args.h3k27ac)
    with temp_bed(tss_windows) as tss_path:
        enhancers, counts = call_enhancers(h3k27ac, tss_path, args.h3k4me3, args.h3k4me1,
                                           args.exclude_h3k4me3, args.require_h3k4me1)

    write_table(promoters, args.promoters_out)
    write_table(enhancers, args.enhancers_out)
    print(f"promoters: {len(promoters)} (H3K27ac+: {int(promoters['h3k27ac'].sum())}) "
          f"-> {args.promoters_out}")
    print(f"enhancers: {counts['distal_h3k27ac']} distal H3K27ac, "
          f"{counts['without_h3k4me3']} without H3K4me3, {counts['kept']} kept "
          f"(H3K4me1 proxy+: {int(enhancers['h3k4me1_proxy'].sum())}) -> {args.enhancers_out}")


if __name__ == "__main__":
    main()
