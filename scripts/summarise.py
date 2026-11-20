#!/usr/bin/env python3
"""Write the QC and summary tables the report reads.

peak_qc.tsv:   one row per cleaned peak file: mark, n_peaks, median_width, total_mb,
               frac_near_tss (share of peaks within --tss-window bp of a TSS)
summary.tsv:   one row per headline number (metric, value), so the report can
               quote every count inline without recomputing it

Example:
  python scripts/summarise.py \
      --peaks bmec_h3k27ac=data/processed/peaks/bmec_h3k27ac_peaks.bed ... \
      --tss data/processed/annotation/tss.bed --tss-window 2500 \
      --promoters data/processed/elements/promoters.tsv \
      --enhancers data/processed/elements/enhancers.tsv \
      --pairs results/tables/ep_pairs.tsv \
      --context results/tables/element_context.tsv \
      --tfs results/tables/candidate_tfs.tsv \
      --peak-qc-out results/tables/peak_qc.tsv \
      --summary-out results/tables/summary.tsv
"""

import argparse

import pandas as pd

from bed_utils import NARROWPEAK, bedtools_df, read_bed, temp_bed, windows, write_table

TSS_COLUMNS = ["chrom", "start", "end", "gene_name", "gene_id", "strand"]


def peak_qc(name, path, tss_windows_path):
    peaks = read_bed(path, names=NARROWPEAK)
    width = peaks["end"] - peaks["start"]
    near = bedtools_df("intersect", "-u", "-a", path, "-b", tss_windows_path)
    return {"mark": name, "n_peaks": len(peaks), "median_width": int(width.median()),
            "total_mb": round(width.sum() / 1e6, 1), "frac_near_tss": round(len(near) / len(peaks), 3)}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--peaks", nargs="+", required=True, help="name=path pairs")
    parser.add_argument("--tss", required=True)
    parser.add_argument("--tss-window", type=int, required=True)
    parser.add_argument("--promoters", required=True)
    parser.add_argument("--enhancers", required=True)
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--context", required=True)
    parser.add_argument("--tfs", required=True)
    parser.add_argument("--peak-qc-out", required=True)
    parser.add_argument("--summary-out", required=True)
    args = parser.parse_args()

    tss_windows = windows(read_bed(args.tss, names=TSS_COLUMNS), args.tss_window)
    with temp_bed(tss_windows) as tss_path:
        qc = pd.DataFrame([peak_qc(*item.split("=", 1), tss_path) for item in args.peaks])

    prom = pd.read_csv(args.promoters, sep="\t")
    enh = pd.read_csv(args.enhancers, sep="\t")
    pairs = pd.read_csv(args.pairs, sep="\t")
    ctx = pd.read_csv(args.context, sep="\t", keep_default_na=False)
    tfs = pd.read_csv(args.tfs, sep="\t", keep_default_na=False)
    linked = enh["enhancer_id"].isin(pairs["enhancer_id"])

    summary = pd.DataFrame([
        ("n_promoters", len(prom)),
        ("n_promoters_h3k27ac", int(prom["h3k27ac"].sum())),
        ("n_enhancers", len(enh)),
        ("n_enhancers_h3k4me1_proxy", int(enh["h3k4me1_proxy"].sum())),
        ("n_ep_pairs", len(pairs)),
        ("n_loop_enhancers", int(linked.sum())),
        ("n_loop_enhancers_h3k4me1_proxy", int(enh.loc[linked, "h3k4me1_proxy"].sum())),
        ("n_linked_promoters", pairs["promoter_id"].nunique()),
        ("n_target_genes", pairs["genes"].str.split(",").explode().nunique()),
        ("median_ep_distance", int(pairs["distance"].median()) if len(pairs) else 0),
        ("frac_enhancers_in_A", round((ctx.loc[ctx["type"] == "enhancer", "compartment"] == "A").mean(), 3)),
        ("frac_promoters_in_A", round((ctx.loc[ctx["type"] == "promoter", "compartment"] == "A").mean(), 3)),
        ("n_candidate_tfs", len(tfs)),
        ("n_candidate_tf_families", tfs.loc[tfs["family"] != "", "family"].nunique()),
    ], columns=["metric", "value"], dtype=object)  # object keeps counts as integers

    write_table(qc, args.peak_qc_out)
    write_table(summary, args.summary_out)
    print(qc.to_string(index=False))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
