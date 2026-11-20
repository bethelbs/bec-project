#!/usr/bin/env python3
"""Rank candidate transcription factors and target genes.

Candidate TFs: known motifs enriched in loop-supported enhancers (HOMER
  knownResults.txt), kept if their Benjamini q-value is at most --max-qvalue,
  ranked by p-value. The TF name and family are parsed from HOMER's motif name,
  e.g. "ERG(ETS)/VCaP-ERG-ChIP-Seq(GSE14097)/Homer" -> tf ERG, family ETS.

Target genes: genes whose promoters are linked to at least one enhancer by a
  loop, ranked by the number of distinct linked enhancers. Marker genes are
  flagged in an annotation column only.

Example:
  python scripts/rank_candidates.py \
      --known-results results/motifs/homer/knownResults.txt \
      --pairs results/tables/ep_pairs.tsv \
      --max-qvalue 0.05 --markers SLC2A1 CLDN5 ... \
      --tfs-out results/tables/candidate_tfs.tsv \
      --genes-out results/tables/target_genes.tsv
"""

import argparse

import pandas as pd

from bed_utils import write_table


def read_known_results(path):
    """Read HOMER knownResults.txt into tidy columns."""
    raw = pd.read_csv(path, sep="\t")
    raw.columns = [c.strip() for c in raw.columns]
    pct = lambda col: raw[col].astype(str).str.rstrip("%").astype(float)
    target_col = next(c for c in raw.columns if c.startswith("% of Target"))
    bg_col = next(c for c in raw.columns if c.startswith("% of Background"))
    df = pd.DataFrame({
        "motif_name": raw["Motif Name"],
        "consensus": raw["Consensus"],
        "pvalue": raw["P-value"],
        "log_pvalue": raw["Log P-value"],
        "qvalue": raw["q-value (Benjamini)"],
        "pct_targets": pct(target_col),
        "pct_background": pct(bg_col),
    })
    df["tf"] = df["motif_name"].str.split("/").str[0].str.replace(r"\(.*\)$", "", regex=True)
    df["family"] = df["motif_name"].str.split("/").str[0].str.extract(r"\(([^)]*)\)$", expand=False).fillna("")
    df["fold_enrichment"] = (df["pct_targets"] / df["pct_background"].where(df["pct_background"] > 0)).round(2)
    return df


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--known-results", required=True)
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--max-qvalue", type=float, required=True)
    parser.add_argument("--markers", nargs="*", default=[])
    parser.add_argument("--tfs-out", required=True)
    parser.add_argument("--genes-out", required=True)
    args = parser.parse_args()

    motifs = read_known_results(args.known_results)
    tfs = motifs[motifs["qvalue"] <= args.max_qvalue].sort_values("log_pvalue").reset_index(drop=True)
    tfs.insert(0, "rank", range(1, len(tfs) + 1))
    tfs = tfs[["rank", "tf", "family", "motif_name", "consensus", "pct_targets",
               "pct_background", "fold_enrichment", "pvalue", "qvalue"]]

    pairs = pd.read_csv(args.pairs, sep="\t")
    genes = pairs.assign(gene=pairs["genes"].str.split(",")).explode("gene")
    genes = (genes.groupby("gene")
             .agg(n_enhancers=("enhancer_id", "nunique"),
                  n_loops=("loop_ids", lambda s: len(set(",".join(s).split(",")))),
                  min_distance=("distance", "min"))
             .reset_index()
             .sort_values(["n_enhancers", "n_loops", "gene"], ascending=[False, False, True]))
    genes["is_marker"] = genes["gene"].isin(args.markers)

    write_table(tfs, args.tfs_out)
    write_table(genes, args.genes_out)
    families = ", ".join(tfs["family"].replace("", "?").head(10).unique())
    print(f"{len(tfs)} motifs with q <= {args.max_qvalue} (top families: {families}) -> {args.tfs_out}")
    print(f"{len(genes)} target genes, {int(genes['is_marker'].sum())} of them markers -> {args.genes_out}")


if __name__ == "__main__":
    main()
