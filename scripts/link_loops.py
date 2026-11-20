#!/usr/bin/env python3
"""Link enhancers to promoters through Hi-C loops.

Each loop has two anchors. A loop links an enhancer to a promoter when the
enhancer overlaps one anchor and the promoter overlaps the other. Anchors can be
widened with --anchor-padding.

Outputs:
  --pairs-out      one row per enhancer-promoter pair (TSV): enhancer_id, promoter_id,
                   genes, n_loops, distance (bp between element midpoints), loop_ids
  --enhancers-out  BED6 of the loop-supported enhancers, ready for HOMER
                   (chrom, start, end, enhancer_id, ".", "+")

Example:
  python scripts/link_loops.py \
      --loops data/raw/hic/bmec_hic_loops.bedpe.gz \
      --promoters data/processed/elements/promoters.tsv \
      --enhancers data/processed/elements/enhancers.tsv \
      --chromosomes chr1 chr2 ... --anchor-padding 0 \
      --pairs-out results/tables/ep_pairs.tsv \
      --enhancers-out data/processed/elements/loop_enhancers.bed
"""

import argparse

import pandas as pd

from bed_utils import bedtools_df, read_bed, sort_bed, temp_bed, windows, write_bed, write_table


def read_anchors(path, chromosomes, padding):
    """Return one row per loop anchor: chrom, start, end, loop_id, side (A or B)."""
    loops = read_bed(path).iloc[:, :6]
    loops.columns = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
    loops = loops[loops["chrom1"].isin(chromosomes) & (loops["chrom1"] == loops["chrom2"])]
    loops = loops.reset_index(drop=True)
    loops["loop_id"] = [f"L{i:05d}" for i in range(1, len(loops) + 1)]
    a = loops[["chrom1", "start1", "end1", "loop_id"]].set_axis(["chrom", "start", "end", "loop_id"], axis=1)
    b = loops[["chrom2", "start2", "end2", "loop_id"]].set_axis(["chrom", "start", "end", "loop_id"], axis=1)
    anchors = pd.concat([a.assign(side="A"), b.assign(side="B")], ignore_index=True)
    return windows(anchors, padding), len(loops)


def on_anchors(anchors, elements, id_col):
    """Which elements overlap which anchors: DataFrame with loop_id, side, element id."""
    with temp_bed(anchors) as a, temp_bed(elements[["chrom", "start", "end", id_col]]) as b:
        hits = bedtools_df("intersect", "-wa", "-wb", "-a", a, "-b", b)
    if hits.empty:
        return pd.DataFrame(columns=["loop_id", "side", id_col])
    # columns: 0-2 anchor, 3 loop_id, 4 side, 5-7 element, 8 element id
    return hits[[3, 4, 8]].set_axis(["loop_id", "side", id_col], axis=1).drop_duplicates()


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--loops", required=True)
    parser.add_argument("--promoters", required=True)
    parser.add_argument("--enhancers", required=True)
    parser.add_argument("--chromosomes", nargs="+", required=True)
    parser.add_argument("--anchor-padding", type=int, default=0)
    parser.add_argument("--pairs-out", required=True)
    parser.add_argument("--enhancers-out", required=True)
    args = parser.parse_args()

    anchors, n_loops = read_anchors(args.loops, args.chromosomes, args.anchor_padding)
    promoters = pd.read_csv(args.promoters, sep="\t")
    enhancers = pd.read_csv(args.enhancers, sep="\t")

    enh_hits = on_anchors(anchors, enhancers, "enhancer_id")
    prom_hits = on_anchors(anchors, promoters, "promoter_id")

    # Same loop, opposite anchors.
    pairs = enh_hits.merge(prom_hits, on="loop_id", suffixes=("_enh", "_prom"))
    pairs = pairs[pairs["side_enh"] != pairs["side_prom"]]

    pairs = (pairs.groupby(["enhancer_id", "promoter_id"])["loop_id"]
             .agg(n_loops="nunique", loop_ids=lambda s: ",".join(sorted(set(s))))
             .reset_index())

    mid = lambda df: (df["start"] + df["end"]) // 2
    enh_mid = enhancers.set_index("enhancer_id").pipe(mid)
    prom_mid = promoters.set_index("promoter_id").pipe(mid)
    pairs["distance"] = (pairs["enhancer_id"].map(enh_mid) - pairs["promoter_id"].map(prom_mid)).abs()
    pairs = pairs.merge(promoters[["promoter_id", "genes"]], on="promoter_id")
    pairs = pairs[["enhancer_id", "promoter_id", "genes", "n_loops", "distance", "loop_ids"]]
    pairs = pairs.sort_values(["enhancer_id", "promoter_id"]).reset_index(drop=True)

    linked = enhancers[enhancers["enhancer_id"].isin(pairs["enhancer_id"])]
    homer_bed = sort_bed(linked[["chrom", "start", "end", "enhancer_id"]]).assign(score=".", strand="+")

    write_table(pairs, args.pairs_out)
    write_bed(homer_bed, args.enhancers_out)
    n_loops_used = pairs["loop_ids"].str.split(",").explode().nunique() if len(pairs) else 0
    print(f"{n_loops} loops; {len(pairs)} enhancer-promoter pairs from {n_loops_used} loops; "
          f"{len(homer_bed)} loop-supported enhancers, "
          f"{pairs['promoter_id'].nunique()} promoters -> {args.pairs_out}")


if __name__ == "__main__":
    main()
