#!/usr/bin/env python3
"""Make a reproducible random genomic background for motif enrichment.

The target regions are copied until there are at least --n regions, then each
copy is placed at a random position on the same chromosome with
`bedtools shuffle`, using a fixed --seed so every run gives the same background.
Blacklisted regions and the targets themselves are excluded. Output is BED6
(chrom, start, end, id, ".", "+"), ready for HOMER's -bg option.

HOMER would otherwise draw its own random background each run, so results
would differ slightly between runs; a fixed background makes them identical.

Example:
  python scripts/make_background.py \
      --targets data/processed/elements/loop_enhancers.bed \
      --chrom-sizes data/raw/reference/grch38_chrom_sizes.tsv \
      --blacklist data/raw/reference/grch38_blacklist.bed.gz \
      --chromosomes chr1 chr2 ... --n 10000 --seed 42 \
      --out data/processed/motifs/background.bed
"""

import argparse
import math

import pandas as pd

from bed_utils import bedtools_df, read_bed, sort_bed, temp_bed, write_bed


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--targets", required=True)
    parser.add_argument("--chrom-sizes", required=True)
    parser.add_argument("--blacklist", required=True)
    parser.add_argument("--chromosomes", nargs="+", required=True)
    parser.add_argument("--n", type=int, required=True, help="minimum number of background regions")
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    targets = read_bed(args.targets).iloc[:, :3]
    targets.columns = ["chrom", "start", "end"]
    copies = math.ceil(args.n / len(targets))
    template = pd.concat([targets] * copies, ignore_index=True)

    sizes = pd.read_csv(args.chrom_sizes, sep="\t", header=None, names=["chrom", "size"])
    sizes = sizes[sizes["chrom"].isin(args.chromosomes)]
    exclude = pd.concat([read_bed(args.blacklist).iloc[:, :3].set_axis(["chrom", "start", "end"], axis=1),
                         targets], ignore_index=True)

    with temp_bed(template) as t, temp_bed(sizes) as g, temp_bed(exclude) as x:
        bg = bedtools_df("shuffle", "-i", t, "-g", g, "-excl", x,
                         "-chrom", "-noOverlapping", "-seed", args.seed)
    bg = sort_bed(bg.iloc[:, :3].set_axis(["chrom", "start", "end"], axis=1))
    bg["id"] = [f"BG{i:06d}" for i in range(1, len(bg) + 1)]
    bg["score"], bg["strand"] = ".", "+"

    write_bed(bg, args.out)
    print(f"{len(bg)} background regions ({copies} per target, seed {args.seed}) -> {args.out}")


if __name__ == "__main__":
    main()
