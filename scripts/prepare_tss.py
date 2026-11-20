#!/usr/bin/env python3
"""Write one transcription start site (TSS) per gene from a GENCODE GTF.

The TSS is the gene's first base on its own strand: `start` for + strand genes,
`end` for - strand genes. Output is BED6: chrom, start, end, gene_name, gene_id, strand.

Example:
  python scripts/prepare_tss.py \
      --gtf data/raw/reference/gencode_v29_annotation.gtf.gz \
      --gene-types protein_coding --chromosomes chr1 chr2 ... \
      --out data/processed/annotation/tss.bed
"""

import argparse
import re

import pandas as pd

from bed_utils import sort_bed, write_bed

GTF_COLUMNS = ["chrom", "source", "feature", "start", "end",
               "score", "strand", "frame", "attributes"]


def attribute(attributes, key):
    """Extract one quoted attribute value (e.g. gene_name) from the GTF attribute column."""
    return attributes.str.extract(rf'{key} "([^"]+)"', expand=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--gene-types", nargs="+", required=True)
    parser.add_argument("--chromosomes", nargs="+", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    gtf = pd.read_csv(args.gtf, sep="\t", comment="#", header=None,
                      names=GTF_COLUMNS, dtype={"chrom": str}, low_memory=False)
    genes = gtf[(gtf["feature"] == "gene") & gtf["chrom"].isin(args.chromosomes)].copy()
    genes["gene_type"] = attribute(genes["attributes"], "gene_type")
    genes = genes[genes["gene_type"].isin(args.gene_types)]
    genes["gene_id"] = attribute(genes["attributes"], "gene_id")
    genes["gene_name"] = attribute(genes["attributes"], "gene_name")

    # GTF is 1-based and inclusive; BED is 0-based, half-open.
    plus = genes["strand"] == "+"
    tss = pd.DataFrame({
        "chrom": genes["chrom"],
        "start": genes["start"].where(plus, genes["end"]) - 1,
        "gene_name": genes["gene_name"],
        "gene_id": genes["gene_id"],
        "strand": genes["strand"],
    })
    tss["end"] = tss["start"] + 1
    tss = sort_bed(tss[["chrom", "start", "end", "gene_name", "gene_id", "strand"]])

    write_bed(tss, args.out)
    print(f"{len(tss)} TSSs from {tss['gene_name'].nunique()} genes -> {args.out}")


if __name__ == "__main__":
    main()
