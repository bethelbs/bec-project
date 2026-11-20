#!/usr/bin/env python3
"""Add 3D-genome context to promoters and enhancers: contact domain and A/B compartment.

Domain: the smallest Hi-C contact domain that fully contains the element
  (domains can be nested); empty if none does.

Compartment: the mean compartment-track value over the element. The sign of a
  compartment eigenvector is arbitrary per chromosome, so each chromosome is
  oriented with the data: if most H3K27ac-marked promoters on a chromosome have
  negative values, that chromosome's sign is flipped. After orientation,
  positive = A (active), negative = B (inactive).

Output (TSV): element_id, type, chrom, start, end, domain_id, compartment_value, compartment

Example:
  python scripts/annotate_context.py \
      --promoters data/processed/elements/promoters.tsv \
      --enhancers data/processed/elements/enhancers.tsv \
      --domains data/raw/hic/bmec_hic_domains.bedpe.gz \
      --compartments data/raw/hic/bmec_hic_compartments.bigWig \
      --out results/tables/element_context.tsv
"""

import argparse

import numpy as np
import pandas as pd
import pyBigWig

from bed_utils import bedtools_df, read_bed, temp_bed, write_table


def assign_domains(elements, domains_path):
    """Smallest domain fully containing each element (domain_id or empty)."""
    domains = read_bed(domains_path).iloc[:, :3]
    domains.columns = ["chrom", "start", "end"]
    domains = domains.sort_values(["chrom", "start", "end"]).reset_index(drop=True)
    domains["domain_id"] = [f"D{i:04d}" for i in range(1, len(domains) + 1)]
    with temp_bed(elements[["chrom", "start", "end", "element_id"]]) as a, temp_bed(domains) as b:
        hits = bedtools_df("intersect", "-f", "1.0", "-wa", "-wb", "-a", a, "-b", b)
    if hits.empty:
        return pd.Series("", index=elements["element_id"], name="domain_id")
    hits["size"] = hits[6] - hits[5]
    best = hits.sort_values("size").drop_duplicates(3).set_index(3)[7]
    return best.reindex(elements["element_id"]).fillna("").rename("domain_id")


def compartment_values(elements, bigwig_path):
    """Mean compartment-track value over each element (NaN where the track has no data)."""
    bw = pyBigWig.open(bigwig_path)
    chroms = bw.chroms()
    values = []
    for chrom, start, end in elements[["chrom", "start", "end"]].itertuples(index=False):
        v = bw.stats(chrom, int(start), int(end), type="mean")[0] if chrom in chroms else None
        values.append(np.nan if v is None else v)
    bw.close()
    return np.array(values, dtype=float)


def orient(elements):
    """Flip the compartment sign per chromosome so active promoters sit on the positive side."""
    flipped = []
    active = elements[(elements["type"] == "promoter") & elements["h3k27ac"].fillna(False).astype(bool)]
    for chrom, sub in active.groupby("chrom"):
        vals = sub["compartment_value"].dropna()
        if len(vals) >= 20 and (vals > 0).mean() < 0.5:
            elements.loc[elements["chrom"] == chrom, "compartment_value"] *= -1
            flipped.append(chrom)
    return flipped


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--promoters", required=True)
    parser.add_argument("--enhancers", required=True)
    parser.add_argument("--domains", required=True)
    parser.add_argument("--compartments", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    prom = pd.read_csv(args.promoters, sep="\t").rename(columns={"promoter_id": "element_id"})
    enh = pd.read_csv(args.enhancers, sep="\t").rename(columns={"enhancer_id": "element_id"})
    elements = pd.concat([prom.assign(type="promoter"), enh.assign(type="enhancer")], ignore_index=True)

    elements["domain_id"] = assign_domains(elements, args.domains).values
    elements["compartment_value"] = compartment_values(elements, args.compartments)
    flipped = orient(elements)
    elements["compartment"] = np.select(
        [elements["compartment_value"] > 0, elements["compartment_value"] < 0], ["A", "B"], default="")

    out = elements[["element_id", "type", "chrom", "start", "end",
                    "domain_id", "compartment_value", "compartment"]]
    write_table(out, args.out)
    in_domain = (out["domain_id"] != "").mean()
    frac_a = out.groupby("type")["compartment"].apply(lambda c: (c == "A").mean()).round(3).to_dict()
    print(f"{len(out)} elements; {in_domain:.1%} inside a domain; fraction in A: {frac_a}; "
          f"sign flipped on: {', '.join(flipped) or 'none'} -> {args.out}")


if __name__ == "__main__":
    main()
