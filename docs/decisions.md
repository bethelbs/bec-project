# Decision log

Analysis and project decisions for this repository, newest at the bottom. Each entry gives the decision and the reason for it.

## 2026-11-19 — Scope of the analysis

This analysis covers an element-first catalogue of active promoters and enhancers in brain microvascular endothelial cells (BMEC), linked to target genes through Hi-C loops, ending in a table of candidate transcription factors from motif enrichment. A focused scope keeps every step checkable and rerunnable; further analyses can build on these outputs.

## 2026-11-19 — Element-first design

Regulatory elements are identified first, and candidate genes and transcription factors are derived from them. Starting from known marker genes would bias the analysis toward what is already known about BMEC identity, and could miss elements and regulators that have not been described before.

## 2026-11-19 — ENCODE processed files instead of raw reads

All inputs are processed ENCODE files on GRCh38 (peak calls, Hi-C loops, domains and compartments), not FASTQ files or the full Hi-C contact matrix (ENCFF068LTZ, several GB). These files come from ENCODE's uniform processing pipelines, so they are consistent with each other, and they keep the inputs small enough to download and rerun the whole analysis quickly on a standard machine.

## 2026-11-19 — HUVEC H3K4me1 as a proxy

No BMEC H3K4me1 ChIP-seq is available on ENCODE, so H3K4me1 from HUVEC (ENCSR000AKL) is used as a proxy enhancer mark. HUVEC is also an endothelial cell type, so it is the closest available match, but it is not brain-derived. The proxy is labelled in `metadata/samples.tsv`, the number of enhancers that depend on it is reported, and the caveat is stated in the report.

## 2026-11-19 — Lightweight motif step

Motif enrichment uses known motifs only, run on loop-supported enhancers (a few hundred regions) rather than on every enhancer. HOMER's own hg38 genome package is not installed; instead, a GRCh38 FASTA is downloaded by the download script, checksummed like the other raw files, and passed to the motif tool as a file path. This keeps the run time and the software environment small, and places the genome sequence under the same provenance checks as the other inputs.

## 2026-11-19 — No chromosome subset for CI (for now)

The workflow runs on the whole genome in every environment. With processed files only, the intersection steps take seconds and the motif step is limited to a few hundred regions, so a reduced test run is not needed. This will be revisited if a full run on the CI runner becomes slow.

## 2026-11-19 — Snakemake as the workflow manager

The pipeline is written in Snakemake (`workflow/Snakefile`, settings in `config/config.yaml`) and launched through Pixi. The steps are mostly Python scripts and interval operations on BED files, which map directly onto Snakemake's file-based rules. Snakemake also reruns only the steps affected by a change to an input or a setting, and needs nothing beyond Python.
