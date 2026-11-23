# Regulatory elements of human brain microvascular endothelial cells

[![Reproduce](https://github.com/bethelbs/bec-project/actions/workflows/reproduce.yml/badge.svg)](https://github.com/bethelbs/bec-project/actions/workflows/reproduce.yml)

## What is this project about?

This project asks which regulatory elements define the identity of human brain microvascular endothelial cells (BMEC), and which genes and transcription factors (TFs) those elements connect to. It takes an element-first approach: active promoters and enhancers are catalogued from histone-mark ChIP-seq, enhancers are linked to promoters through Hi-C chromatin loops, and the loop-supported enhancers are searched for enriched TF motifs. Known blood-brain-barrier marker genes are used as a check.

There is no publication associated with this project. The analysis builds on public ENCODE data; see [Data](#data) for the sources and `reports/bec-report.qmd` for the methods and references.

## What is in this repository?

| Folder | Contents |
| --- | --- |
| `config/` | `config.yaml`: every setting the analysis uses |
| `data/raw/` | Downloaded input data (not in Git; read-only) |
| `data/processed/` | Intermediate files made by the workflow (not in Git) |
| `docs/` | `decisions.md`: the analysis decisions and why they were made |
| `metadata/` | Data sources, sample sheet, checksums and data dictionary |
| `reports/` | `bec-report.qmd`: the report, computed from the results |
| `results/` | Tables, motif results and logs made by the workflow (not in Git) |
| `scripts/` | The download script and one script per analysis step |
| `workflow/` | `Snakefile`: the Snakemake workflow that runs the steps |

## What hardware does it need?

- **Operating system:** Linux or macOS 11 or newer (Apple silicon or Intel). Windows is not supported, because bedtools and HOMER have no Windows builds.
- **CPU and time:** the analysis takes about 3 minutes on 2 cores, after the one-time download (0.9 GB) and environment install (a few minutes each).
- **Memory:** about 2.5 GB at peak (HOMER); 4 GB of RAM is enough.
- **Disk:** about 8 GB free: 3 GB for the software environments, 0.9 GB for the raw data, and 3 GB of temporary space while the genome is unpacked for motif analysis.

## How do I reproduce the analyses?

You need [Git](https://git-scm.com) and [Pixi](https://pixi.sh). Then:

```bash
git clone https://github.com/bethelbs/bec-project.git
cd bec-project
pixi install --locked
pixi run reproduce
```

`pixi run reproduce` downloads and verifies the data (`download`), runs the Snakemake workflow (`run-snakemake`), and renders the report (`report`). Each step can also be run on its own with `pixi run <task>`; `pixi task list` shows them all. Every setting is in `config/config.yaml`; after changing one, `pixi run reproduce` reruns only the steps it affects.

## What outputs should I expect?

The report is written to `reports/bec-report.html`. The tables it is built from are in `results/tables/`:

| File | Contents |
| --- | --- |
| `summary.tsv` | Headline counts quoted in the report |
| `peak_qc.tsv` | Peak counts, widths and TSS proximity per mark |
| `ep_pairs.tsv` | Loop-supported enhancer–promoter pairs |
| `element_context.tsv` | Contact domain and A/B compartment of each element |
| `candidate_tfs.tsv` | Enriched known motifs (candidate TFs) |
| `target_genes.tsv` | Genes linked to loop-supported enhancers |

With the default settings you should get these values in `results/tables/summary.tsv`:

| Metric | Value |
| --- | --- |
| Promoters (`n_promoters`) | 17,104 |
| Enhancers (`n_enhancers`) | 43,278 |
| Enhancer–promoter pairs (`n_ep_pairs`) | 840 |
| Loop-supported enhancers (`n_loop_enhancers`) | 540 |
| Target genes (`n_target_genes`) | 290 |
| Candidate TF motifs (`n_candidate_tfs`) | 33 |

The motif step uses a fixed random seed, so these values are identical on every run.

## How can I contact the authors or discuss the project?

Open an [issue](https://github.com/bethelbs/bec-project/issues) on GitHub, or email Bethel Bahiru Sisay at <bethel.sisay@mbzuai.ac.ae>.

## Data

All input data come from the [ENCODE portal](https://www.encodeproject.org) (GRCh38): ChIP-seq peaks for H3K27ac (ENCSR111QCU), H3K4me3 (ENCSR000DTC), H3K27me3 (ENCSR032CAS), H3K36me3 (ENCSR587VHH) and CTCF (ENCSR000DTA) in BMEC; H3K4me1 in HUVEC (ENCSR000AKL), used as a proxy because no BMEC H3K4me1 data exist; Hi-C loops, domains and compartments in BMEC (ENCSR507AHE); and the GRCh38 genome, GENCODE v29 annotation and ENCODE blacklist. The data are not stored in this repository. `scripts/download_data.py` downloads them, checks each file against ENCODE's MD5 and against the SHA-256 checksums in `metadata/checksums.sha256`, and makes them read-only. `metadata/README.md` describes every file.

## Licence

- **Code** (scripts, workflow, configuration): MIT licence; see `LICENSE`.
- **Report text and figures:** [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).
- **Data:** the ENCODE data keep ENCODE's own terms of use and are not relicensed here; cite ENCODE when using them.

## Citation

If you use this work, please cite it using the metadata in `CITATION.cff` (GitHub shows it under "Cite this repository").
