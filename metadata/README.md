# metadata

Metadata

This folder describes every input file used by the analysis: where it comes from, what it contains and how it is checked. These files are small, written by hand and committed to Git. The data files themselves are not in Git; scripts/download_data.py downloads them into data/raw/ and makes them read-only.

All inputs are processed files from the ENCODE portal, on the GRCh38 assembly. Together they are about 0.92 GB, of which 0.87 GB is the reference genome sequence.

Files in this folder
File	Contents
sources.tsv	One row per downloaded file: where to get it, its provider checksum and where it goes. Read by scripts/download_data.py.
samples.tsv	One row per experimental data file (ChIP-seq and Hi-C): what it measures and which ENCODE experiment it comes from. Read by the workflow.
checksums.sha256	SHA-256 checksums of every file in data/raw/, in sha256sum format (hash, two spaces, relative path). Created after the first download and used to verify the stored files.
README.md	This data dictionary.
sources.tsv

Tab-separated, one row per file.

Column	Meaning
name	Short identifier for the file, used throughout the project.
url	Direct download link on the ENCODE portal.
md5	MD5 checksum published by ENCODE for the file. Used only to check that the download is complete and correct.
unpack	none: the file is kept exactly as downloaded (compressed files stay compressed).
dest	Path of the file in the project, relative to the project root.
samples.tsv

Tab-separated, one row per ChIP-seq or Hi-C file. Reference files (genome, annotation, blacklist, chromosome sizes) are listed only in sources.tsv.

Column	Meaning
name	Same identifier as in sources.tsv.
assay	ChIP-seq or Hi-C.
target	Histone mark or protein targeted by the ChIP-seq; none for Hi-C.
cell_type	ENCODE biosample term name.
encode_experiment	ENCODE experiment accession (ENCSR…).
encode_file	ENCODE file accession (ENCFF…) of the file used.
output_type	ENCODE output type of the file (see below).
biological_replicates	Biological replicates the file was derived from, comma-separated.
is_proxy	yes if the file comes from a different cell type and stands in for missing BMEC data; otherwise no.
path	Path of the file in the project, relative to the project root.
Input files
ChIP-seq peaks (data/raw/chip/)

All ChIP-seq files are narrowPeak BED files (gzip-compressed) from the ENCODE4 histone and TF ChIP-seq pipelines (v1.5.1). For each experiment, the file used is the one ENCODE marks as the preferred default of the experiment's current default analysis.

Name	Mark	Cell type	Experiment	File	Output type	Peaks
bmec_h3k27ac_peaks	H3K27ac	BMEC	ENCSR111QCU	ENCFF718IQW	pseudoreplicated peaks	62,171
bmec_h3k4me3_peaks	H3K4me3	BMEC	ENCSR000DTC	ENCFF430BFN	pseudoreplicated peaks	32,093
bmec_h3k27me3_peaks	H3K27me3	BMEC	ENCSR032CAS	ENCFF732VIJ	replicated peaks	30,024
bmec_h3k36me3_peaks	H3K36me3	BMEC	ENCSR587VHH	ENCFF784BVN	replicated peaks	100,404
bmec_ctcf_peaks	CTCF	BMEC	ENCSR000DTA	ENCFF985NAQ	IDR thresholded peaks	62,398
huvec_h3k4me1_peaks	H3K4me1	HUVEC (proxy)	ENCSR000AKL	ENCFF213BAF	pseudoreplicated peaks	112,482

NarrowPeak columns: chromosome, start (0-based), end, name, score, strand, signal value, −log10 p-value, −log10 q-value, summit offset from start.

ENCODE output types: replicated peaks are peaks found in both biological replicates; pseudoreplicated peaks are peaks found in both halves of the pooled reads; IDR thresholded peaks (TF ChIP-seq) pass the irreproducible discovery rate threshold across replicates.

Proxy: no H3K4me1 ChIP-seq exists for BMEC on ENCODE, so H3K4me1 from human umbilical vein endothelial cells (HUVEC) is used as a proxy enhancer mark. It is marked is_proxy = yes in samples.tsv.

Hi-C (data/raw/hic/)

From experiment ENCSR507AHE (BMEC, replicates 1 and 2 combined), ENCODE4 Hi-C pipeline (Juicer tools 2.13.06).

Name	File	Format	Contents
bmec_hic_loops	ENCFF915TBK	BEDPE (gzip)	1,662 chromatin loops. Columns 1–6 give the two anchors (chr, start, end); later columns give observed and expected counts and FDRs.
bmec_hic_domains	ENCFF698ULZ	BEDPE (gzip)	633 contact domains (TADs). Columns 1–3 and 4–6 give the same interval; later columns give domain scores.
bmec_hic_compartments	ENCFF245BDA	bigWig	A/B compartment signal (eigenvector); positive and negative values separate the two compartments.

Lines starting with # in the BEDPE files are headers.

Reference files (data/raw/reference/)
Name	ENCODE file	Contents
grch38_genome	GRCh38_no_alt_analysis_set_GCA_000001405.15	GRCh38 genome sequence (FASTA, gzip), the same build ENCODE aligned all the data above to. Used to extract enhancer sequences for motif analysis.
grch38_chrom_sizes	ENCFF667IGK	Chromosome names and lengths (two-column TSV).
grch38_blacklist	ENCFF356LFX	ENCODE unified GRCh38 exclusion list (BED, gzip): regions with artefactual signal, removed from peak sets.
gencode_v29_annotation	ENCFF159KBI	GENCODE v29 primary gene annotation (GTF, gzip), as used by ENCODE. Source of gene TSSs.
Source and licence

All files were produced by the ENCODE Consortium and are distributed through the ENCODE portal. They keep ENCODE's terms of use and are not relicensed by this project; when using them, cite ENCODE as described on the portal's citation page. Accessions and checksums above were taken from the ENCODE portal on 2026-11-20.
