---
doc_id: workflow-020
domain: transcriptomics
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  RNA-seq differential expression analysis using featureCounts and DESeq2.
  Quantifies mapped reads per gene from BAM files, performs statistical
  testing with negative binomial modeling, and generates MA, Volcano, and
  PCA plots.
tags: [rna-seq, deseq2, featurecounts, differential-expression, transcriptomics, bioconductor]
---

# Workflow 020: RNA-Seq Differential Expression Analysis

RNA-seq differential expression analysis using [featureCounts](https://subread.sourceforge.net/) (Rsubread) and [DESeq2](https://bioconductor.org/packages/DESeq2/). This workflow takes aligned BAM files and a GTF gene annotation, quantifies reads per gene, identifies differentially expressed genes using a negative binomial generalized linear model, and produces MA, Volcano, and PCA plots.

## Overview

The pipeline implements a standard RNA-seq differential expression workflow. featureCounts assigns aligned reads to genomic features (genes) defined in a GTF annotation file, producing a count matrix. DESeq2 then normalizes counts using the median-of-ratios method (size factors), models gene-level counts with a negative binomial distribution, and tests for differential expression between two conditions. Log2 fold changes are shrunk using the adaptive shrinkage (ashr) estimator to reduce noise from low-count genes. Results are visualized with MA, Volcano, and PCA plots (Liao et al., 2014; Love et al., 2014).

## When to use this workflow

Use this workflow when you have aligned RNA-seq reads (BAM files) and want to identify genes that are differentially expressed between two conditions (e.g., treatment vs control). The workflow expects pre-aligned data — it does not perform read alignment.

For metabolomics analysis, use workflow-008. For protein function prediction, use workflow-015 (mDeepFRI). For protein structure prediction, use workflow-012 (Boltz-2).

## Architecture and data flow

```text
[00: Download Inputs] ──> [01: featureCounts] ──> [02: DESeq2] ──> [03: Plots]
         |                        |                      |                |
    *.bam, *.gtf             counts.csv            merged.csv       ma_plot.png
                                                   dds.rds          volcano_plot.png
                                                   pca_data.csv     pca_plot.png
```

Nodes run sequentially: 00 → 01 → 02 → 03.

## Input requirements

- **BAM files:** Aligned RNA-seq reads. One BAM file per sample. The default test data includes 4 BAM files: `Ctrl_1.bam`, `Ctrl_2.bam`, `Treat_1.bam`, `Treat_2.bam`.
- **GTF annotation:** Gene annotation file matching the reference genome used for alignment. The default test data uses `Saccharomyces_cerevisiae.gtf`.
- **Sample data:** Node 00 downloads test data from Zenodo record [18301020](https://zenodo.org/records/18301020) (S. cerevisiae). To use a different dataset, set the `zenodo_record_id` parameter and update `conditions` in Node 02.

### Condition assignment

Sample conditions are specified as a comma-separated string in Node 02, with one label per BAM file in alphabetical filename order:

```json
{ "conditions": "normal,normal,treatment,treatment" }
```

## Workflow nodes

### Node 00: Download Input Files

**Goal:** Fetch BAM and GTF files for analysis.

**Process:** Downloads all files from the specified Zenodo record via the Zenodo API. Additionally downloads the GTF annotation from the collab-workflows GitHub repository.

**Outputs:**
- `*.bam` — aligned read files
- `*.gtf` — gene annotation file

### Node 01: featureCounts

**Goal:** Quantify mapped reads per gene from BAM files.

**Process:** Runs `Rsubread::featureCounts()` with the configured parameters: strandedness, paired-end mode, feature type, attribute type, multi-overlap handling, and fractional counting. Counts reads assigned to each gene (grouped by `attribute_type`, default `gene_id`) across all BAM files. Only primary alignments are counted (`primaryOnly = TRUE`). Supports long-read mode via the `long_read` parameter.

**Scientific notes:** featureCounts assigns reads to genomic features by checking overlap between read coordinates and feature intervals defined in the GTF. When `allowMultiOverlap = TRUE` and `fraction = TRUE`, reads overlapping multiple features are counted fractionally (1/n for n overlapping features) rather than being discarded or double-counted. The `useMetaFeatures = TRUE` setting groups features (e.g., exons) by their parent attribute (e.g., gene_id) for gene-level quantification (Liao et al., 2014).

**Outputs:**
- `counts.csv` — combined annotation and count matrix
- `counts_table.csv` — raw count matrix only
- `counts_annotation.csv` — feature annotation

### Node 02: DESeq2 Analysis

**Goal:** Perform statistical differential expression analysis.

**Process:** Loads the count matrix from Node 01, auto-detects count columns (numeric, non-negative, excluding annotation columns), builds the sample metadata from the `conditions` parameter, and creates a DESeqDataSet with design `~ condition`. Runs the DESeq2 pipeline (estimation of size factors, dispersion estimation, negative binomial GLM fitting, Wald test). Extracts results for the treatment vs control contrast at the configured `alpha` threshold. Performs regularized log (rlog) transformation for PCA, extracts PC1/PC2 coordinates and variance percentages.

**Scientific notes:** DESeq2 normalizes for sequencing depth differences using the median-of-ratios method rather than simple total-count normalization, which is more robust to highly expressed genes. The negative binomial model accounts for both biological and technical variability (overdispersion). The rlog transformation stabilizes variance across the mean for visualization and clustering; for large datasets (>30 samples), `vst()` would be more computationally efficient (Love et al., 2014).

**Outputs:**
- `merged.csv` — DESeq2 results (baseMean, log2FoldChange, padj) merged with featureCounts annotation
- `dds.rds` — serialized DESeqDataSet R object for downstream analysis
- `res.rds` — serialized DESeq2 results object
- `pca_data.csv` — PC1/PC2 coordinates and condition labels

### Node 03: Visualization Plots

**Goal:** Generate publication-quality MA, Volcano, and PCA plots.

**Process:** Loads the DESeq2 objects and merged results. Performs log2 fold change shrinkage using `lfcShrink(type = "ashr")` for the MA plot. Classifies genes as Up, Down, or None based on adjusted p-value < `alpha`. Generates three ggplot2 figures with configurable colors and dimensions:
- **MA plot:** log10(mean expression) vs shrunken log2 fold change, with top 20 significant genes labeled
- **Volcano plot:** shrunken log2 fold change vs -log10(adjusted p-value), with significance threshold line
- **PCA plot:** PC1 vs PC2 colored by condition

Uses `ggrastr` for rasterized point rendering when available (improves performance with many genes). Gene labels use `ggrepel` for non-overlapping placement.

**Scientific notes:** The ashr (Adaptive Shrinkage) estimator shrinks noisy log2 fold change estimates from low-count genes toward zero, producing more reliable effect size estimates for visualization and ranking. This is an alternative to the `apeglm` method (DESeq2's primary recommendation); ashr is particularly useful when specifying contrasts not directly available to apeglm (Stephens, 2017).

**Outputs:**
- `ma_plot.png` — MA plot (mean expression vs fold change)
- `volcano_plot.png` — Volcano plot (fold change vs significance)
- `pca_plot.png` — PCA plot (sample clustering)
- `merged_with_regulation.csv` — results with Up/Down/None regulation flags

## Parameters

### zenodo_record_id

- **Type:** string
- **Default:** `"18301020"`
- **Node:** 00
- **Description:** Zenodo record ID to download BAM and GTF files from. The default record contains S. cerevisiae test data.

### threads

- **Type:** integer
- **Default:** `8`
- **Node:** 01
- **Description:** Number of CPU threads for featureCounts read counting.

### strandedness

- **Type:** integer
- **Default:** `1`
- **Node:** 01
- **Description:** Library strandedness for read counting.

| Value | Description |
|-------|-------------|
| `0` | Unstranded |
| `1` (default) | Stranded (sense) |
| `2` | Reverse stranded (antisense) |

### paired_end

- **Type:** boolean
- **Default:** `true`
- **Node:** 01
- **Description:** Set to `true` for paired-end sequencing data, `false` for single-end.

### feature_type

- **Type:** string
- **Default:** `"exon"`
- **Node:** 01
- **Description:** GTF feature type to count. Must match the third column of the GTF file (e.g., `exon`, `gene`, `CDS`).

### attribute_type

- **Type:** string
- **Default:** `"gene_id"`
- **Node:** 01
- **Description:** GTF attribute to group features by for gene-level quantification (e.g., `gene_id`, `locus_tag`).

### long_read

- **Type:** boolean
- **Default:** `false`
- **Node:** 01
- **Description:** Enable long-read mode (adds `-L` flag). Set to `true` for Oxford Nanopore or PacBio reads.

### allow_multi_overlap

- **Type:** boolean
- **Default:** `true`
- **Node:** 01
- **Description:** Allow reads to be assigned to multiple overlapping features. When combined with `use_fraction`, reads are counted fractionally.

### use_fraction

- **Type:** boolean
- **Default:** `true`
- **Node:** 01
- **Description:** Use fractional counting (1/n) when a read overlaps n features. Only applies when `allow_multi_overlap` is enabled.

### conditions

- **Type:** string
- **Default:** `"normal,normal,treatment,treatment"`
- **Node:** 02
- **Description:** Comma-separated condition labels for each sample, in alphabetical BAM filename order. Must have exactly as many labels as there are BAM files.

### control_condition

- **Type:** string
- **Default:** `"normal"`
- **Node:** 02
- **Description:** Label for the reference/control group in DESeq2 contrast.

### treatment_condition

- **Type:** string
- **Default:** `"treatment"`
- **Node:** 02
- **Description:** Label for the treatment/experimental group in DESeq2 contrast.

### alpha

- **Type:** float
- **Default:** `0.05`
- **Node:** 02, 03
- **Description:** Adjusted p-value threshold for calling differential expression. Genes with padj < alpha are flagged as significant.

### Plot parameters (Node 03)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `plot_width` | `10` | Plot width in inches |
| `plot_height` | `8` | Plot height in inches |
| `plot_dpi` | `300` | Plot resolution in DPI |
| `up_color` | `"firebrick"` | Color for upregulated genes |
| `down_color` | `"dodgerblue4"` | Color for downregulated genes |
| `none_color` | `"honeydew4"` | Color for non-significant genes |
| `highlight_color` | `"gold1"` | Color for labeled gene highlights |

## Outputs and interpretation

### merged_with_regulation.csv

The primary results file. Key columns:

| Column | Description |
|--------|-------------|
| `baseMean` | Mean normalized count across all samples |
| `log2FoldChange` | Estimated log2 fold change (treatment vs control) |
| `log2FoldChange.shrink` | Shrunken log2 fold change (ashr) |
| `padj` | Benjamini-Hochberg adjusted p-value |
| `Regulation` | `Up`, `Down`, or `None` based on padj < alpha |

### MA plot

Displays the relationship between mean expression level and fold change. Points spread along the y-axis at low expression indicate noisy estimates; shrunken fold changes reduce this noise. Significant genes are colored; top 20 are labeled.

### Volcano plot

Highlights genes that are both statistically significant (high on y-axis) and biologically meaningful (far from center on x-axis). The horizontal dashed line marks the alpha threshold. Genes in the upper corners are the strongest candidates for follow-up.

### PCA plot

Shows sample clustering based on overall gene expression. Samples should cluster by condition. Outlier samples or unexpected clustering may indicate batch effects or sample quality issues.

## Quick start

### Running with Docker

```bash
docker pull python:3.11-slim                                    # Node 00
docker pull chiral.sakuracr.jp/bioconductor:rna_seq_r_2025_12_27_v1  # Nodes 01–03
```

### Running on Silva

1. Select "RNA-seq Differential Expression Analysis" from the workflow list
2. Set `zenodo_record_id` to your data record (or use the default S. cerevisiae test)
3. Update `conditions` to match your experimental design
4. Set `strandedness` and `paired_end` to match your library prep
5. Click Run

### Test run

The default settings run the S. cerevisiae test dataset (4 BAM files, 2 conditions x 2 replicates). A successful run produces count matrices, differential expression results, and three PNG plots.

## References

- Liao, Y., Smyth, G. K. & Shi, W. "featureCounts: an efficient general purpose program for assigning sequence reads to genomic features." *Bioinformatics* 30(7):923–930, 2014. DOI: https://doi.org/10.1093/bioinformatics/btt656
- Love, M. I., Huber, W. & Anders, S. "Moderated estimation of fold changes and dispersion for RNA-seq data with DESeq2." *Genome Biology* 15:550, 2014. DOI: https://doi.org/10.1186/s13059-014-0550-8
- Stephens, M. "False discovery rates: a new deal." *Biostatistics* 18(2):275–294, 2017. DOI: https://doi.org/10.1093/biostatistics/kxw041
- [DESeq2 Bioconductor vignette](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [Rsubread/featureCounts](https://subread.sourceforge.net/)
