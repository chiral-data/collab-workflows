# RNA-Seq Differential Expression Analysis

Quantify mapped reads and identify differentially expressed genes from RNA-Seq data.

## Pipeline

```
00_download_inputs  →  01_featurecounts  →  02_deseq2  →  03_plots
```

1. **Download inputs** — fetches BAM and GTF files from a Zenodo record
2. **featureCounts** — counts reads per gene from BAM files against a GTF annotation
3. **DESeq2** — statistical differential expression analysis (size-factor normalization, negative binomial model)
4. **Plots** — MA plot, Volcano plot, PCA plot

Nodes 01–03 use `chiral.sakuracr.jp/bioconductor:rna_seq_r_2025_12_27_v1`. Node 00 uses `python:3.11-slim`.

## Input Files

The `00_download_inputs` node downloads files automatically from Zenodo. The default record is [18301020](https://zenodo.org/records/18301020) (Saccharomyces cerevisiae test data):

- 4 BAM files: `Ctrl_1.bam`, `Ctrl_2.bam`, `Treat_1.bam`, `Treat_2.bam`
- 1 GTF annotation: `Saccharomyces_cerevisiae.gtf`

To use a different dataset, set `zenodo_record_id` in `00_download_inputs/params.json`:

```json
{ "zenodo_record_id": "YOUR_RECORD_ID" }
```

## Outputs

| File                          | Description                                      |
| ----------------------------- | ------------------------------------------------ |
| `counts.csv`                  | Raw gene counts per sample                       |
| `merged.csv`                  | DESeq2 results merged with annotation            |
| `merged_with_regulation.csv`  | Same as above with Up/Down/None regulation flags |
| `dds.rds` / `res.rds`         | DESeq2 R objects for downstream analysis         |
| `ma_plot.png`                 | MA plot (mean expression vs log fold change)     |
| `volcano_plot.png`            | Volcano plot (fold change vs significance)       |
| `pca_plot.png`                | PCA plot of sample clustering                    |

## Key Parameters

| Parameter             | Default   | Description                                         |
| --------------------- | --------- | --------------------------------------------------- |
| `strandedness`        | 1         | 0 = unstranded, 1 = stranded, 2 = reverse stranded |
| `paired_end`          | true      | Set to false for single-end reads                   |
| `feature_type`        | exon      | GTF feature type to count (must match GTF column 3) |
| `attribute_type`      | gene_id   | GTF attribute to group features by                  |
| `control_condition`   | normal    | Label for the reference/control group               |
| `treatment_condition` | treatment | Label for the treatment group                       |
| `alpha`               | 0.05      | Adjusted p-value significance threshold             |

### Sample conditions

`02_deseq2/params.json` sets the condition for each sample. List one label per BAM file in alphabetical filename order:

```json
{ "conditions": "normal,normal,treatment,treatment" }
```

Update this when using different samples or a different group layout.
