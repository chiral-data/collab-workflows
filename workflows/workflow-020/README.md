# RNA-Seq Differential Expression Analysis

Quantify mapped reads and identify differentially expressed genes from RNA-Seq data.

## Pipeline

```
01_featurecounts  →  02_deseq2  →  03_plots
```

1. **featureCounts** — counts reads per gene from BAM files against a GTF annotation
2. **DESeq2** — statistical differential expression analysis (size-factor normalization, negative binomial model)
3. **Plots** — MA plot, Volcano plot, PCA plot

## Input Files

- One or more sorted BAM files (aligned reads)
- A GTF annotation file matching the reference genome used for alignment

### Test data (Saccharomyces cerevisiae)

A test dataset is hosted on Zenodo (record [18301020](https://zenodo.org/records/18301020)):

- 4 BAM files: `Ctrl_1.bam`, `Ctrl_2.bam`, `Treat_1.bam`, `Treat_2.bam`
- 1 GTF annotation: `Saccharomyces_cerevisiae.gtf`

Download with:

```bash
mkdir -p input_files && cd input_files

curl -s "https://zenodo.org/api/records/18301020" \
  | jq -r '.files[].links.self' \
  | while read -r url; do
      filename=$(basename "$(dirname "$url")")
      echo "Downloading $filename..."
      curl -L -o "$filename" "$url"
    done
```

Requires `curl` and `jq`.

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
