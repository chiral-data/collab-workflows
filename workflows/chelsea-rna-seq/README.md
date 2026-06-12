---
doc_id: chelsea-rna-seq
domain: transcriptomics
doc_type: workflow
version: "0.1.0"
deprecated: true
description: >
  Legacy RNA-seq input data directory. Contains S. cerevisiae BAM files and GTF
  annotation. Superseded by workflow-020.
tags: [rna-seq, legacy, deprecated]
---

# chelsea-rna-seq (Deprecated)

This directory contains input data files for RNA-seq analysis but is **not a functioning workflow** — it has no `workflow.toml`, no node directories, and no processing scripts.

## Current status

The directory holds Saccharomyces cerevisiae test data (4 BAM files + GTF annotation) and a download script that fetches the same files from Zenodo record [18301020](https://zenodo.org/records/18301020).

## Use workflow-020 instead

For RNA-seq differential expression analysis, use **workflow-020** which implements a complete 4-node pipeline (download, featureCounts, DESeq2, visualization) using the same S. cerevisiae test dataset.

## Files

- `input_files/Ctrl_1.bam`, `Ctrl_2.bam` — control condition BAM files
- `input_files/Treat_1.bam`, `Treat_2.bam` — treatment condition BAM files
- `input_files/Saccharomyces_cerevisiae.gtf` — gene annotation
- `input_files/download_data.sh` — shell script to re-download files from Zenodo
