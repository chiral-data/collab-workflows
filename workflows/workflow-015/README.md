# Workflow 015: mDeepFRI Protein Function Prediction

AI-based protein function annotation using [Metagenomic-DeepFRI](https://github.com/bioinf-mcb/Metagenomic-DeepFRI).

## Overview

This workflow takes a FASTA file of protein sequences and annotates each with Gene Ontology (GO) terms using mDeepFRI's deep learning pipeline. For proteins with structural matches in PDB100, a Graph Convolutional Network (GCN) makes structure-aware predictions; sequences without hits fall back to a CNN sequence-only model. It produces a GO term predictions TSV and a self-contained HTML report.

## Nodes

| Node | Name | Description |
|------|------|-------------|
| 01 | Validate Inputs | Validate FASTA format, confirm protein sequences, filter by length |
| 02 | Predict Protein Function | Run mDeepFRI (MMseqs2 alignment + GCN/CNN inference) to assign GO terms |
| 03 | Generate Visualization | Build a self-contained HTML report with GO term tables and score distributions |

Nodes run sequentially: 01 → 02 → 03.

## Parameters

| Parameter | Node | Type | Default | Description |
|-----------|------|------|---------|-------------|
| `min_length` | 01 | integer | 10 | Minimum sequence length (aa) to include |
| `max_length` | 01 | integer | 5000 | Maximum sequence length (aa) to include |
| `prediction_modes` | 02 | string | `mf bp cc` | Space-separated prediction modes: `mf` (Molecular Function), `bp` (Biological Process), `cc` (Cellular Component) |
| `threads` | 02 | integer | 1 | CPU threads for MMseqs2 search |
| `mmseqs_sensitivity` | 02 | float | 5.7 | MMseqs2 sensitivity (1.0–7.5; higher = more sensitive) |
| `min_score` | 03 | float | 0.3 | Minimum confidence score threshold for GO terms shown in report |
| `top_n_terms` | 03 | integer | 10 | Maximum GO terms to display per protein per mode |

## Input Format

A FASTA file with one or more protein sequences. Swiss-Prot/TrEMBL headers (`sp|ACC|NAME` or `tr|ACC|NAME`) are normalized to the accession ID automatically.

```
>sp|P62988|UBIQ_HUMAN Ubiquitin OS=Homo sapiens
MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG
```

## Outputs

- **`validated.fasta`** — Filtered, normalized FASTA passed to prediction
- **`validation_report.json`** — Counts of sequences filtered and reasons
- **`results.tsv`** — Per-protein GO term predictions with confidence scores and network type (GCN/CNN)
- **`alignment_summary.tsv`** — Per-protein structural alignment statistics against PDB100
- **`report.html`** — Self-contained HTML dashboard with GO term tables and score distribution chart

## Containers

| Node | Image |
|------|-------|
| 01, 02, 03 | `mdeepfri:latest` |

The image installs `mdeepfri` via pip (Python 3.12), patches known upstream bugs, and downloads model weights (~1 GB) and the PDB100 MMseqs2 index at build time.

## Example

The test input (`input_files/sample_proteins.fasta`) contains 10 well-characterized proteins (ubiquitin, thioredoxin, GFP, albumin, etc.). A typical run aligns ~70% of sequences against PDB100 for GCN prediction and falls back to CNN for the remainder.
