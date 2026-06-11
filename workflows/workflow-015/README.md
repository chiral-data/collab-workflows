---
doc_id: workflow-015
domain: protein-function
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  AI-based protein function annotation using Metagenomic-DeepFRI. Takes protein
  sequences in FASTA format and assigns Gene Ontology terms using graph
  convolutional networks (structure-aware) or convolutional neural networks
  (sequence-only fallback).
tags: [mdeepfri, deepfri, protein-function, gene-ontology, gcn, cnn, mmseqs2]
---

# Workflow 015: mDeepFRI Protein Function Prediction

AI-based protein function annotation using [Metagenomic-DeepFRI](https://github.com/bioinf-mcb/Metagenomic-DeepFRI). This workflow takes protein sequences in FASTA format, searches for structural homologs in PDB100 via MMseqs2, and predicts Gene Ontology (GO) terms using either a Graph Convolutional Network (GCN) for proteins with structural matches or a Convolutional Neural Network (CNN) for sequence-only prediction. Results are presented in an interactive HTML dashboard.

## Overview

mDeepFRI is a high-throughput pipeline built on [DeepFRI](https://github.com/flatironinstitute/DeepFRI) (Deep Functional Residue Identification), optimized for metagenomic-scale datasets. The pipeline first uses MMseqs2 to search each query protein against PDB100 (a clustered PDB structure database compressed with FoldComp). Proteins with structural hits above a similarity threshold get contact maps derived from the matched structure, which feed into a GCN that makes structure-aware GO term predictions. Proteins without structural matches fall back to a CNN that predicts function from sequence alone. The mDeepFRI implementation achieves 2–12x speedup over standard DeepFRI through ONNX model optimization (Gligorijevic et al., 2021).

The workflow uses model weights v1.1, which are finetuned on AlphaFold-predicted structures and machine-generated GO annotations from UniProt, improving accuracy over the original DeepFRI publication weights (v1.0).

## When to use this workflow

Use this workflow when you have one or more protein sequences and want to predict their biological function in terms of Gene Ontology annotations — Molecular Function (MF), Biological Process (BP), and Cellular Component (CC). It is particularly suited for uncharacterized or metagenomic proteins where experimental annotation is unavailable.

Do not use this workflow for protein structure prediction — use workflow-012 (Boltz-2). For protein-ligand binding analysis, use workflow-004 (AutoDock Vina) or workflow-002/003 (Smina). For ADMET property prediction of small molecules, use workflow-014 (ADMET-AI).

## Architecture and data flow

```text
[00: Download] ──> [01: Validate Inputs] ──> [02: Predict] ──> [03: Visualize]
       |                    |                       |                  |
  sample_proteins     validated.fasta          results.tsv        report.html
  .fasta           validation_report.json   alignment_summary.tsv
```

Nodes run sequentially: 00 → 01 → 02 → 03.

## Input requirements

- **Protein sequences in FASTA format:** One or more protein sequences with standard amino acid characters. Files with `.fasta`, `.fa`, or `.faa` extensions are accepted.
- **Header format:** Swiss-Prot/TrEMBL headers (`sp|ACC|NAME` or `tr|ACC|NAME`) are automatically normalized to the accession ID for compatibility with MMseqs2.
- **Sample data:** `input_files/sample_proteins.fasta` contains 10 well-characterized proteins (ubiquitin, thioredoxin, GFP, hemoglobin, lysozyme, beta-lactamase, beta-2-microglobulin, and others).

```
>sp|P62988|UBIQ_HUMAN Ubiquitin OS=Homo sapiens
MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG
```

## Workflow nodes

### Node 00: Download Sample Inputs

**Goal:** Fetch the sample protein FASTA file for testing.

**Process:** Downloads `sample_proteins.fasta` from the collab-workflows GitHub repository using urllib.

**Outputs:**
- `sample_proteins.fasta` — 10 well-characterized protein sequences

### Node 01: Validate Inputs

**Goal:** Validate FASTA format, confirm sequences are protein (not nucleotide), and filter by length.

**Process:** Parses all `.fasta`/`.fa`/`.faa` files from the inputs directory. For each sequence: checks that it contains valid amino acid characters (standard 20 + IUPAC ambiguity codes), verifies it is not a nucleotide sequence (rejects sequences composed solely of A/C/G/T/U), applies length filters (`min_length` to `max_length`), deduplicates by header ID, and normalizes Swiss-Prot/TrEMBL headers to accession IDs. Writes passing sequences to `validated.fasta` in 80-character line format.

**Scientific notes:** The nucleotide detection heuristic flags sequences where the unique character set is a subset of {A, C, G, T, U} or where >90% of unique characters are nucleotide bases with 5 or fewer unique characters total. This prevents accidental submission of DNA/RNA sequences.

**Outputs:**
- `validated.fasta` — filtered, normalized FASTA
- `validation_report.json` — counts of total, valid, and filtered sequences with reasons (not_protein, invalid_chars, too_short, too_long, duplicates)

### Node 02: Predict Protein Function

**Goal:** Run mDeepFRI to assign GO terms to each validated protein sequence.

**Process:** Invokes `mDeepFRI predict-function` with the validated FASTA, model weights from `/opt/mdeepfri-weights`, and the configured prediction modes. The `--skip-matrix` flag is used to skip contact map matrix output. For each prediction mode (`-p mf`, `-p bp`, `-p cc`), mDeepFRI:
1. Searches the query against PDB100 via MMseqs2 at the configured sensitivity
2. For structural hits: derives a contact map from the matched PDB structure and runs the GCN model
3. For sequences without hits: runs the CNN sequence-only model
4. Outputs per-protein GO term predictions with confidence scores

**Scientific notes:** The GCN model takes a protein contact map (residue-residue distances from 3D structure) as a graph and uses graph convolutions to predict function, achieving higher accuracy than sequence-only methods. The CNN fallback uses 1D convolutions over the raw amino acid sequence. Model weights v1.1 are finetuned on AlphaFold structures and UniProt machine-generated GO annotations. MMseqs2 sensitivity of 5.7 (default) balances speed and structural hit rate; increasing to 7.5 finds more remote homologs at the cost of longer search time.

**Outputs:**
- `results.tsv` — per-protein GO term predictions with columns: protein, prediction_mode, go_term, go_name, score, network_type (GCN or CNN)
- `alignment_summary.tsv` — per-protein structural alignment statistics against PDB100

### Node 03: Generate Visualization

**Goal:** Build a self-contained HTML dashboard summarizing GO term predictions.

**Process:** Reads `results.tsv` and `alignment_summary.tsv`. Filters predictions by `min_score` threshold, sorts by score descending, and caps at `top_n_terms` per protein per mode. Generates an HTML report with Bootstrap 5 styling and Plotly.js charts including: summary cards (proteins annotated, total predictions, score threshold), a score distribution histogram, and per-protein expandable cards with GO term tables showing mode, GO ID, name, confidence score (color-coded), and model type (GCN/CNN badge).

**Scientific notes:** Score color coding in the report: green (score >= 0.7, high confidence), orange (0.4–0.7, moderate), red (< 0.4, low confidence). These thresholds are visual guides — the `min_score` parameter controls which predictions appear at all.

**Outputs:**
- `report.html` — self-contained HTML dashboard with GO term tables and score distribution chart

## Parameters

### min_length

- **Type:** integer
- **Default:** `10`
- **Node:** 01
- **Description:** Minimum protein sequence length (amino acids) to include after validation.

### max_length

- **Type:** integer
- **Default:** `5000`
- **Node:** 01
- **Description:** Maximum protein sequence length (amino acids) to include after validation.

### prediction_modes

- **Type:** string
- **Default:** `"mf bp cc"`
- **Node:** 02
- **Description:** Space-separated GO prediction modes to run.

| Mode | Description |
|------|-------------|
| `mf` | Molecular Function — what the protein does biochemically |
| `bp` | Biological Process — the larger biological pathway or process |
| `cc` | Cellular Component — where the protein localizes in the cell |

### threads

- **Type:** integer
- **Default:** `1`
- **Node:** 02
- **Description:** Number of CPU threads for the MMseqs2 structural search step. Increase for large input sets.

### mmseqs_sensitivity

- **Type:** float
- **Default:** `5.7`
- **Node:** 02
- **Description:** MMseqs2 search sensitivity (range 1.0–7.5). Higher values find more remote structural homologs but increase search time.

| Value | Use case |
|-------|----------|
| `1.0`–`4.0` | Fast search, close homologs only |
| `5.7` (default) | Balanced sensitivity and speed |
| `7.0`–`7.5` | Maximum sensitivity, slowest |

### min_score

- **Type:** float
- **Default:** `0.3`
- **Node:** 03
- **Description:** Minimum confidence score threshold for GO terms to include in the HTML report. Predictions below this threshold are excluded from the visualization.

### top_n_terms

- **Type:** integer
- **Default:** `10`
- **Node:** 03
- **Description:** Maximum number of GO terms to display per protein per prediction mode in the report.

## Outputs and interpretation

### Confidence scores

mDeepFRI outputs a confidence score (0–1) for each GO term prediction, reflecting the model's certainty that the term applies to the query protein. Higher scores indicate stronger predictions.

| Range | Interpretation |
|-------|---------------|
| >= 0.7 | High confidence — strong functional prediction |
| 0.4–0.7 | Moderate confidence — plausible annotation, consider validation |
| 0.3–0.4 | Low confidence — included by default threshold but uncertain |
| < 0.3 | Below default threshold — filtered out unless `min_score` is lowered |

### Network type (GCN vs CNN)

Each prediction is labeled with the model that produced it:
- **GCN** — structure-aware prediction using a contact map from a PDB100 structural match. Generally more accurate.
- **CNN** — sequence-only fallback when no structural homolog is found. Still informative but less precise for structure-dependent functions.

A typical run with diverse proteins aligns ~70% of sequences against PDB100 for GCN prediction, with the remainder falling back to CNN.

### results.tsv

Tab-separated file with columns: `protein`, `prediction_mode`, `go_term`, `go_name`, `score`, `network_type`. One row per GO term prediction.

### alignment_summary.tsv

Tab-separated file with per-protein structural alignment statistics from the MMseqs2 search against PDB100.

### report.html

Self-contained HTML dashboard with interactive score distribution histogram and per-protein GO term tables. Can be opened directly in a browser without a server.

## Quick start

### Running with Docker

```bash
docker build -t mdeepfri:latest .
```

The Docker build downloads model weights (~1 GB) and pre-builds the PDB100 MMseqs2 index, so the first build takes several minutes.

### Running on Silva

1. Select "mDeepFRI Protein Function Prediction" from the workflow list
2. Upload your protein FASTA file (or use the default sample_proteins.fasta)
3. Adjust `prediction_modes` if you only need specific GO categories
4. Click Run

### Test vs production settings

| Setting | Test (default) | Production |
|---------|---------------|------------|
| `prediction_modes` | `mf bp cc` | As needed |
| `threads` | `1` | `4`–`8` |
| `mmseqs_sensitivity` | `5.7` | `5.7`–`7.5` |
| `min_score` | `0.3` | `0.3`–`0.5` |

A successful test run with the 10-protein sample completes in a few minutes and produces GO term annotations for all three categories.

## Troubleshooting

### Low GCN hit rate

If most proteins fall back to CNN, the input sequences may be distant from known PDB structures. Try increasing `mmseqs_sensitivity` to 7.0–7.5 to find more remote homologs.

### Empty results

Ensure input sequences are protein (not nucleotide) and within the length range. Check `validation_report.json` for filter statistics.

## References

- Gligorijevic, V. et al. "Structure-based protein function prediction using graph convolutional networks." *Nature Communications* 12:3168, 2021. DOI: https://doi.org/10.1038/s41467-021-23303-9
- [Metagenomic-DeepFRI GitHub](https://github.com/bioinf-mcb/Metagenomic-DeepFRI)
- [mDeepFRI documentation](https://metagenomic-deepfri.readthedocs.io/)
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
