---
doc_id: workflow-010
domain: structure-prediction
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  Predicts antibody Fv structure from paired heavy/light chain FASTA
  sequences using ABodyBuilder3 (ABB3).
tags: [antibody, structure-prediction, abodybuilder3, protein-folding]
---

# Workflow 010: Antibody Structure Prediction (ABodyBuilder3)

Predicts the three-dimensional Fv structure of an antibody from paired
heavy and light chain amino-acid sequences using ABodyBuilder3 (ABB3).
Supports both the standard ABB3 model and the language-model enhanced
ABB3-LM variant with ProtT5 embeddings.

## Overview

ABodyBuilder3 is a deep-learning method for antibody structure prediction
built on the ImmuneBuilder framework. Given a pair of heavy and light chain
sequences, it predicts full-atom coordinates for the antibody variable
fragment (Fv), including the complementarity-determining regions (CDRs)
that determine antigen binding specificity.

The workflow accepts standard FASTA files as input and produces PDB
structure files together with an interactive HTML report for visual
inspection.

**Citation:** Kenlay, H. et al. "ABodyBuilder3: improved and scalable
antibody structure predictions." *Bioinformatics*, 40(10):btae576, 2024.

## When to use this workflow

Use this workflow when you have paired heavy and light chain sequences for
one or more antibodies and need predicted Fv structures. Input sequences
must use standard 20 amino-acid single-letter codes. The workflow handles
both single pairs and batches (one entry per FASTA record, matched by
position).

This workflow predicts the **Fv region only** (variable heavy + variable
light). It does not model full-length antibodies, constant domains, or
Fc regions. For general protein structure prediction (non-antibody), use
workflow-026 (ColabFold) or workflow-018 (Boltz) instead.

## Architecture and data flow

```text
heavy.fasta  --+
               +--> [01: Input Preparation] --> *.pt
light.fasta  --+                                 |
                                      (optional) |
                                                 v
                                   [02: PLM Embedding] --> *.pt (+ plm_embedding)
                                                              |
                        (checkpoint baked in image) ----------+
                                                              v
                                          [03: Structure Prediction] --> *.pdb
                                                                          |
                                                                          v
                                                  [04: Visualization Report] --> report.html
```

Nodes run sequentially: 01 -> 03 -> 04 (plain ABB3) or 01 -> 02 -> 03 -> 04 (ABB3-LM).

## Input requirements

- **heavy.fasta** -- one or more heavy chain sequences in FASTA format
- **light.fasta** -- matching light chain sequences in FASTA format
- Sequences must pair by position (1st heavy with 1st light, etc.)
- Only standard 20 amino-acid single-letter codes (ACDEFGHIKLMNPQRSTVWY)
- Place files in `input_files/`

Test data is included: `input_files/heavy.fasta` and `input_files/light.fasta`
contain the 6yio antibody Fv sequences from the ABB3 example notebook.

## Workflow nodes

### Node 01: Input Preparation

**Goal:** Validate and pair heavy/light chain sequences for downstream prediction.

**Process:** Reads both FASTA files, validates that all residues are standard
amino acids, pairs chains by position, and converts each pair to an ABB3
`ab_input` tensor dictionary using `string_to_input()`. Saves one `.pt` file
per pair.

**Scientific notes:** Strict validation rejects ambiguous residues (B, Z, X, U)
because ABB3 was trained on canonical amino acids only. Order-based pairing
assumes the user has pre-aligned their FASTA records.

**Outputs:**
- `<heavy_id>-<light_id>.pt` -- serialized dict with `id`, `heavy`, `light`, and `ab_input` tensors

### Node 02: PLM Embedding (Optional)

**Goal:** Generate ProtT5 protein language model embeddings to enable ABB3-LM mode.

**Process:** Loads each `.pt` file from Node 01, runs the concatenated
heavy+light sequence through ProtT5-XL-U50 to produce per-residue embeddings
of dimension 1024, and saves an enriched `.pt` file with the `plm_embedding`
key added. Supports pre-computed embedding caching.

**Scientific notes:** ProtT5 embeddings capture evolutionary and structural
information learned from billions of protein sequences. ABB3-LM uses these
as additional input features, which can improve prediction accuracy for
antibodies with unusual CDR sequences. The model is ~3 GB and requires a GPU
with sufficient VRAM.

**Outputs:**
- `<pair_id>.pt` -- same as Node 01 output, extended with `plm_embedding` tensor (L x 1024)

### Node 03: Structure Prediction

**Goal:** Run ABB3 inference to predict the 3D Fv structure.

**Process:** Loads the ABB3 checkpoint, batches each input, runs the forward
pass, reconstructs atom37 coordinates via `add_atom37_to_output()`, and
writes one PDB file per antibody pair. Per-residue pLDDT confidence scores
are extracted from the model output and written into the B-factor column.
Automatically detects whether PLM embeddings are present (ABB3-LM) or absent
(plain ABB3).

**Scientific notes:** ABB3 uses a structure module similar to AlphaFold2's
IPA (Invariant Point Attention) architecture, specialized for antibody Fv
regions. The model predicts backbone and side-chain atom positions. The
`plddt-loss` checkpoint is used for standard mode; the `language-loss`
checkpoint for ABB3-LM mode. pLDDT (predicted local-distance difference test)
is the model's per-residue confidence score (0–100).

**Outputs:**
- `<pair_id>.pdb` -- predicted Fv structure in PDB format (B-factor = pLDDT)

### Node 04: Visualization Report

**Goal:** Generate an interactive HTML report with 3D viewer and pLDDT analysis.

**Process:** Copies all PDB files to the output directory and generates a
single HTML file with a Mol\* (PDBe) 3D viewer, a per-chain summary table
(residue count, average/minimum pLDDT, confidence badge), a per-residue pLDDT
line chart with confidence bands, and a pLDDT color legend. The report header
indicates whether ABB3-LM or plain ABB3 mode was used. An internet connection
is required to load the Mol\* viewer from CDN at view time.

**Scientific notes:** pLDDT ≥ 90 (dark blue) indicates very high confidence;
70–90 (cyan) is generally accurate; 50–70 (yellow) should be treated with
caution; < 50 (orange) indicates low confidence, often disordered loops. CDR
loops, especially CDR H3, tend to have lower pLDDT than framework regions.

**Outputs:**
- `report.html` -- interactive HTML report (requires internet for Mol\* viewer)
- `<pair_id>.pdb` -- copied PDB files

## Parameters

### DEVICE

| Value | Description |
|-------|-------------|
| `cuda` (default) | Use GPU for inference. Required for reasonable speed. |
| `cpu` | CPU-only mode. Much slower but works without a GPU. |

**Trade-off:** GPU is 10-50x faster. Use `cpu` only for debugging or when no GPU is available.

### USE_PLM

| Value | Description |
|-------|-------------|
| `0` | Plain ABB3 mode. Skip Node 02. |
| `1` (default) | ABB3-LM mode. Run Node 02 to generate ProtT5 embeddings. |

**Trade-off:** ABB3-LM may produce better predictions for unusual CDR sequences
but requires additional GPU memory (~3 GB for ProtT5) and computation time.

**Test vs production:** Default (`1`) enables ABB3-LM for best accuracy.
Set to `0` to skip ProtT5 embeddings and run faster with plain ABB3.

### REPORT_TITLE

- **Type:** string
- **Default:** `ABB3 Structure Predictions`
- **Description:** Title displayed in the HTML report header.
- **Guidance:** Change to describe your experiment, e.g., "Anti-VEGF Fv Predictions".

## Outputs and interpretation

### PDB structure files

Each `<pair_id>.pdb` contains the predicted full-atom Fv structure. The
heavy chain is listed first, followed by the light chain. The B-factor
column contains per-residue pLDDT scores (0–100). These files can be
opened in PyMOL, ChimeraX, or any PDB viewer for detailed analysis.
To color by pLDDT in PyMOL: `spectrum b, blue_white_red, minimum=50, maximum=100`.

**Caveats:** ABB3 predicts the Fv region only. Coordinates are not
experimentally determined -- use predicted structures for hypothesis
generation, not as ground truth. For high-confidence applications,
validate predictions against experimental data when available.

### HTML report

The `report.html` file provides a quick visual overview with:
- **Mol\* 3D viewer** -- interactive structure visualization
- **Summary table** -- per-chain residue count, avg/min pLDDT, confidence badge
- **Per-residue pLDDT chart** -- line plot with confidence bands per chain
- **Mode indicator** -- header badge showing ABB3-LM or plain ABB3

An internet connection is required to load the Mol\* viewer library from CDN.

## Quick start

### Running with Docker

Build the image from the workflow directory:

```bash
docker build -t abodybuilder3:2026_06_15 .
```

Run each node sequentially:

```bash
# Node 01 -- Input Preparation
docker run --rm --gpus all \
  -v $(pwd)/input_files:/workflow/01-input-preparation/inputs \
  -v $(pwd)/results/01:/workflow/01-input-preparation/outputs \
  -w /workflow/01-input-preparation \
  abodybuilder3:2026_06_15 bash run.sh

# Node 03 -- Structure Prediction (plain ABB3, skip Node 02)
docker run --rm --gpus all \
  -v $(pwd)/results/01:/workflow/03-structure-prediction/inputs \
  -v $(pwd)/results/03:/workflow/03-structure-prediction/outputs \
  -w /workflow/03-structure-prediction \
  abodybuilder3:2026_06_15 bash run.sh

# Node 04 -- Visualization Report
docker run --rm --gpus all \
  -v $(pwd)/results/03:/workflow/04-visualization-report/inputs \
  -v $(pwd)/results/04:/workflow/04-visualization-report/outputs \
  -w /workflow/04-visualization-report \
  abodybuilder3:2026_06_15 bash run.sh
```

### Running on Silva

1. Select "Antibody Structure Prediction (ABB3)" from the workflow list
2. Upload your `heavy.fasta` and `light.fasta` files
3. Adjust parameters if needed (see Parameters section)
4. Click Run

### Test vs production settings

| Setting | Test (default) | Production |
|---------|---------------|------------|
| DEVICE | `cuda` | `cuda` |
| USE_PLM | `1` | `1` (or `0` for faster plain ABB3) |
| REPORT_TITLE | `ABB3 Structure Predictions` | Descriptive experiment name |

Test data (6yio Fv) produces a single PDB in under 2 minutes on GPU.
A successful run produces `results/04/report.html` with one structure card.

## Troubleshooting

**CUDA out of memory (Node 02)**
ProtT5 requires ~3 GB VRAM. If running ABB3-LM mode on a GPU with limited
memory, try setting `DEVICE=cpu` for Node 02 only (slower but avoids OOM).

**Empty outputs from Node 01**
Check that both FASTA files have the same number of sequences.
Sequences must be paired by position.

**"invalid characters" error**
Input sequences contain non-standard residues. Remove or replace ambiguous
codes (B, Z, X, U) with standard amino acids before running.

## References

- Kenlay, H. et al. "ABodyBuilder3: improved and scalable antibody structure predictions." *Bioinformatics*, 40(10):btae576, 2024. DOI: https://doi.org/10.1093/bioinformatics/btae576
- [ABodyBuilder3 GitHub repository](https://github.com/Exscientia/abodybuilder3)
- Elnaggar, A. et al. "ProtTrans: Toward Understanding the Language of Life Through Self-Supervised Learning." *IEEE TPAMI*, 44(10):7112-7127, 2022.
