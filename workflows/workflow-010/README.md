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
writes one PDB file per antibody pair using `output_to_pdb()`. Automatically
detects whether PLM embeddings are present (ABB3-LM) or absent (plain ABB3).

**Scientific notes:** ABB3 uses a structure module similar to AlphaFold2's
IPA (Invariant Point Attention) architecture, specialized for antibody Fv
regions. The model predicts backbone and side-chain atom positions. The
`plddt-loss` checkpoint is used for standard mode; the `language-loss`
checkpoint for ABB3-LM mode.

**Outputs:**
- `<pair_id>.pdb` -- predicted Fv structure in PDB format

### Node 04: Visualization Report

**Goal:** Generate a self-contained interactive HTML report with 3D viewers.

**Process:** Copies all PDB files to the output directory and generates a
single HTML file embedding each structure as a py3Dmol viewer. PDB data is
inlined as JavaScript template literals so the report works offline in any
modern browser. Each structure card includes a 3D viewer (cartoon coloring
by spectrum) and a download link.

**Scientific notes:** Spectrum coloring maps the rainbow from N-terminus (blue)
to C-terminus (red), making it easy to identify heavy vs. light chain regions
and locate CDR loops.

**Outputs:**
- `report.html` -- interactive HTML report
- `<pair_id>.pdb` -- copied PDB files (for download links)

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
| `0` (default) | Plain ABB3 mode. Skip Node 02. |
| `1` | ABB3-LM mode. Run Node 02 to generate ProtT5 embeddings. |

**Trade-off:** ABB3-LM may produce better predictions for unusual CDR sequences
but requires additional GPU memory (~3 GB for ProtT5) and computation time.

**Test vs production:** Default (`0`) is fine for both testing and production.
Set to `1` when working with non-standard or engineered antibodies where
extra accuracy matters.

### REPORT_TITLE

- **Type:** string
- **Default:** `ABB3 Structure Predictions`
- **Description:** Title displayed in the HTML report header.
- **Guidance:** Change to describe your experiment, e.g., "Anti-VEGF Fv Predictions".

## Outputs and interpretation

### PDB structure files

Each `<pair_id>.pdb` contains the predicted full-atom Fv structure. The
heavy chain is listed first, followed by the light chain. These files can
be opened in PyMOL, ChimeraX, or any PDB viewer for detailed analysis.

**Caveats:** ABB3 predicts the Fv region only. Coordinates are not
experimentally determined -- use predicted structures for hypothesis
generation, not as ground truth. For high-confidence applications,
validate predictions against experimental data when available.

### HTML report

The `report.html` file provides a quick visual overview of all predicted
structures in a browser. Each card shows a 3D viewer and a PDB preview.
No internet connection is required to view the report after generation
(py3Dmol is loaded from CDN at generation time and embedded).

## Quick start

### Running with Docker

Build the image from the workflow directory:

```bash
docker build -t abodybuilder3:latest .
```

Run each node sequentially:

```bash
# Node 01 -- Input Preparation
docker run --rm --gpus all \
  -v $(pwd)/input_files:/workflow/01-input-preparation/inputs \
  -v $(pwd)/results/01:/workflow/01-input-preparation/outputs \
  -w /workflow/01-input-preparation \
  abodybuilder3:latest bash run.sh

# Node 03 -- Structure Prediction (plain ABB3, skip Node 02)
docker run --rm --gpus all \
  -v $(pwd)/results/01:/workflow/03-structure-prediction/inputs \
  -v $(pwd)/results/03:/workflow/03-structure-prediction/outputs \
  -w /workflow/03-structure-prediction \
  abodybuilder3:latest bash run.sh

# Node 04 -- Visualization Report
docker run --rm --gpus all \
  -v $(pwd)/results/03:/workflow/04-visualization-report/inputs \
  -v $(pwd)/results/04:/workflow/04-visualization-report/outputs \
  -w /workflow/04-visualization-report \
  abodybuilder3:latest bash run.sh
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
| USE_PLM | `0` | `0` or `1` depending on need |
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
