---
doc_id: workflow-012
domain: structure-prediction
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  Predicts 3D biomolecular structures using Boltz-2 diffusion models,
  with confidence metrics and an interactive HTML dashboard.
tags: [structure-prediction, boltz-2, protein-folding, diffusion-model, plddt]
---

# Workflow 012: Boltz-2 Structure Prediction

AI-based 3D biomolecular structure prediction using Boltz-2, a diffusion generative model that predicts structures for proteins, nucleic acids, small molecule ligands, and their complexes. The workflow validates input sequences, runs GPU-accelerated structure prediction, and generates an interactive HTML dashboard with confidence metrics and 3D visualization.

## Overview

Boltz-2 uses a diffusion-based architecture to generate 3D structures by iteratively denoising atomic coordinates from random noise. It produces multiple structural models ranked by confidence, along with per-residue quality metrics (pLDDT, pTM, iPTM) and error estimates (PAE, PDE). The model supports proteins, RNA, DNA, small molecule ligands (via SMILES or CCD codes), and multi-chain complexes in a single unified framework (Wohlwend et al., 2024; Passaro et al., 2025).

## When to use this workflow

Use this workflow when you have one or more biomolecular sequences (protein, RNA, DNA, or ligand) and want to predict their 3D structure or complex geometry. Input is provided as a Boltz-2 YAML file specifying sequences and entity types. The workflow is suitable for single-chain fold prediction, multi-chain complex assembly, and protein-ligand binding pose prediction.

Do not use this workflow if you want to compare predictions from multiple tools — use workflow-018 (Boltz-2 vs Chai-1 comparison) instead. Do not use this workflow for molecular docking against a known protein structure — use workflow-004 (AutoDock Vina) or workflow-016 (DiffDock-PP). For ADMET property prediction from SMILES, use workflow-014.

## Architecture and data flow

```text
input.yaml ──> [01: Validate] ──> [02: Predict] ──> [03: Report]
                                       |                  |
                                  *.pdb, *.json       dashboard.html
                                  pae/pde/plddt.npz
```

Nodes run sequentially: 01 → 02 → 03.

## Input requirements

- **Format:** Boltz-2 YAML file specifying sequences with entity types and chain IDs.
- **Supported entity types:** `protein`, `rna`, `dna`, `ligand` (SMILES), `ccd` (Chemical Component Dictionary).
- **Placement:** Place the YAML file in `input_files/`.
- **Sample data:** `input_files/prot.yaml` — human hemoglobin alpha subunit (141 residues, single chain).

Example input format:
```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: MVLSPADKTNVKAAWGKVGA...
```

## Workflow nodes

### Node 01: Sequence Upload

**Goal:** Validate the input YAML file and extract sequence metadata.

**Process:** Parses the YAML file, checks for required fields (`version`, `sequences`), validates that each entity has a non-empty sequence and recognized type (protein, rna, dna, ligand, ccd). Outputs the validated YAML and a JSON summary with entity metadata (type, ID, sequence length).

**Scientific notes:** Boltz-2 requires structured input specifying entity types explicitly, unlike FASTA-only tools, because the model uses type-specific featurization for proteins vs. nucleic acids vs. small molecules.

**Outputs:**
- `*.yaml` — validated input file
- `input_summary.json` — entity metadata (types, chain IDs, lengths)

### Node 02: Structure Prediction

**Goal:** Run Boltz-2 structure prediction with GPU acceleration.

**Process:** Invokes the `boltz predict` CLI with configurable parameters: `diffusion_samples` (number of independent models to generate), `recycling_steps` (refinement iterations through the model trunk), and optional MSA server for evolutionary information. Each diffusion sample produces an independent structure through iterative denoising. Outputs PDB structure files, confidence metrics (JSON), and quality matrices (PAE, PDE, pLDDT as NPZ files).

**Scientific notes:** The diffusion process starts from random atomic coordinates and iteratively refines them through learned denoising steps. Multiple diffusion samples explore different regions of the structural landscape — the model may converge to different conformations, especially for flexible regions. Recycling steps pass the output back through the model trunk for further refinement, similar to AlphaFold2's recycling mechanism. The MSA server provides multiple sequence alignment information from evolutionary databases, which significantly improves prediction accuracy for proteins with many homologs.

**Outputs:**
- `*_model_N.pdb` — predicted 3D structures (one per diffusion sample)
- `confidence_*.json` — confidence score, pLDDT, pTM, iPTM, PDE per model
- `pae_*.npz`, `pde_*.npz`, `plddt_*.npz` — per-residue quality matrices

### Node 03: Report Generation

**Goal:** Generate an interactive HTML dashboard with quality metrics and 3D structure visualization.

**Process:** Aggregates confidence metrics from all models, ranks them by confidence score, generates Plotly charts (confidence per model, pLDDT distribution, quality breakdown), embeds a 3Dmol.js interactive 3D viewer with multiple display styles (cartoon, stick, sphere, surface), and produces evidence-based recommendations based on quality thresholds.

**Scientific notes:** The dashboard enables rapid assessment of prediction reliability. Model consistency (coefficient of variation across samples) indicates whether the prediction is well-determined or if the model samples diverse conformations, which may suggest intrinsic flexibility or prediction uncertainty.

**Outputs:**
- `boltz_dashboard_*.html` — self-contained interactive HTML dashboard

## Parameters

### diffusion_samples

| Value | Description |
|-------|-------------|
| `10` (default) | Generates 10 independent structural models. Suitable for most predictions. |
| `2`–`5` | Faster testing with fewer models. |
| `20`+ | More thorough sampling for difficult targets or flexible complexes. |

**Trade-off:** More samples increase the chance of finding the best conformation but increase runtime linearly. Each sample is an independent diffusion trajectory.

**Test vs production:** Use 2–5 for testing, 10+ for production.

### recycling_steps

| Value | Description |
|-------|-------------|
| `5` (default) | Number of refinement iterations through the model trunk. |
| `3` | Faster, suitable for well-folded single-domain proteins. |
| `7`–`10` | More refinement for large or multi-domain complexes. |

**Trade-off:** More recycling steps improve structure quality but increase compute time per sample.

### use_msa_server

| Value | Description |
|-------|-------------|
| `true` (default) | Uses ColabFold MSA server for evolutionary information. Recommended for best accuracy. |
| `false` | Runs without MSA. Faster but may produce lower-quality predictions, especially for proteins with few homologs. |

**Trade-off:** MSA provides evolutionary covariance signals that improve contact prediction. Disabling it is faster but sacrifices accuracy, particularly for novel folds.

## Outputs and interpretation

### Confidence score

Aggregated ranking metric: 0.8 × complex_pLDDT + 0.2 × iPTM. Range 0–1, higher is better. Used to rank models from the same prediction run. Scores > 0.7 indicate a reliable prediction.

### pLDDT (predicted Local Distance Difference Test)

Per-residue confidence in local structure accuracy. Scale 0–1 (some tools report 0–100).

| Range | Interpretation |
|-------|---------------|
| > 0.9 | Very high confidence — backbone and side chains accurately predicted |
| 0.7–0.9 | Confident — correct backbone, some side chain uncertainty |
| 0.5–0.7 | Low confidence — uncertain local structure |
| < 0.5 | Very low confidence — likely intrinsically disordered region |

### pTM (predicted TM-score)

Confidence in the overall fold topology. Range 0–1. Values > 0.5 indicate the predicted fold is roughly correct; > 0.8 indicates high-confidence fold prediction. Relevant for single-chain predictions.

### iPTM (interface predicted TM-score)

Confidence in the predicted interface geometry for multi-chain complexes. Range 0–1. Higher values indicate reliable inter-chain contact prediction. Contributes to the ranking score alongside pLDDT (see confidence score formula above) for assessing complex modeling quality.

### PAE (Predicted Aligned Error)

Matrix of predicted positional errors (in Å) between all residue pairs. Low PAE between two residues indicates high confidence in their relative positioning. Useful for assessing domain-domain relationships — high pLDDT within each domain but high PAE between domains suggests uncertain relative orientation.

### PDE (Predicted Distance Error)

Predicted error in pairwise distances (in Å). Lower values indicate more accurate distance predictions. Summarized as a per-complex mean.

## Quick start

### Running with Docker

| Node | Image |
|------|-------|
| 01, 02 | `ghcr.io/chiral-data/boltz:2025_09_05` |
| 03 | `ghcr.io/chiral-data/boltz_report:2026_02_13` |

GPU is required for Node 02 (structure prediction).

### Running on Silva

1. Select workflow-012 from the workflow list
2. Upload your Boltz-2 YAML file to `input_files/`
3. Adjust parameters if needed (see Parameters section)
4. Click Run

### Test vs production settings

| Setting | Test | Production |
|---------|------|------------|
| `diffusion_samples` | `2` | `10`+ |
| `recycling_steps` | `3` | `5`+ |
| `use_msa_server` | `true` | `true` |

A successful test run with the hemoglobin alpha subunit (141 residues) generates 10 PDB models, confidence metrics, and an interactive dashboard in approximately 4 minutes on GPU.

## References

- Wohlwend J, Corso G, Passaro S, Getz M, Reveiz M, Leidal K, Swiderski W, Atkinson J, Portnoi T, Chinn I, Silterra J, Jaakkola T, Barzilay R. "Boltz-1: Democratizing Biomolecular Interaction Modeling." *bioRxiv*, 2024. DOI: https://doi.org/10.1101/2024.11.19.624167
- Passaro S, Corso G, Wohlwend J, Reveiz M, Thaler S, Somnath VR, Getz M, Portnoi T, Roy H, Stark H, Kwabi-Addo B, Beaini D, Jaakkola T, Barzilay R. "Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction." *bioRxiv*, 2025. DOI: https://doi.org/10.1101/2025.06.14.659707
- [Boltz source code](https://github.com/jwohlwend/boltz)
