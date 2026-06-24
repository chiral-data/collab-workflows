---
doc_id: workflow-028
domain: molecular-docking
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  End-to-end AI-driven pipeline that takes a protein sequence and ligand SMILES
  and produces docked 3D poses — without requiring a crystal structure.
  Chains Boltz-2 structure prediction, P2Rank pocket detection, and
  Uni-Mol Docking V2.
tags: [molecular-docking, structure-prediction, boltz2, p2rank, unimol, virtual-screening, pocket-detection]
---

# Workflow 028: AI-Driven Structure Prediction → Pocket Detection → Docking

An end-to-end pipeline that produces docked small-molecule poses directly from a protein sequence, with no experimental structure required. **Boltz-2** predicts the 3D protein fold in holo mode (with the ligand SMILES included), **P2Rank** identifies druggable binding pockets using an AlphaFold-aware pLDDT profile, and **Uni-Mol Docking V2** docks the ligand into the highest-confidence pocket. A final HTML report integrates structure confidence metrics, pocket rankings, and interactive 3D visualization.

## Overview

Boltz-2 (Wohlwend et al. 2024; Passaro et al. 2025) is a diffusion-based co-folding model that predicts protein structure, ligand pose, binding affinity, and confidence metrics (pLDDT, PAE) in a single forward pass. Including the ligand SMILES in the input YAML enables *holo-mode* prediction, which constrains the pocket geometry toward a ligand-bound conformation and substantially improves downstream docking accuracy over an apo fold.

P2Rank (Jakubec et al. 2025) predicts pocket centers from 3D structure using a random-forest model. The `-c alphafold` profile is required when the input comes from a predicted structure: it interprets the B-factor column as pLDDT rather than thermal displacement.

Uni-Mol Docking V2 (Zhou et al. 2022, ICLR 2023) achieves 77.6% RMSD < 2 Å on PoseBusters — significantly outperforming traditional physics-based methods on predicted structures. It outputs 3D poses only; no binding affinity score is written to the SDF file.

## When to use this workflow

Use this workflow when you have a protein sequence and a ligand SMILES and want to generate docked poses without an experimental structure. It is particularly well suited for novel targets where no PDB entry exists, or for comparing predicted-structure docking against a known crystal structure as a validation exercise.

**When NOT to use:** If you already have a high-quality crystal structure (≤ 2.5 Å, co-crystallized ligand), use **workflow-025** (Vina vs GNINA comparison) instead — it skips the expensive structure-prediction step and uses a more mature docking pipeline. This workflow is also not appropriate for large-scale virtual screening (hundreds of ligands): Uni-Mol Docking V2 is slower than AutoDock Vina/GNINA and requires GPU. For high-throughput SMILES-based property prediction without docking, see **workflow-014** (ADMET-AI).

## Architecture and data flow

```
input_files/protein.yaml  (sequence + ligand SMILES + affinity block)
         │
         ▼
01_boltz2_predict ──► *_model_*.cif + plddt_*.npz + pae_*.npz + affinity_*.json
         │
         ▼
02_p2rank_pocket_find  (-c alphafold) ──► predictions.csv + residues.csv
         │
         ▼
03_pocket_qc_grid ──► receptor.pdb + ligand.sdf + grid.json + pocket_qc.json
         │
         ▼
04_unimol_docking ──► docked_poses.sdf + docking_summary.json
         │
         ▼
05_generate_report ──► report.html
```

Nodes run sequentially: 01 → 02 → 03 → 04 → 05. Node 05 reads outputs from all prior nodes.

## Input requirements

A single **YAML file** at `input_files/protein.yaml` with three sections:

- **`sequences`** — one protein entry with a `protein` key containing the amino acid sequence.
- **`ligand`** — SMILES string for the molecule to dock.
- **`properties`** — `affinity: true` to enable Boltz-2 affinity prediction.

Input must be YAML — FASTA mode does not support ligand specification or affinity prediction. A test input for EGFR kinase domain (PDB: 1XKK) with Lapatinib is included at `input_files/protein.yaml`.

## Workflow nodes

### Node 01: Boltz-2 Structure Prediction

**Goal:** Predict the 3D protein–ligand co-structure and confidence metrics from the input YAML.

**Process:** Validates the YAML (checks sequence length, ligand SMILES parsability, required keys) and runs Boltz-2 with the specified number of diffusion samples and recycling steps. Output is mmCIF format — the `--output_format pdb` flag is broken for protein-ligand complexes (boltz #298) and must not be used. Node selects the model with the highest mean pLDDT and writes `selected_model_id.txt`.

**Scientific notes:** Including the ligand SMILES in the YAML places Boltz-2 in holo mode, which conditions the predicted fold on the ligand geometry and produces a more accurate binding-site conformation than apo prediction alone. P2Rank can accept mmCIF directly, so no format conversion is needed at this stage.

**Outputs:**
- `*_model_*.cif` — predicted structure(s) in mmCIF format
- `plddt_*.npz` — per-residue pLDDT confidence (0–1 scale)
- `pae_*.npz` — predicted aligned error matrix
- `confidence_*.json` — overall confidence metrics
- `affinity_*.json` — predicted binding affinity and confidence
- `input_summary.json` — parsed input metadata (sequence length, SMILES, etc.)

### Node 02: P2Rank Pocket Detection

**Goal:** Identify ligand-binding pocket candidates on the best predicted structure.

**Process:** Runs P2Rank with the `-c alphafold` configuration profile, which treats the B-factor column as pLDDT rather than thermal displacement. Outputs a ranked list of pocket centers with scores and the residue-level pocket membership used for pLDDT QC in node 03.

**Scientific notes:** The `-c alphafold` profile is required for any structure from Boltz-2 or AlphaFold. Without it, P2Rank misinterprets pLDDT values as thermal B-factors, which degrades pocket scoring.

**Outputs:**
- `predictions.csv` — ranked pocket centers (`center_x/y/z`), `score`, `probability`
- `residues.csv` — per-residue pocket membership assignments
- `selected_structure.cif` — the mmCIF file passed to P2Rank
- `selected_model_id.txt` — identifier of the selected Boltz-2 model

### Node 03: Pocket QC and Grid Preparation

**Goal:** Filter pockets by structural confidence, convert file formats, and write the docking grid.

**Process:** Reads per-residue pLDDT from `plddt_*.npz` and computes the mean, min, and std pLDDT across the residues of each P2Rank pocket. The selected pocket (by `pocket_rank`) is flagged if its mean pLDDT falls below `plddt_threshold`, but docking proceeds regardless so the report can show the confidence context. Converts the Boltz-2 mmCIF to PDB using `gemmi` (BioPython's `MMCIFIO` is incompatible with Boltz-2's CIF format). Extracts the ligand SMILES from the input YAML, validates it with RDKit, and writes a 3D-embedded SDF. Writes `grid.json` with the pocket center and the fixed `box_size` applied uniformly to all three dimensions.

**Scientific notes:** The ≥ 70 pLDDT threshold is well-supported: AlphaFold's official scale defines ≥ 70 as "confident, backbone correct," and PrankWeb 4 uses it as the default filtering cutoff. However, Eguida & Rognan (2023, *JCIM*) found that high pLDDT is necessary but not sufficient — four of five worst-performing docking targets still had binding-site pLDDT ≥ 70. The per-pocket distribution (mean, min, std) is reported rather than a binary pass/fail.

**Outputs:**
- `receptor.pdb` — protein structure in PDB format (required by Uni-Mol)
- `ligand.sdf` — 3D-embedded ligand conformer
- `grid.json` — docking grid: `{center_x, center_y, center_z, size_x, size_y, size_z}`
- `pocket_qc.json` — pLDDT statistics per pocket and QC pass/fail flag

### Node 04: Uni-Mol Docking V2

**Goal:** Generate docked 3D poses of the ligand in the selected pocket.

**Process:** Validates the receptor PDB, ligand SDF (pre-validating with RDKit to catch silently-dropped molecules; Uni-Mol #281), and model weight files. Invokes `interface/demo.py` from the cloned Uni-Mol repository with `--mode single`, passing `grid.json` directly via `--input-docking-grid`. The `--steric-clash-fix` and `--cluster` flags are enabled for improved pose quality. Collected SDF poses are concatenated into `docked_poses.sdf`.

**Scientific notes:** Uni-Mol Docking V2 outputs 3D poses only — the internal `prmsd_score` used for pose ranking is not written to the output file. Binding affinity estimation should use Boltz-2's `affinity_*.json` as a complementary signal.

**Outputs:**
- `docked_poses.sdf` — concatenated docked conformations
- `docking_summary.json` — receptor path, ligand SMILES, grid, number of poses requested and generated

### Node 05: Pipeline Report

**Goal:** Produce a self-contained HTML dashboard summarizing the full pipeline run.

**Process:** Reads outputs from all prior nodes and renders a Bootstrap 5 + Plotly dashboard with an embedded Mol* 3D viewer. Includes a pipeline summary, Boltz-2 confidence panels, P2Rank pocket table, and docking results section with a caveats note on the absence of Uni-Mol affinity scores.

**Outputs:**
- `report.html` — self-contained HTML; can be opened in any browser without re-running code

## Parameters

### `diffusion_samples`

- **Type:** integer
- **Default:** `2`
- **Description:** Number of Boltz-2 structural models to generate.
- **Guidance:** Use 2–3 for testing. For publication-quality results or when exploring conformational diversity, use 10 or more. Node 02 selects the model with the highest mean pLDDT.

### `recycling_steps`

- **Type:** integer
- **Default:** `3`
- **Description:** Number of Boltz-2 recycling iterations through the model trunk.
- **Guidance:** Use 3 for testing. Increase to 5–10 for production runs; more recycling steps generally improve structure quality at the cost of runtime.

### `use_msa_server`

- **Type:** boolean
- **Default:** `true`
- **Description:** Whether to query the ColabFold MSA server for evolutionary co-variation information.
- **Guidance:** Keep `true` for production. Set `false` only for rapid local testing on sequences with known structure — MSA significantly improves prediction accuracy.

### `plddt_threshold`

- **Type:** float
- **Default:** `70.0`
- **Description:** Minimum mean pLDDT (0–100 scale) for a pocket to pass QC. Pockets below this threshold are flagged in the report but docking still proceeds.
- **Guidance:** 70 is the standard AlphaFold "confident backbone" cutoff. Lower to 60 for exploratory runs on difficult targets; raise to 80 when structural confidence must be high before acting on docking results.

### `box_size`

- **Type:** float
- **Default:** `22.5`
- **Description:** Docking grid box size in Å, applied uniformly to all three dimensions (`size_x = size_y = size_z = box_size`). This is a fixed user-configurable parameter — it is not computed from the pocket spatial extent.
- **Guidance:** 20–25 Å covers most drug-like molecules in a typical pocket. Increase to 28–30 Å for larger ligands or if the pocket is shallow and wide. Reducing below 18 Å risks clipping part of the binding site.

### `pocket_rank`

- **Type:** integer
- **Default:** `1`
- **Description:** Which P2Rank pocket to dock into (1 = highest P2Rank score).
- **Guidance:** Increase if the top-ranked pocket fails pLDDT QC or is in a structurally uncertain region. Check `predictions.csv` and `pocket_qc.json` to compare pocket scores and confidence.

### `num_poses`

- **Type:** integer
- **Default:** `10`
- **Description:** Number of docked conformations to generate per ligand (maps to `--conf-size` in `interface/demo.py`).
- **Guidance:** 10 is sufficient for inspecting the top poses. Increase to 20–50 when clustering or ensemble docking is needed.

## Outputs and interpretation

### `report.html`

Self-contained HTML dashboard with: (1) pipeline summary table, (2) Boltz-2 pLDDT histogram and PAE heatmap, (3) P2Rank pocket table with pLDDT statistics per pocket, (4) Mol* 3D viewer with protein colored by pLDDT and the top docking pose overlaid, (5) methods and caveats section. Open in any browser.

### `docked_poses.sdf`

Concatenated SDF of all generated docking poses. Poses are ordered by Uni-Mol's internal `prmsd_score` (lower = better predicted RMSD from true binding mode), but this score is **not written to the file**. Inspect pose geometry in the report's 3D viewer or with a molecular visualization tool; do not rank poses by SDF record order alone.

### `docking_summary.json`

Metadata record: receptor path, ligand SMILES, grid center and size, number of poses requested and generated. Useful for provenance tracking and debugging.

### `pocket_qc.json`

Per-pocket pLDDT statistics (mean, min, std across pocket residues) and a `selected_pocket_passes_qc` boolean. A `false` value means the selected pocket has lower structural confidence; interpret docking results cautiously and consider increasing `pocket_rank` to try a higher-confidence pocket.

### `affinity_*.json` (from node 01)

Boltz-2's predicted binding affinity and confidence score. This is the only affinity signal in the pipeline — Uni-Mol Docking V2 does not output binding energy. Use in combination with P2Rank `score` and `probability` to prioritize results.

## Quick start

### Running on Silva

1. Select workflow **028** from the workflow list.
2. Upload `input_files/protein.yaml` with your protein sequence and ligand SMILES.
3. Adjust parameters (see Parameters section) — defaults are configured for a fast test run.
4. Click **Run**.

### Running with Docker

Each node has its own Dockerfile. Node 01 reuses the pre-built image `ghcr.io/chiral-data/boltz:2025_09_05`. Node 04 builds from `dptechnology/unicore:latest-pytorch1.12.1-cuda11.6-rdma`.

**Note:** Uni-Mol Docking V2 model weights (464 MB) are not distributed via Docker, pip, or HuggingFace. Download them from the Dropbox link in the [Uni-Mol Docking V2 README](https://github.com/deepmodeling/Uni-Mol/tree/main/unimol_docking_v2) and place them at the path set by `PARAM_WEIGHTS_PATH` (default: `/opt/unimol_weights`).

### Test vs production settings

| Setting | Test (default) | Production |
|---------|---------------|------------|
| `diffusion_samples` | `2` | `10`+ |
| `recycling_steps` | `3` | `5`–`10` |
| `use_msa_server` | `true` | `true` |
| `num_poses` | `10` | `20`–`50` |

The included test input (EGFR kinase domain + Lapatinib, PDB: 1XKK) has a known crystal structure at 2.40 Å, making it suitable for validating that Boltz-2 recovers the correct fold and Uni-Mol places poses in the correct pocket.

## Troubleshooting

**Uni-Mol model weights not found**
Node 04 will exit with a clear error if the weights directory is missing or empty. Download the weights from the Dropbox link in the Uni-Mol Docking V2 README and set `PARAM_WEIGHTS_PATH` to the directory containing the `.pt`/`.pkl` files.

**Node 04 produces an empty `docked_poses.sdf`**
Usually caused by an invalid ligand SDF. RDKit sanitization inside Uni-Mol can silently drop molecules that fail valence checks (Uni-Mol #281). Node 04 pre-validates the ligand with RDKit before invoking Uni-Mol — check the error output for the SMILES that failed parsing.

**Node 02 pocket scores are all low or pocket center is far from the expected site**
Confirm that P2Rank is running with `-c alphafold`. Without this flag, pLDDT values in the B-factor column are misread as thermal displacement factors, degrading pocket scoring.

**Node 01 runs out of VRAM**
Full EGFR sequence (~1186 residues) requires ≥ 16 GB VRAM. Set `PARAM_ACCELERATOR=cpu` for a slower but memory-unconstrained run, or truncate the input to the kinase domain only.

## References

- Wohlwend et al. "Boltz-1: Democratizing Biomolecular Interaction Modeling." *bioRxiv*, 2024.
- Passaro et al. "Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction." *bioRxiv*, 2025.
- Jakubec et al. "PrankWeb 4: Neural Networks, Evolutionary Information, and AlphaFold Structures." *Nucleic Acids Research*, 2025.
- Zhou et al. "Uni-Mol: A Universal 3D Molecular Representation Learning Framework." *ICLR*, 2023. https://openreview.net/forum?id=6K2RM6wVqKu
- Eguida & Rognan. "Estimating the Ease of Protein–Ligand Docking with Predicted Structures." *JCIM*, 2023. PMC9852548.
- [Boltz GitHub](https://github.com/jwohlwend/boltz)
- [P2Rank GitHub](https://github.com/rdk/p2rank)
- [Uni-Mol GitHub](https://github.com/deepmodeling/Uni-Mol)
