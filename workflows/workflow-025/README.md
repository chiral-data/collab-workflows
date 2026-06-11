---
doc_id: workflow-025
domain: molecular-docking
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  Compares AutoDock Vina and GNINA docking performance on Carbonic Anhydrase II (PDB: 1OKL),
  a Zn²⁺-binding target, using the same MCMC pose generation engine to isolate
  scoring function differences between empirical and CNN-based methods.
tags: [molecular-docking, vina, gnina, virtual-screening, carbonic-anhydrase, metal-binding]
---

# Workflow 025: Docking Comparison of AutoDock Vina vs GNINA

This workflow compares two molecular docking tools — **AutoDock Vina** and **GNINA** — on a metal-coordinating target, **Carbonic Anhydrase II (PDB: 1OKL)**, which binds Zn²⁺. Both tools share the same MCMC sampling engine (GNINA is a fork of Smina, which forks Vina), so any observed differences in ranking, enrichment, or predicted affinity are attributable solely to the **scoring function**: Vina's empirical physics-based energy approximation vs. GNINA's 3D CNN.

## Overview

Vina's additive atom-contact terms cannot model coordinate covalent bonds to transition metal ions (Zn²⁺), which skews ΔG predictions and can overscore large, heavy compounds with more pocket contacts. GNINA's CNN implicitly learns the geometric "picture" of a metal coordination pocket but can occasionally score chemically implausible poses highly due to training data bias.

Buccheri et al. (2025) found that for CA-II, GNINA achieved a pose RMSD of 1.37 Å vs. Vina's 6.78 Å, and an enrichment factor (EF1%) of 20.75 vs. 0 — making the comparison scientifically clean and interpretable.

The workflow validates inputs, prepares receptor and ligand structures, runs a redocking QC check to confirm the binding pocket, runs parallel Vina and GNINA screens, and generates a side-by-side HTML report with score correlations, ranking histograms, and runtime metrics.

## When to use this workflow

Use this workflow when:
- You are screening a ligand library against a **metal-coordinating target** and want to compare empirical vs. CNN-based scoring.
- You want a reproducible benchmark run on CA-II (PDB: 1OKL) before applying to a novel target.
- You need a self-contained HTML report comparing Vina and GNINA performance metrics.

**When NOT to use:**
- For non-metal targets without coordination chemistry: standard Vina screening is sufficient (see `workflow-014`).
- For GNINA-only runs without Vina comparison: see `workflow-018`.

## Architecture and data flow

```
input_files/ligands.smiles
        │
        ▼
01_validate_inputs ──► validated receptor + ligand files
        │
        ▼
02_prepare_structures ──► receptor.pdbqt, optimized_screening_library.pdbqt, pocket_config.txt
        │
        ▼
03_target_redocking ──► qc_validation_results.json (CNN score + RMSD pass/fail)
        │
        ├──────────────────────┐
        ▼                      ▼
04a_run_vina           04b_run_gnina
vina_screening_poses   gnina_screening_poses
        │                      │
        └──────────┬───────────┘
                   ▼
         05_generate_report
         comparative_screening_metrics.csv
         docking_performance_report.html
```

## Input requirements

- **`input_files/ligands.smiles`**: SMILES file with one compound per line; identifier in the second column.
- Receptor is fetched automatically from RCSB PDB using PDB ID `1OKL`; no manual download required.
- Test inputs are provided at `01_validate_inputs/.chiral/test_inputs/sample_ligands.smiles`.

## Workflow nodes

### Node 01: Validate Inputs

**Goal:** Confirm all required inputs are present and fetch the CA-II receptor from RCSB PDB.

**Process:** Checks that the ligand SMILES file is present and well-formed, queries RCSB PDB for structure `1OKL`, validates resolution against a threshold, and writes the receptor PDB file for downstream processing.

**Scientific notes:** Resolution filtering ensures the receptor structure is suitable for docking; CA-II at 1OKL resolves to ≤2 Å, which is considered high quality for docking.

**Outputs:**
- Validated receptor PDB file
- Validated ligand SMILES file

### Node 02: Prepare Structures

**Goal:** Convert receptor to PDBQT and generate 3D-optimized ligand conformers.

**Process:** Uses Open Babel to convert the receptor PDB to PDBQT format and generate initial 3D ligand geometries from SMILES, then runs xTB for semi-empirical geometry optimization. A pocket configuration file is derived from the receptor's zinc coordination center.

**Scientific notes:** xTB optimization at the GFN2-xTB level corrects strained geometries that Open Babel's rule-based 3D builder can introduce, producing more physically realistic starting conformations for docking.

**Outputs:**
- `receptor.pdbqt`
- `optimized_screening_library.pdbqt`
- `pocket_config.txt`

### Node 03: Target Redocking (QC)

**Goal:** Confirm the binding pocket is well-defined by redocking the co-crystallized ligand.

**Process:** Runs GNINA redocking of the native ligand against the prepared receptor. Parses the CNN score and pose RMSD from GNINA output and writes a pass/fail JSON result. Downstream nodes depend on this QC gate.

**Scientific notes:** A redocking RMSD < 2 Å indicates the receptor preparation and pocket definition are suitable for virtual screening. This step catches preparation errors before an expensive full screen.

**Outputs:**
- `qc_validation_results.json` (CNN score + RMSD pass/fail)

### Node 04a: Run Vina

**Goal:** Screen the ligand library using AutoDock Vina.

**Process:** Invokes the Vina Python API with the prepared receptor PDBQT and pocket configuration, writing docked poses and a runtime log.

**Scientific notes:** Vina uses an empirical scoring function; scores are in kcal/mol (more negative = more favorable). Metal coordination interactions are not explicitly modeled.

**Outputs:**
- `vina_screening_poses.pdbqt`
- `vina_runtime_log.txt`

### Node 04b: Run GNINA

**Goal:** Screen the ligand library using GNINA's CNN scoring function.

**Process:** Invokes the GNINA prebuilt binary with the same pocket configuration as node 04a. GNINA rescores each pose with a 3D CNN trained on PDBbind data.

**Scientific notes:** GNINA's CNN score captures geometric features of metal coordination pockets not encoded in additive contact terms. CNN scores > 0.5 are considered confident binders.

**Outputs:**
- `gnina_screening_poses.pdbqt`
- `gnina_runtime_log.txt`

### Node 05: Generate Report

**Goal:** Produce a side-by-side HTML report comparing Vina and GNINA performance.

**Process:** Parses output from nodes 04a and 04b, computes score correlations, ranking histograms, pose RMSD per compound, and runtime metrics, then renders a self-contained HTML report.

**Outputs:**
- `comparative_screening_metrics.csv`
- `docking_performance_report.html`

## Parameters

### `pdb_id`

| Value / Range | Description |
|---------------|-------------|
| `1OKL` (default) | CA-II structure with Zn²⁺; used for all benchmark runs. |
| Any 4-character PDB ID | Alternative target; pocket config must be re-derived manually. |

**Trade-off:** Using a non-default PDB ID requires verifying the pocket definition; automated pocket detection is only validated for `1OKL`.

### `resolution_threshold`

| Value / Range | Description |
|---------------|-------------|
| `2.0` Å (default) | Rejects receptor structures with resolution > 2 Å. |
| `2.5` Å | More permissive; acceptable for exploratory runs. |

**Trade-off:** Stricter thresholds improve docking reliability but may exclude valid structures.

## Outputs and interpretation

### `docking_performance_report.html`

Self-contained HTML report with score correlation scatter plots, ranking histograms, per-compound pose RMSD comparisons, and runtime metrics. Can be opened in any browser without re-running code.

### `comparative_screening_metrics.csv`

Per-compound Vina ΔG (kcal/mol), GNINA CNN score, Vina rank, GNINA rank, and pose RMSD. Use for downstream statistical analysis or filtering.

### `qc_validation_results.json`

Redocking QC result containing `rmsd`, `cnn_score`, and `pass` (boolean). If `pass` is `false`, inspect receptor preparation in node 02 before proceeding.

## Quick start

### Running on Silva

1. Select workflow **025** from the workflow list.
2. Upload `input_files/ligands.smiles`.
3. Adjust parameters (PDB ID, resolution threshold) if needed.
4. Click **Run**.

### Running with Docker

- AutoDock Vina Dockerfile: `workflows/workflow-004/Dockerfile`
- GNINA Dockerfile: `apps/g/gnina_2025_12_04/Dockerfile`

### Test vs production settings

| Setting | Test (default) | Production |
|---------|----------------|------------|
| Ligand library | `sample_ligands.smiles` (5 compounds) | Full screening library |
| `resolution_threshold` | `2.0` Å | `2.0` Å |

## Troubleshooting

**GNINA score parsing returns None**
GNINA output format may vary by binary version. Ensure you are using the binary from `apps/g/gnina_2025_12_04/Dockerfile`. Check `gnina_runtime_log.txt` for raw output.

**Node 03 QC fails (RMSD > 2 Å)**
Likely caused by incorrect pocket center coordinates in `pocket_config.txt`. Re-examine the zinc coordination site in the prepared receptor PDBQT and adjust the box center accordingly.

**xTB optimization fails for a ligand**
Some SMILES strings produce geometries that xTB cannot optimize. Check Open Babel's 3D output for that ligand and consider providing a pre-optimized SDF instead.

## References

- Buccheri et al. "Benchmarking GNINA and AutoDock Vina for metal-coordinating targets." *PMC*, 2025. https://pmc.ncbi.nlm.nih.gov/articles/PMC12388557/
- [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina) — Apache 2.0; `pip install vina`
- [GNINA](https://github.com/gnina/gnina) — Apache 2.0; prebuilt binary
- [Open Babel](https://github.com/openbabel/openbabel) — GPL 2.0; `conda install -c conda-forge openbabel`
- [xTB](https://github.com/grimme-lab/xtb) — LGPL 3.0; `conda install -c conda-forge xtb`
