---
doc_id: workflow-017
domain: protein-protein-docking
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  Protein-protein docking using LightDock's glowworm swarm optimization.
  Takes two PDB structures and predicts the binding interface and top-scoring
  complex models with an interactive HTML report.
tags: [lightdock, protein-protein-docking, gso, dfire, swarm-optimization]
---

# Workflow 017: LightDock Protein-Protein Docking

Structure-based protein-protein docking using [LightDock](https://github.com/lightdock/lightdock). This workflow predicts the 3D structure of a protein complex from two input PDB files using glowworm swarm optimization (GSO) and produces a ranked list of predicted complex structures with an interactive HTML report.

## Overview

LightDock is a flexible docking framework that uses the Glowworm Swarm Optimization (GSO) algorithm to sample and rank possible binding modes between two protein structures. Multiple independent "swarms" of "glowworms" explore the conformational space in parallel, each swarm starting from a different position on the receptor surface. After a configurable number of optimization steps, conformations are generated, clustered using BSAS (Basic Sequential Algorithmic Scheme), and ranked by the chosen scoring function. The default scoring function is fastDFIRE, an efficient implementation of the DFIRE (Distance-scaled, Finite Ideal-gas Reference) statistical potential.

The workflow includes a download node that fetches the barnase-barstar complex (PDB 1BRS) as a built-in test case. For production use, provide your own receptor and ligand PDB files (Jiménez-García et al., 2018).

## When to use this workflow

Use this workflow when you have two protein structures in PDB format and want to predict their binding interface and top-scoring complex models. It is suited for binary protein-protein interactions where both partners have known experimental structures. The included barnase-barstar test case (PDB 1BRS, chains A and D) provides a well-characterized benchmark for validation.

Do not use this workflow for protein-ligand (small molecule) docking — use workflow-004 (AutoDock Vina) or workflow-002/003 (Smina). For protein-protein docking with deep learning scoring, use workflow-016 (DiffDock-PP). For predicting protein complex structure from sequence when experimental structures are unavailable, use workflow-012 (Boltz-2). LightDock does not model large conformational changes — both input structures should be close to their bound-state conformations.

## Architecture and data flow

```text
[00: Download] ──> [01: Validate Inputs] ──> [02: Prepare Swarms] ──> [03: Run LightDock] ──> [04: Generate Report]
       |                    |                        |                        |                        |
  protein1.pdb        receptor.pdb          lightdock_workspace       top_predictions.pdb        report.html
  protein2.pdb        ligand.pdb               .tar.gz               rank_results.json
                  validation_report.json
```

Nodes run sequentially: 00 → 01 → 02 → 03 → 04.

## Input requirements

- **Two protein structures in PDB format:**
  - Receptor PDB file (default: `protein2_barnase.pdb`)
  - Ligand PDB file (default: `protein1_barstar.pdb`)
- **Requirements:** Each PDB must contain ATOM records with valid coordinates and at least one CA (alpha-carbon) atom per chain. HETATM records are accepted.
- **Sample data:** Node 00 downloads PDB 1BRS from RCSB and extracts chain A (barnase) and chain D (barstar) as the test inputs. To use your own structures, place them in `input_files/` and set the `receptor_file` and `ligand_file` parameters.

## Workflow nodes

### Node 00: Download Sample Inputs

**Goal:** Fetch the barnase-barstar test structures from RCSB PDB.

**Process:** Downloads PDB 1BRS, then extracts chain A (barnase, saved as `protein2_barnase.pdb`) and chain D (barstar, saved as `protein1_barstar.pdb`) by filtering ATOM records on the chain ID column. Only ATOM/HETATM and END lines are retained.

**Scientific notes:** Barnase-barstar is a well-characterized protein-protein interaction commonly used as a docking benchmark. The complex has a high binding affinity (Kd ~ 10⁻¹⁴ M) with a well-defined interface.

**Outputs:**
- `protein2_barnase.pdb` — barnase (chain A from 1BRS)
- `protein1_barstar.pdb` — barstar (chain D from 1BRS)

### Node 01: Validate Inputs

**Goal:** Check that both PDB files are valid protein structures.

**Process:** Parses each PDB file to count ATOM/HETATM records, unique residues, and chains. Verifies that at least one CA atom is present (confirming the file contains protein). Copies validated files to standardized names (`receptor.pdb`, `ligand.pdb`) and writes a validation summary JSON.

**Scientific notes:** The CA atom check ensures the input is a protein structure rather than a small molecule or empty file. No structural cleanup (hydrogen addition, missing residue modeling) is performed — LightDock's setup handles that.

**Outputs:**
- `receptor.pdb` — validated receptor
- `ligand.pdb` — validated ligand
- `validation_report.json` — atom/residue/chain counts and validation status

### Node 02: Prepare Docking Swarms

**Goal:** Generate initial swarm positions and glowworm configurations for LightDock.

**Process:** Runs `lightdock3_setup.py` with the validated receptor and ligand PDBs. Uses flags `--noxt` (no extra terminal oxygen atoms), `--noh` (no hydrogen atoms), and `--now` (no water molecules). Creates the configured number of swarm directories, each containing initial glowworm positions. Archives the entire workspace (setup.json, init/, swarm directories, lightdock_ files) into `lightdock_workspace.tar.gz` for transfer to Node 03.

**Scientific notes:** Each swarm is placed at a different position on the receptor surface, providing diverse starting orientations for the optimization. The number of swarms controls spatial coverage — 10 swarms is suitable for testing, while 400 swarms provides thorough coverage for production. Each glowworm within a swarm represents a candidate binding pose defined by translation and rotation parameters.

**Outputs:**
- `receptor.pdb`, `ligand.pdb` — copies for downstream use
- `lightdock_workspace.tar.gz` — archived swarm configurations

### Node 03: Run LightDock Docking

**Goal:** Run the GSO simulation, generate conformations, cluster poses, and rank results.

**Process:** Extracts the workspace archive, then executes the LightDock pipeline:
1. `lightdock3.py` — runs GSO optimization for the configured number of steps with the chosen scoring function
2. `lgd_generate_conformations.py` — generates PDB conformations from the optimized glowworm positions (per swarm)
3. `lgd_cluster_bsas.py` — clusters similar poses within each swarm using BSAS
4. `lgd_rank.py` — ranks all poses across swarms by scoring function value

Parses `rank_by_scoring.list` to extract the top N poses (falls back to per-swarm best gso scores if the rank file is unavailable). Collects top poses into a multi-model PDB file.

**Scientific notes:** GSO is a swarm intelligence algorithm where glowworms move toward brighter neighbors (higher-scoring poses). The fastDFIRE scoring function is a statistical potential derived from known protein-protein interfaces — higher scores indicate more favorable interactions. BSAS clustering removes near-duplicate poses using an RMSD threshold. The number of optimization steps controls convergence — 50 steps is suitable for quick testing, while 200+ steps allows more thorough optimization.

**Outputs:**
- `top_predictions.pdb` — multi-model PDB with top N predicted complexes (each model includes a REMARK with rank and score)
- `rank_results.json` — docking parameters and ranked pose list with scores

### Node 04: Generate Docking Report

**Goal:** Build an interactive HTML report with 3D visualization and ranked pose table.

**Process:** Reads `rank_results.json` and `top_predictions.pdb`, generates a self-contained HTML report using Mol* (Molstar) for 3D molecular visualization. The report includes a ranked table of poses with scores and an interactive 3D viewer showing the predicted complex structures.

**Scientific notes:** The report visualizes the top-ranked complex model by default. Examining multiple top poses can reveal whether the predicted interface is consistent across independent swarms, which increases confidence in the prediction.

**Outputs:**
- `report.html` — self-contained HTML report with 3D viewer and pose ranking table

## Parameters

### receptor_file

- **Type:** string
- **Default:** `"inputs/protein2_barnase.pdb"`
- **Node:** 01
- **Description:** Path to the receptor PDB file.

### ligand_file

- **Type:** string
- **Default:** `"inputs/protein1_barstar.pdb"`
- **Node:** 01
- **Description:** Path to the ligand (second protein) PDB file.

### num_swarms

- **Type:** integer
- **Default:** `10`
- **Node:** 02
- **Description:** Number of independent swarms placed on the receptor surface.

| Value | Use case |
|-------|----------|
| `5`–`10` | Quick testing |
| `10` (default) | Default test run |
| `400` | Production — thorough surface coverage |

**Trade-off:** More swarms provide better spatial coverage but increase runtime linearly.

### num_glowworms

- **Type:** integer
- **Default:** `200`
- **Node:** 02
- **Description:** Number of glowworms (candidate poses) per swarm.
- **Guidance:** 200 is the standard default. Reducing to 50–100 speeds up testing; increasing beyond 200 has diminishing returns.

### steps

- **Type:** integer
- **Default:** `50`
- **Node:** 03
- **Description:** Number of GSO optimization steps.

| Value | Use case |
|-------|----------|
| `50` (default) | Testing |
| `100` | Moderate quality |
| `200`+ | Production — thorough optimization |

**Trade-off:** More steps allow better convergence at the cost of proportional runtime increase.

### scoring_function

- **Type:** enum
- **Default:** `"fastdfire"`
- **Node:** 03
- **Description:** Scoring function used to evaluate binding poses.

| Value | Description |
|-------|-------------|
| `fastdfire` (default) | Fast implementation of the DFIRE statistical potential. Good general-purpose choice. |
| `dfire` | Original DFIRE implementation. Slower but identical scoring. |
| `dfire2` | Updated DFIRE parameterization. |
| `pydock` | PyDock scoring function combining electrostatics, desolvation, and van der Waals. |

### num_conformations

- **Type:** integer
- **Default:** `200`
- **Node:** 03
- **Description:** Number of conformations to generate per swarm for clustering and ranking.

### top_n

- **Type:** integer
- **Default:** `10`
- **Node:** 03
- **Description:** Number of top-ranked poses to collect into `top_predictions.pdb`.

## Outputs and interpretation

### Scoring function values

LightDock scores are **not** binding free energies in kcal/mol. They are unitless scoring function values where **higher (more positive) = better predicted binding**. This is the opposite convention from Vina/Smina-based workflows.

For fastDFIRE/DFIRE scoring, typical ranges for well-docked complexes are in the hundreds (e.g., 300–600). Values vary significantly by system size and interface area, so absolute thresholds are not meaningful — compare scores within the same docking run.

### top_predictions.pdb

A multi-model PDB file containing the top N predicted complex structures. Each MODEL block includes a REMARK line with the rank and score. Load in PyMOL, ChimeraX, or the included HTML report to inspect individual poses.

### rank_results.json

Machine-readable JSON with docking parameters (num_swarms, steps, scoring_function) and the ranked pose list (rank and score for each).

## Quick start

### Running with Docker

All nodes use the same container image:

```bash
docker build -t lightdock-pipeline:latest .
```

### Running on Silva

1. Select "LightDock Protein-Protein Docking" from the workflow list
2. Upload your two protein PDB files (or use the built-in barnase-barstar test)
3. Adjust `num_swarms` and `steps` for your desired quality level
4. Click Run

### Test vs production settings

| Setting | Test (default) | Production |
|---------|---------------|------------|
| `num_swarms` | `10` | `400` |
| `num_glowworms` | `200` | `200` |
| `steps` | `50` | `200`+ |
| `scoring_function` | `fastdfire` | `fastdfire` |
| `top_n` | `10` | `10`–`50` |

A successful test run with the barnase-barstar defaults completes in a few minutes and produces 10 ranked complex models.

## References

- Jiménez-García, B., Roel-Touris, J., Romero-Durana, M., Vidal, M., Jiménez-González, D. & Fernández-Recio, J. "LightDock: a new multi-scale approach to protein-protein docking." *Bioinformatics* 34(1):49–55, 2018. DOI: https://doi.org/10.1093/bioinformatics/btx555
- Zhou, H. & Zhou, Y. "Distance-scaled, finite ideal-gas reference state improves structure-derived potentials of mean force for structure selection and stability prediction." *Protein Sci.* 11(11):2714–2726, 2002. DOI: https://doi.org/10.1110/ps.0217002
- [LightDock documentation](https://lightdock.org/)
- [LightDock GitHub](https://github.com/lightdock/lightdock)
