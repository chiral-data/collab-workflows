# Workflow 017: LightDock Protein-Protein Docking

Structure-based protein-protein docking using [LightDock](https://github.com/lightdock/lightdock). This workflow predicts the 3D structure of a protein complex from two input PDB files, using swarm-based optimization and scoring. It is suitable for modeling binary protein interactions when both partners have known structures.

## Overview

LightDock is a flexible, swarm-based docking framework for predicting the structure of protein-protein complexes. It uses glowworm swarm optimization to efficiently sample and rank possible binding modes. Use this workflow when you have two protein structures (in PDB format) and want to predict their likely binding interface and top-scoring complex models.

## Input Requirements

- **Two protein structures in PDB format:**
    - `protein1.pdb` (ligand)
    - `protein2.pdb` (receptor)

Sample files for the barnase-barstar protein complex are included in `input_files/` for testing and demonstration.

## Nodes

| Node | Name                   | Description                                                        |
|------|------------------------|--------------------------------------------------------------------|
| 01   | Validate Inputs        | Checks PDB files for completeness, CA atoms, and chain integrity   |
| 02   | Prepare Docking Swarms | Generates initial swarm positions and glowworm configurations      |
| 03   | Run LightDock Docking  | Runs LightDock simulation, clusters, and ranks docking poses       |
| 04   | Generate Docking Report| Builds a self-contained HTML report with 3D visualization          |

Nodes run sequentially: 01 → 02 → 03 → 04.

## Parameters

| Parameter         | Node | Type    | Default | Description                                                        |
|-------------------|------|---------|---------|--------------------------------------------------------------------|
| receptor_file     | 01   | string  | inputs/protein2_barnase.pdb | Path to receptor PDB file                |
| ligand_file       | 01   | string  | inputs/protein1_barstar.pdb | Path to ligand PDB file                  |
| num_swarms        | 02   | integer | 10      | Number of independent swarms (400 for production, 5 for testing)   |
| num_glowworms     | 02   | integer | 200     | Number of glowworms per swarm                                      |
| steps             | 03   | integer | 50      | Number of GSO optimization steps (50 for testing, 200+ for prod)   |
| scoring_function  | 03   | enum    | fastdfire | LightDock scoring function (fastdfire, dfire, dfire2, pydock)   |
| num_conformations | 03   | integer | 200     | Number of conformations to generate per swarm                      |
| top_n             | 03   | integer | 10      | Number of top-ranked poses to collect into output                  |

## Output Files

- **receptor.pdb, ligand.pdb** — Validated input structures (copied to working dir)
- **validation_report.json** — Input validation summary (atom/residue counts, chains)
- **lightdock_workspace.tar.gz** — Swarm and glowworm configuration archive
- **top_predictions.pdb** — Top N predicted complex structures (multi-model PDB)
- **rank_results.json** — Docking scores and ranking for all poses
- **report.html** — Interactive HTML report with 3D viewer and ranked pose table

Sample output files are provided in `results/` and `sample_outputs/`.

## Example

To test the workflow, use the provided barnase-barstar PDB files in `input_files/`. The workflow will validate the structures, generate swarms, run docking, and produce a ranked list of predicted complexes with an interactive report.

---
For more details on LightDock, see the [LightDock documentation](https://lightdock.org/). For workflow internals, see the scripts in each node directory.