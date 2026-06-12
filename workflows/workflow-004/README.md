---
doc_id: workflow-004
domain: molecular-docking
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  End-to-end protein-ligand virtual screening using AutoDock Vina with
  P2Rank pocket prediction. Takes a PDB ID and PubChem compound CIDs as
  input and produces ranked binding affinities with an interactive dashboard.
tags: [autodock-vina, virtual-screening, p2rank, protein-ligand-docking, pocket-prediction]
---

# Workflow 004: AutoDock Vina Virtual Screening

A fully modular pipeline for protein-ligand virtual screening using AutoDock Vina with integrated P2Rank pocket prediction. The workflow downloads a receptor from RCSB PDB and ligands from PubChem, prepares both for docking, predicts binding pockets, and runs high-throughput docking with interactive visualization of results.

## Overview

AutoDock Vina is a widely used molecular docking program that predicts the binding affinity and pose of small-molecule ligands to a protein target. Vina uses an empirical scoring function combining steric (Gaussian and repulsion), hydrophobic, and hydrogen-bonding terms, optimized via iterated local search with BFGS quasi-Newton optimization. This workflow wraps Vina in a six-node pipeline that automates receptor and ligand acquisition, structure preparation (PDBFixer, OpenBabel), ML-based binding pocket prediction (P2Rank), and docking with ranked output. Each node generates an interactive HTML report with 3D molecular visualization using NGL.js.

The pipeline is designed for structure-based virtual screening where the researcher has a known protein target (PDB ID) and a set of candidate compounds (PubChem CIDs) and wants to rank them by predicted binding affinity (Eberhardt et al., 2021; Trott & Olson, 2010; Krivák & Hoksza, 2018).

## When to use this workflow

Use this workflow when you have a protein target with a known PDB structure and a set of candidate small-molecule ligands identified by PubChem CID. It is best suited for screening a small to moderate number of compounds (up to ~50) with automated pocket detection — you do not need to know the binding site in advance. The workflow handles the full preparation pipeline from raw PDB/SDF files to docking-ready PDBQT format.

Do not use this workflow for screening large compound libraries from external sources such as Zenodo — use workflow-002 (Smina library screening) instead. If you want to generate and screen chemical variants of a known co-crystallized ligand, use workflow-003 (Smina ligand variant screening). For protein-protein docking, use workflow-016 (DiffDock-PP) or workflow-017 (LightDock). If you need ADMET property predictions rather than binding affinity, use workflow-014.

## Architecture and data flow

```text
┌─────────────────────┐
│ 01: Receptor        │
│    Acquisition      │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│ 03: Receptor        │
│   Preparation       │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐       ┌─────────────────────┐
│ 05: Pocket          │       │ 02: Ligand          │
│   Discovery         │       │   Collection        │
└──────────┬──────────┘       └──────────┬──────────┘
           │                             │
           │                             ▼
           │                  ┌─────────────────────┐
           │                  │ 04: Ligand          │
           │                  │   Preparation       │
           │                  └──────────┬──────────┘
           │                             │
           └──────────┬──────────────────┘
                      ▼
           ┌─────────────────────┐
           │ 06: VS Analytics    │
           │  (Vina + Dashboard) │
           └─────────────────────┘
```

Two independent branches (receptor: 01 → 03 → 05, ligand: 02 → 04) converge at Node 06 for docking. Nodes 01 and 02 can run in parallel.

## Input requirements

- **Receptor:** A valid PDB ID (e.g., `5kir`, `1hsg`). The structure is downloaded from RCSB PDB automatically.
- **Ligands:** A JSON array of PubChem Compound IDs (CIDs) as strings (e.g., `["3672", "2662"]`). 3D conformers are downloaded from PubChem.
- **No local files required** — all inputs are fetched by the first two nodes. The `input_files/` directory is reserved for future use.

## Workflow nodes

### Node 01: Receptor Acquisition

**Goal:** Download and validate a receptor structure from RCSB PDB.

**Process:** Uses BioPython's `PDBList` to fetch the PDB file for the given `pdb_id` parameter, renames it to `{pdb_id}.pdb`, then parses the structure with `PDBParser` to count atoms, residues, and chains. Generates an interactive HTML report with an NGL.js 3D viewer.

**Scientific notes:** The downloaded PDB file contains the deposited crystal structure including water molecules, heterogens (ligands, cofactors), and all chains. These are retained at this stage and cleaned in Node 03. Structure quality depends on the experimental resolution of the deposited structure.

**Outputs:**
- `{pdb_id}.pdb` — raw crystal structure
- `receptor_metadata.json` — atom/residue/chain counts and download timestamp
- `report.html` — interactive 3D structure viewer

### Node 02: Ligand Collection

**Goal:** Download ligand 3D conformers from PubChem.

**Process:** For each CID in the `ligand_ids` parameter, fetches the 3D SDF conformer from PubChem's PUG REST API (`/rest/pug/compound/CID/{cid}/SDF?record_type=3d`). Validates download status and file sizes. Generates a multi-ligand HTML report.

**Scientific notes:** PubChem provides pre-generated 3D conformers. These are reasonable starting geometries but will be further optimized during docking. If a CID has no 3D conformer available, the download will fail for that compound and the workflow will exit if zero ligands succeed.

**Outputs:**
- `{cid}.sdf` — 3D conformer for each ligand
- `ligands_metadata.json` — download status per compound
- `report.html` — grid of ligand 3D viewers

### Node 03: Receptor Preparation

**Goal:** Clean the receptor structure and convert to docking-ready PDBQT format.

**Process:** Applies PDBFixer to find and replace non-standard residues, add missing atoms, and remove heterogens (while keeping crystallographic waters via `removeHeterogens(False)`). Adds hydrogens at pH 7.0 using OpenMM. Converts the cleaned PDB to PDBQT format using OpenBabel with `-xr` (rigid molecule, no torsion tree) and `-xh` (preserve explicit hydrogens) flags.

**Scientific notes:** The `-xr` flag produces a rigid receptor PDBQT without rotatable bonds, which is standard for Vina receptor files. Protonation at pH 7.0 assumes physiological conditions — adjust the source code if your target operates at a different pH (e.g., lysosomal proteins at pH 4.5–5.0). Removing heterogens strips co-crystallized ligands and cofactors, which is appropriate for blind docking but may lose important context for allosteric sites.

**Outputs:**
- `protein_fixed.pdb` — cleaned PDB structure
- `receptor.pdbqt` — docking-ready rigid receptor
- `refined_receptor_metadata.json` — preparation steps and file size statistics
- `report.html` — preparation pipeline visualization

### Node 04: Ligand Preparation

**Goal:** Convert ligands from SDF to docking-ready PDBQT format.

**Process:** Converts each SDF file to PDBQT using OpenBabel with the `-xh` flag (preserve explicit hydrogens). Counts atoms and detects the number of rotatable bonds from the TORSDOF field in the output PDBQT.

**Scientific notes:** Unlike the receptor (`-xr`), ligand PDBQT files include a torsion tree that defines which bonds are sampled during docking. Vina handles up to ~10–12 rotatable bonds well; accuracy degrades for highly flexible ligands with more. The TORSDOF (torsional degrees of freedom) value in the PDBQT affects the entropy penalty in the scoring function.

**Outputs:**
- `{ligand}.pdbqt` — docking-ready ligand files with torsion tree
- `refined_ligands_metadata.json` — conversion statistics and rotatable bond counts per ligand
- `report.html` — multi-ligand preparation report

### Node 05: Pocket Discovery

**Goal:** Predict binding pockets on the receptor and define the docking grid box.

**Process:** Runs P2Rank (v2.4.2) on the cleaned receptor PDB. Parses the `_predictions.csv` output to extract pocket rankings (score, probability, center coordinates, residue count, surface atoms) and the `_residues.csv` for per-pocket residue lists. Automatically places a fixed 20 × 20 × 20 Å grid box centered on the top-ranked pocket and writes `grid_config.json`.

**Scientific notes:** P2Rank is a machine-learning pocket predictor that uses a random forest classifier trained on geometric and physicochemical features of the protein's solvent-accessible surface. It evaluates individual surface points and clusters them into pockets. The "probability" output is a calibrated probability of the pocket being a genuine ligand-binding site. The grid box size is fixed at 20 Å per side — this is suitable for typical drug-like binding sites but may be too small for large, shallow binding grooves or too large for very compact pockets. Verify the selected pocket is biologically relevant; P2Rank may rank a crystallographic artifact above the true active site.

**Outputs:**
- `pockets.pdb` — predicted pocket locations
- `protein.pdb` — receptor copy for visualization
- `grid_config.json` — grid box center (x, y, z) and size (20 × 20 × 20 Å)
- `pocket_discovery_metadata.json` — ranked pocket list with scores, coordinates, and residues
- `report.html` — interactive pocket visualization with grid box overlay

### Node 06: VS Analytics

**Goal:** Run AutoDock Vina docking for all ligands and generate a ranked results dashboard.

**Process:** Identifies the receptor as the largest PDBQT file (by file size) and remaining PDBQT files as ligands. Loads the grid configuration from Node 05. Docks each ligand by invoking `vina` with the receptor, grid center/size, and user-configurable exhaustiveness/num_modes/energy_range parameters. Splits multi-model output PDBQT files into individual pose files, extracts binding affinities from `REMARK VINA RESULT` lines, ranks ligands by best affinity, and generates an interactive HTML dashboard with 3D visualization and distance-based interaction analysis (contacts < 4.0 Å).

**Scientific notes:** Vina's search algorithm uses iterated local search with BFGS quasi-Newton optimization — it generates random initial conformations, applies gradient-based local optimization, then perturbs and repeats. The `exhaustiveness` parameter controls how many independent runs are performed; higher values reduce the chance of missing the global minimum but increase runtime linearly. The scoring function has a reported standard error of ~2.85 kcal/mol against experimental binding data (Trott & Olson, 2010), so differences smaller than ~1 kcal/mol between compounds should not be considered meaningful. GPU acceleration is auto-detected if the installed Vina binary supports the `--gpu` flag.

**Outputs:**
- `{ligand}_docked.pdbqt` — multi-pose docking results
- `{ligand}_pose{N}.pdbqt` — individual extracted poses
- `{ligand}_docking.log` — Vina scoring output
- `receptor.pdbqt` — receptor copy for reference
- `virtual_screening_metadata.json` — ranked results summary with affinities and parameters
- `report.html` — interactive dashboard with 3D viewer, affinity ranking, and interaction analysis

## Parameters

### pdb_id

- **Type:** string
- **Default:** `"5kir"`
- **Node:** 01
- **Description:** RCSB PDB identifier for the target receptor.
- **Guidance:** Use a structure with good resolution (< 2.5 Å) and minimal missing residues in the binding site region.

### ligand_ids

- **Type:** string (JSON array)
- **Default:** `'["3672", "2662"]'`
- **Node:** 02
- **Description:** PubChem Compound IDs to download and dock.
- **Guidance:** Use CIDs for compounds with available 3D conformers. The default pair (Sunitinib CID 3672 and Imatinib CID 2662) are known kinase inhibitors suitable for testing with the 5KIR target.

### exhaustiveness

- **Type:** integer
- **Default:** `32`
- **Node:** 06
- **Description:** Number of independent Vina search runs. Each run starts from a random initial conformation and performs iterated local search.

| Value | Use case |
|-------|----------|
| `8` | Quick test (~1 min per ligand) |
| `32` (default) | Standard screening |
| `64`–`128` | Publication-quality results |

**Trade-off:** Higher exhaustiveness improves the probability of finding the global minimum pose at the cost of linear runtime increase.

**Test vs production:** The default of 32 is suitable for initial screening. For publication results, use 64 or higher.

### num_modes

- **Type:** integer
- **Default:** `10`
- **Node:** 06
- **Description:** Maximum number of binding poses to generate per ligand.
- **Guidance:** 10 is sufficient for most screening purposes. Increase to 20 if you need to analyze alternative binding modes or if the pocket is large and shallow.

### energy_range

- **Type:** float
- **Default:** `5.0`
- **Node:** 06
- **Description:** Energy window in kcal/mol — only poses within this range of the best pose are retained.
- **Guidance:** 5.0 kcal/mol is a reasonable default. Narrowing to 3.0 filters out weak alternative poses; widening to 10.0 captures more diverse binding modes.

## Outputs and interpretation

### Binding affinity (kcal/mol)

The primary output is the predicted binding free energy for each ligand, extracted from Vina's `REMARK VINA RESULT` lines. More negative values indicate stronger predicted binding.

| Range | Interpretation |
|-------|---------------|
| < −10 kcal/mol | Very strong binding (rare, verify carefully) |
| −8 to −10 kcal/mol | Strong binding — promising hit |
| −6 to −8 kcal/mol | Moderate binding — worth investigating |
| > −6 kcal/mol | Weak binding — likely not a viable candidate |

Vina's scoring function has a standard error of ~2.85 kcal/mol against experimental data, so treat rankings as relative rather than absolute. Differences of < 1 kcal/mol between compounds are not statistically meaningful.

### P2Rank pocket scores

Each predicted pocket has two key metrics:
- **Score:** Raw ligandability score from the random forest classifier (higher = more likely binding site).
- **Probability:** Calibrated probability (0–1) of being a genuine ligand-binding site. Values > 0.5 generally indicate druggable pockets.

The top-ranked pocket is automatically used for grid box placement.

### Interaction analysis

The dashboard reports distance-based contacts (< 4.0 Å) between the docked ligand and receptor residues. These are purely geometric measurements, not verified hydrogen bonds or hydrophobic interactions — use them as a starting point for detailed interaction analysis in specialized tools (e.g., PLIP, ProLIF).

## Quick start

### Running with Docker

All six nodes use a single unified container image built from the included Dockerfile.

```bash
docker build -t vina-workflow:latest .
```

Key tools installed: Python 3.11, BioPython, PDBFixer, OpenMM, RDKit, OpenBabel, AutoDock Vina, Meeko, P2Rank 2.4.2 (Java), Gemmi.

### Running on Silva

1. Select "AutoDock Vina Virtual Screening Workflow" from the workflow list
2. Set `pdb_id` to your target PDB identifier
3. Set `ligand_ids` to a JSON array of PubChem CIDs
4. Adjust `exhaustiveness` if needed (see Parameters section)
5. Click Run

### Test vs production settings

| Setting | Test (default) | Production |
|---------|---------------|------------|
| `pdb_id` | `"5kir"` | Your target PDB ID |
| `ligand_ids` | `["3672", "2662"]` | Your compound CIDs |
| `exhaustiveness` | `32` | `64`–`128` |
| `num_modes` | `10` | `10`–`20` |

A successful test run with the defaults docks Sunitinib and Imatinib against the 5KIR kinase structure and should complete in under 10 minutes.

## References

- Trott, O. & Olson, A.J. "AutoDock Vina: Improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading." *J. Comput. Chem.* 31(2):455–461, 2010. DOI: https://doi.org/10.1002/jcc.21334
- Eberhardt, J., Santos-Martins, D., Tillack, A.F. & Forli, S. "AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings." *J. Chem. Inf. Model.* 61(8):3891–3898, 2021. DOI: https://doi.org/10.1021/acs.jcim.1c00203
- Krivák, R. & Hoksza, D. "P2Rank: machine learning based tool for rapid and accurate prediction of ligand binding sites from protein structure." *J. Cheminform.* 10:39, 2018. DOI: https://doi.org/10.1186/s13321-018-0285-8
- [AutoDock Vina documentation](https://autodock-vina.readthedocs.io/)
- [P2Rank](https://github.com/rdk/p2rank)
- [Open Babel](https://github.com/openbabel/openbabel)
- [PDBFixer](https://github.com/openmm/pdbfixer)
