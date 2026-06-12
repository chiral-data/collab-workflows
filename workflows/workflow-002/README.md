---
doc_id: workflow-002
domain: molecular-docking
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  Virtual screening of a pre-built compound library against a protein target
  using Smina (AutoDock Vina fork). Downloads the library from Zenodo, prepares
  the receptor, docks all compounds, and ranks them by binding affinity.
tags: [smina, virtual-screening, molecular-docking, compound-library, vina]
---

# Workflow 002: Smina Virtual Screening (Library)

A four-node pipeline that screens a pre-built compound library against a protein target using Smina, a fork of AutoDock Vina with extended scoring function support. The workflow downloads a compound library from Zenodo, prepares the target protein with binding site detection from a co-crystallized ligand, docks all compounds, and produces a ranked list of hits by binding affinity.

## Overview

Smina is a fork of AutoDock Vina (Koes et al., 2013) developed at the University of Pittsburgh that adds support for custom scoring functions, faster energy minimization, and broader ligand format handling via OpenBabel. This workflow uses Smina with the original Vina scoring function (`--scoring vina`) to dock a library of pre-enumerated compounds against a target protein. The binding site is automatically defined by calculating the geometric center of a co-crystallized reference ligand, and the receptor is prepared with PDBFixer (structure cleanup) and PDB2PQR (AMBER force field charges).

The pipeline is designed for high-throughput screening of existing compound libraries where the researcher has a known protein target with a co-crystallized ligand that defines the binding site (Trott & Olson, 2010).

## When to use this workflow

Use this workflow when you want to screen a large compound library (from Zenodo) against a protein target with a known binding site defined by a co-crystallized ligand. It is suited for early-stage hit identification from pre-enumerated chemical libraries. The default library and target (PDB 5Y7J with ligand 8OL) are included for testing, and a test mode limits screening to 11 compounds for quick validation runs.

Do not use this workflow if you want to dock specific compounds by PubChem CID — use workflow-004 (AutoDock Vina) instead, which fetches individual ligands from PubChem and includes automated pocket prediction. If you want to generate and screen chemical variants of a known ligand rather than screening an external library, use workflow-003 (Smina ligand variant screening). For protein-protein docking, use workflow-016 (DiffDock-PP) or workflow-017 (LightDock).

## Architecture and data flow

```text
[01: Download Library] ──> [02: Prepare Protein] ──> [03: Virtual Screening] ──> [04: Generate Report]
        |                         |                          |                          |
  constructed_library/    {pdb_id}_AB_chains_fixed.pdb   docking_results/       results/docking_ranking.txt
                          config.txt                     *_docked.sdf
                          *.pqr                          *_log.txt
```

Nodes run sequentially: 01 → 02 → 03 → 04.

## Input requirements

- **No local files required** — the compound library is downloaded from Zenodo and the protein structure from RCSB PDB automatically.
- **Target protein:** Specified by PDB ID (default: `5Y7J`). Must contain a co-crystallized reference ligand for binding site definition.
- **Reference ligand:** Specified by three-letter residue name (default: `8OL`). Used to calculate the binding site center coordinates; must be present in the PDB structure.
- **Compound library:** Downloaded as `constructed_library.zip` from Zenodo (record 17374422), containing SDF files named `clean_drug*.sdf`.

## Workflow nodes

### Node 01: Download Library

**Goal:** Download the pre-built compound library from Zenodo.

**Process:** Fetches `constructed_library.zip` from `https://zenodo.org/records/17374422` via HTTP streaming with progress reporting, then extracts the ZIP archive to the working directory. The archive contains a `constructed_library/` folder with individual SDF ligand files.

**Scientific notes:** The library consists of pre-enumerated drug-like compounds in SDF format with 3D coordinates. These are ready for docking without further conformer generation.

**Outputs:**
- `constructed_library/` — directory containing `clean_drug*.sdf` files

### Node 02: Prepare Protein

**Goal:** Download, clean, and prepare the target protein structure for docking.

**Process:** Downloads the PDB file from RCSB, extracts chains A and B (if present) along with the reference ligand, calculates the binding site center as the geometric centroid of the reference ligand's atom coordinates using NumPy, and writes a `config.txt` with the docking grid parameters (center coordinates, 15 × 15 × 15 Å grid, exhaustiveness 8). Then applies PDBFixer to replace non-standard residues, remove all heterogens including waters (`removeHeterogens(keepWater=False)`), and add missing hydrogens at default pH 7.0. Finally, assigns AMBER force field charges and radii using PDB2PQR (`--ff=AMBER`).

**Scientific notes:** The binding site is defined entirely from the co-crystallized ligand position — this requires that the reference ligand is present in the deposited PDB structure. The grid size is fixed at 15 × 15 × 15 Å, which is appropriate for typical small-molecule binding pockets but may be too small for large or shallow binding sites. Unlike workflow-004 (which uses P2Rank for blind pocket prediction), this workflow requires prior knowledge of where the ligand binds. PDB2PQR assigns AMBER94 partial charges, which are used by Smina's scoring function.

**Outputs:**
- `{pdb_id}_AB_chains_fixed.pdb` — cleaned receptor structure
- `{pdb_id}_AB_chains_fixed.pqr` — receptor with AMBER charges
- `config.txt` — docking grid configuration (center, size, exhaustiveness)
- `results/{pdb_id}_AB_chains_fixed.pdb` — copy for downstream use

### Node 03: Virtual Screening

**Goal:** Dock all library compounds against the prepared receptor using Smina.

**Process:** Iterates over SDF files in the library directory, invoking Smina for each compound with the receptor, config file, and `--scoring vina` flag. In test mode (default), only compounds matching `clean_drug108*.sdf` are docked (~11 compounds); in full mode, all `clean_drug*.sdf` files are processed. Each docking has a 5-minute timeout. Binding affinities are extracted from log files (mode 1 / best pose).

**Scientific notes:** Smina uses the Vina scoring function here, which combines steric (Gaussian + repulsion), hydrophobic, and hydrogen-bonding terms optimized via iterated local search with BFGS. The exhaustiveness parameter (default 8, read from config.txt but overridable via environment variable) controls the number of independent search runs. Unlike workflow-004 which uses exhaustiveness 32, this workflow defaults to 8 for faster throughput across many compounds — suitable for initial screening but less thorough per compound.

**Outputs:**
- `docking_results/{compound}_docked.sdf` — docked poses for each compound
- `docking_results/{compound}_log.txt` — Smina scoring output with binding affinities

### Node 04: Generate Report

**Goal:** Parse docking results and produce a ranked compound list.

**Process:** Reads all `*_log.txt` files from the docking results, extracts the mode 1 (best pose) binding affinity from each, sorts compounds by binding energy in ascending order (strongest binding first), writes the ranking to a text file, and copies the top-ranked compound's docked SDF to the results directory.

**Scientific notes:** Only the best pose (mode 1) per compound is used for ranking. This is standard for initial screening but means that alternative binding modes are not considered. For a more detailed analysis of promising hits, re-dock with higher exhaustiveness and num_modes in workflow-004.

**Outputs:**
- `results/docking_ranking.txt` — ranked list of all compounds with binding energies
- `results/{top_compound}_docked.sdf` — docked structure of the best-scoring compound

## Parameters

### pdb_id

- **Type:** string
- **Default:** `"5Y7J"`
- **Node:** 02 (workflow-level)
- **Description:** RCSB PDB identifier for the target protein.
- **Guidance:** The structure must contain the reference ligand specified by `ligand_name`. Use a structure with good resolution (< 2.5 Å).

### ligand_name

- **Type:** string
- **Default:** `"8OL"`
- **Node:** 02 (workflow-level)
- **Description:** Three-letter residue name of the co-crystallized reference ligand used to define the binding site center.
- **Guidance:** Must match a HETATM residue in the PDB file. The geometric center of this ligand's atoms becomes the docking grid center.

### test_mode

- **Type:** boolean
- **Default:** `true`
- **Node:** 03
- **Description:** When true, screens only ~11 compounds (matching `clean_drug108*.sdf`); when false, screens the full library.

**Trade-off:** Test mode completes in minutes; full library screening can take hours depending on library size and exhaustiveness.

**Test vs production:** Set to `false` for actual screening campaigns.

### exhaustiveness

- **Type:** integer
- **Default:** `8`
- **Node:** 03
- **Description:** Number of independent Smina search runs per compound.
- **Guidance:** The default of 8 prioritizes throughput for library screening. For more thorough docking of selected hits, increase to 32 or re-dock promising compounds in workflow-004.

### top_n

- **Type:** integer
- **Default:** `10`
- **Node:** 04
- **Description:** Number of top-ranked compounds to display in the report output.

## Outputs and interpretation

### Binding affinity (kcal/mol)

Compounds are ranked by predicted binding free energy from Smina's Vina scoring function. More negative values indicate stronger predicted binding.

| Range | Interpretation |
|-------|---------------|
| < −10 kcal/mol | Very strong binding (rare, verify carefully) |
| −8 to −10 kcal/mol | Strong binding — promising hit |
| −6 to −8 kcal/mol | Moderate binding — worth investigating |
| > −6 kcal/mol | Weak binding — likely not viable |

The Vina scoring function has a standard error of ~2.85 kcal/mol against experimental data. Differences < 1 kcal/mol between compounds are not meaningful. For high-confidence results, re-dock top hits with higher exhaustiveness.

### docking_ranking.txt

A text file listing all successfully docked compounds sorted by binding energy (strongest first), with rank number, compound name, and affinity in kcal/mol.

## Quick start

### Running with Docker

All nodes use the same container image:

```bash
docker pull ghcr.io/chiral-data/smina:2025_10_17_v2
```

### Running on Silva

1. Select "In-Silico Virtual Screening" from the workflow list
2. Set `pdb_id` and `ligand_name` for your target (or keep defaults for testing)
3. Set `test_mode` to `false` for full library screening
4. Adjust `exhaustiveness` if needed
5. Click Run

### Test vs production settings

| Setting | Test (default) | Production |
|---------|---------------|------------|
| `pdb_id` | `"5Y7J"` | Your target PDB ID |
| `ligand_name` | `"8OL"` | Your reference ligand residue name |
| `test_mode` | `true` (~11 compounds) | `false` (full library) |
| `exhaustiveness` | `8` | `16`–`32` |

A successful test run docks ~11 compounds against PDB 5Y7J and completes in a few minutes.

## References

- Koes, D.R., Baumgartner, M.P. & Camacho, C.J. "Lessons Learned in Empirical Scoring with smina from the CSAR 2011 Benchmarking Exercise." *J. Chem. Inf. Model.* 53(8):1893–1904, 2013. DOI: https://doi.org/10.1021/ci300604z
- Trott, O. & Olson, A.J. "AutoDock Vina: Improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading." *J. Comput. Chem.* 31(2):455–461, 2010. DOI: https://doi.org/10.1002/jcc.21334
- Dolinsky, T.J. et al. "PDB2PQR: an automated pipeline for the setup of Poisson-Boltzmann electrostatics calculations." *Nucleic Acids Res.* 32(suppl_2):W665–W667, 2004. DOI: https://doi.org/10.1093/nar/gkh381
- [Smina](https://sourceforge.net/projects/smina/)
- [PDBFixer](https://github.com/openmm/pdbfixer)
- [PDB2PQR](https://pdb2pqr.readthedocs.io/)
