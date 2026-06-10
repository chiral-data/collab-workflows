---
doc_id: workflow-003
domain: molecular-docking
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  Virtual screening of systematically generated ligand variants against a
  protein target using Smina. Extracts a co-crystallized ligand, applies
  functional group substitutions via RDKit, and docks all variants.
tags: [smina, virtual-screening, molecular-docking, ligand-modification, rdkit, lead-optimization]
---

# Workflow 003: Smina Virtual Screening (Ligand Variants)

A six-node pipeline that extracts a co-crystallized ligand from a PDB structure, generates chemical variants through systematic functional group substitutions using RDKit, and docks all variants against the prepared receptor using Smina. Designed for lead optimization scenarios where a researcher wants to explore how modifications to a known binder affect binding affinity.

## Overview

This workflow automates a common medicinal chemistry exploration: starting from a known protein-ligand complex, systematically modifying the ligand's functional groups and predicting how each modification affects binding. RDKit performs substructure-based replacements (e.g., hydroxyl → amine, halogen → hydroxyl) to generate a focused variant library, then generates 3D conformers using the ETKDGv3 algorithm with MMFF force field optimization. Smina (a fork of AutoDock Vina) docks each variant using the Vina scoring function, and results are ranked by predicted binding energy.

The default system is PDB 4OHU (chain A with NAD cofactor) and ligand 2TK, but any PDB with a co-crystallized ligand can be used (Koes et al., 2013; Trott & Olson, 2010).

## When to use this workflow

Use this workflow when you have a known protein-ligand complex and want to explore how systematic chemical modifications to the ligand affect binding affinity. It is designed for lead optimization — you already know a ligand that binds, and you want to test functional group replacements. The workflow handles ligand extraction, variant generation, receptor preparation (including cofactor preservation), and docking automatically.

Do not use this workflow for screening an external compound library — use workflow-002 (Smina library screening from Zenodo) instead. If you want to dock specific compounds by PubChem CID without generating variants, use workflow-004 (AutoDock Vina with P2Rank pocket prediction). For protein-protein docking, use workflow-016 (DiffDock-PP) or workflow-017 (LightDock).

## Architecture and data flow

```text
[01: Download PDB]
       │
       ▼
[02: Protein Preparation] ──────────────────────────┐
       │                                             │
       ▼                                             │
[03: Ligand Extraction]                              │
       │                                             │
       ▼                                             │
[04: Ligand Modification]                            │
       │                                             │
       ▼                                             ▼
[05: Virtual Screening] ◄── receptor + config + ligand_library
       │
       ▼
[06: Report]
```

Dependencies: 02 ← 01, 03 ← 02, 04 ← 03, 05 ← 02 + 04, 06 ← 05. Nodes 02 and 04 feed into Node 05 from separate branches.

## Input requirements

- **No local files required** — the PDB structure is downloaded from RCSB automatically.
- **Target protein:** Specified by PDB ID (default: `4OHU`). Must contain a co-crystallized ligand for extraction.
- **Reference ligand:** Specified by three-letter residue name (default: `2TK`). This ligand is extracted, used to define the binding site, and serves as the base molecule for variant generation.
- **Constraints:** The ligand must be present in the PDB file and parseable by RDKit. The SMILES derived from the extracted ligand must contain at least one functional group matching the substitution patterns (hydroxyl, amine, carboxyl, halogen, or nitro).

## Workflow nodes

### Node 01: Download PDB

**Goal:** Download the target protein structure from RCSB PDB.

**Process:** Fetches the PDB file for the given `pdb_id` parameter from `https://files.rcsb.org/download/{pdb_id}.pdb` and saves it to the outputs directory.

**Scientific notes:** The full deposited structure is downloaded, including all chains, cofactors, ligands, and waters. Chain and ligand selection happens in subsequent nodes.

**Outputs:**
- `{pdb_id}.pdb` — raw crystal structure

### Node 02: Protein Preparation

**Goal:** Extract the relevant chain and cofactor, calculate the binding site, fix the structure, and assign charges.

**Process:** This node runs `screening_preparation.py`, which performs several steps:
1. Extracts chain A and the NAD cofactor from the PDB (excluding the target ligand 2TK) using BioPython
2. Calculates the binding site center as the geometric centroid of the target ligand's atom coordinates from the original PDB
3. Computes a dynamic grid size from the ligand's spatial extent plus 8 Å padding (minimum 20 Å per side) and writes `config.txt` with center, size, exhaustiveness, num_modes, and energy_range
4. Fixes the structure with PDBFixer (replace non-standard residues, remove heterogens except water, add hydrogens)
5. Reattaches NAD cofactor residues if PDBFixer removed them (parses original PDB for HETATM/ATOM lines containing "NAD")
6. Optionally assigns AMBER charges via PDB2PQR (`--ff=AMBER`); gracefully skips if PDB2PQR fails on non-standard residues like NAD

**Scientific notes:** The NAD cofactor is essential for this target (4OHU is an NAD-dependent enzyme). PDBFixer's `removeHeterogens(keepWater=False)` removes all heterogens including cofactors, so the explicit reattach step preserves the biologically relevant NAD. The grid size is dynamically computed from the ligand geometry, unlike workflow-002 which uses a fixed 15 × 15 × 15 Å grid.

**Outputs:**
- `{pdb_id}_A_NAD_fixed_with_NAD.pdb` — cleaned receptor with NAD cofactor
- `config.txt` — docking grid configuration (dynamic center and size)
- `results/{pdb_id}_A_NAD_fixed_with_NAD.pdb` — receptor copy for downstream use

### Node 03: Ligand Extraction

**Goal:** Extract the target ligand from the PDB and save as SDF with original crystallographic coordinates.

**Process:** Uses BioPython to isolate residues matching the `ligand_name` parameter from the chain A + NAD structure, saves to a temporary PDB, then reads with RDKit (`MolFromPDBFile` with `removeHs=False`) to convert to SDF format preserving the original 3D coordinates.

**Scientific notes:** Preserving the crystallographic coordinates is important — this conformation represents the experimentally observed binding pose. The extracted SDF serves as the base molecule for variant generation in Node 04. If RDKit cannot parse the extracted ligand PDB (e.g., unusual bonding or missing atoms), the node will fail.

**Outputs:**
- `{ligand_name}.sdf` — extracted ligand with crystallographic coordinates

### Node 04: Ligand Modification

**Goal:** Generate a focused library of chemical variants by systematic functional group substitution.

**Process:** Reads the extracted ligand SDF, converts to canonical SMILES, then applies 9 predefined substructure-based replacements using RDKit's `ReplaceSubstructs`:

| Substitution | SMARTS pattern | Replacement |
|-------------|----------------|-------------|
| hydroxyl → amine | `[OX2H]` | `N` |
| hydroxyl → thiol | `[OX2H]` | `S` |
| hydroxyl → halogen | `[OX2H]` | `Cl` |
| amine → hydroxyl | `[NX3;H2]` | `O` |
| carboxyl → amide | `C(=O)[OH]` | `C(=O)N` |
| carboxyl → ester | `C(=O)[OH]` | `C(=O)OC` |
| halogen → hydroxyl | `[F,Cl,Br,I]` | `O` |
| halogen → amine | `[F,Cl,Br,I]` | `N` |
| nitro → amine | `[NX3](=O)=O` | `N` |

Each unique variant is embedded in 3D using `AllChem.EmbedMolecule` with ETKDGv3, then optimized with the MMFF (Merck Molecular Force Field). A 2D SVG overview (`variants.svg`) is also generated.

**Scientific notes:** The substitutions follow standard bioisosteric replacement strategies used in medicinal chemistry. Only functional groups present in the base molecule will produce variants — if the ligand has no hydroxyl groups, the hydroxyl→X substitutions yield nothing. The number of variants depends entirely on the base molecule's functional group inventory. ETKDGv3 is RDKit's distance geometry method with experimental torsion angle preferences for generating realistic 3D conformers.

**Outputs:**
- `ligand_library/original.sdf` — base molecule with generated 3D coordinates
- `ligand_library/{substitution_name}.sdf` — one SDF per successful variant
- `variants.svg` — 2D depiction grid of all variants

### Node 05: Virtual Screening

**Goal:** Dock all ligand variants against the prepared receptor using Smina.

**Process:** Runs in two stages. First, `screening_preparation.py` generates the docking config and prepares the receptor (see Node 02 — this script is re-run in Node 05's context to produce `config.txt` and the receptor PDBQT). Then `in_silico_screening.py` converts the receptor to PDBQT format using MGLTools' `prepare_receptor4.py` (with `-A hydrogens` flag), and docks each SDF file in `ligand_library/` against the receptor using Smina with `--scoring vina` and `--num_modes 1`. Each docking has a 10-minute timeout.

**Scientific notes:** Only the single best pose (`--num_modes 1`) is generated per variant, which is appropriate for comparing relative affinities across a focused library. The Vina scoring function is used explicitly. The receptor is converted to PDBQT using MGLTools' `prepare_receptor4.py` rather than OpenBabel (as in workflow-004), which better handles non-standard residues like NAD. Grid parameters (center, size, exhaustiveness, num_modes, energy_range) come from config.txt.

**Outputs:**
- `docking_results/{variant}_docked.sdf` — docked pose for each variant
- `docking_results/{variant}_docking.log` — Smina scoring output

### Node 06: Report

**Goal:** Parse docking results and produce a ranked variant list.

**Process:** Reads all `*_docking.log` files, extracts the mode 1 binding affinity from each, sorts variants by binding energy (ascending = strongest binding first), writes the ranking to a text file, and copies the top-ranked variant's docked SDF to the results directory.

**Scientific notes:** Comparing variants against the original ligand's docking score indicates whether the modification is predicted to improve or worsen binding. However, Vina's scoring function has a standard error of ~2.85 kcal/mol, so only differences > 1 kcal/mol should be considered meaningful. For promising modifications, validate with more rigorous methods (e.g., FEP or experimental assays).

**Outputs:**
- `results/docking_ranking.txt` — ranked list of all variants with binding energies
- `results/{top_variant}_docked.sdf` — docked structure of the best-scoring variant

## Parameters

### pdb_id

- **Type:** string
- **Default:** `"4OHU"`
- **Node:** workflow-level (used in Nodes 01, 02, 05)
- **Description:** RCSB PDB identifier for the target protein.
- **Guidance:** The structure must contain the reference ligand specified by `ligand_name`. The code extracts chain A and the NAD cofactor specifically, so modify the source code if your target uses a different chain or cofactor.

### ligand_name

- **Type:** string
- **Default:** `"2TK"`
- **Node:** workflow-level (used in Nodes 03, 04)
- **Description:** Three-letter residue name of the co-crystallized ligand to extract and modify.
- **Guidance:** Must match a HETATM residue in the PDB file. The ligand must be parseable by RDKit and contain functional groups matching at least some of the substitution patterns.

### exhaustiveness

- **Type:** integer
- **Default:** `8`
- **Node:** 05
- **Description:** Number of independent Smina search runs per variant.
- **Guidance:** Default of 8 is suitable for initial screening of a small variant library. For more thorough sampling, increase to 16–32.

### num_modes

- **Type:** integer
- **Default:** `5`
- **Node:** 05
- **Description:** Maximum number of binding poses to generate per variant.
- **Guidance:** The actual docking command overrides this with `--num_modes 1` (only best pose). To generate multiple poses, modify `in_silico_screening.py`.

### energy_range

- **Type:** integer
- **Default:** `4`
- **Node:** 05
- **Description:** Energy window in kcal/mol — only poses within this range of the best pose are retained.

## Outputs and interpretation

### Binding affinity (kcal/mol)

Variants are ranked by predicted binding free energy. More negative = stronger binding.

| Range | Interpretation |
|-------|---------------|
| < −10 kcal/mol | Very strong binding (rare, verify carefully) |
| −8 to −10 kcal/mol | Strong binding — promising modification |
| −6 to −8 kcal/mol | Moderate binding |
| > −6 kcal/mol | Weak binding |

Compare each variant's score to the original ligand's score. If a substitution improves the score by > 1 kcal/mol, it is a candidate for further investigation.

### variants.svg

A 2D grid depiction of the original molecule and all generated variants, labeled by substitution type. Useful for quickly inspecting which modifications were applied.

### docking_ranking.txt

A text file listing all successfully docked variants sorted by binding energy (strongest first), with rank, compound name, and affinity in kcal/mol.

## Quick start

### Running with Docker

All nodes use the same container image:

```bash
docker pull ghcr.io/chiral-data/smina:2025_11_06
```

### Running on Silva

1. Select "In-silico Screening by Diego" from the workflow list
2. Set `pdb_id` and `ligand_name` for your target (or keep defaults: 4OHU / 2TK)
3. Adjust `exhaustiveness` if needed
4. Click Run

### Test vs production settings

| Setting | Test (default) | Production |
|---------|---------------|------------|
| `pdb_id` | `"4OHU"` | Your target PDB ID |
| `ligand_name` | `"2TK"` | Your ligand residue name |
| `exhaustiveness` | `8` | `16`–`32` |
| `energy_range` | `4` | `4`–`5` |

A successful test run generates ~7–10 variants of 2TK, docks them against 4OHU with NAD cofactor, and completes in under 30 minutes.

## References

- Koes, D.R., Baumgartner, M.P. & Camacho, C.J. "Lessons Learned in Empirical Scoring with smina from the CSAR 2011 Benchmarking Exercise." *J. Chem. Inf. Model.* 53(8):1893–1904, 2013. DOI: https://doi.org/10.1021/ci300604z
- Trott, O. & Olson, A.J. "AutoDock Vina: Improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading." *J. Comput. Chem.* 31(2):455–461, 2010. DOI: https://doi.org/10.1002/jcc.21334
- Riniker, S. & Landrum, G.A. "Better Informed Distance Geometry: Using What We Know To Improve Conformation Generation." *J. Chem. Inf. Model.* 55(12):2562–2574, 2015. DOI: https://doi.org/10.1021/acs.jcim.5b00654
- [Smina](https://sourceforge.net/projects/smina/)
- [RDKit](https://www.rdkit.org/)
- [PDBFixer](https://github.com/openmm/pdbfixer)
