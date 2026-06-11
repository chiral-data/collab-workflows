---
doc_id: workflow-009
domain: molecular-docking
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  Active learning-guided lead optimization using FEGrow for functional group
  elaboration and gnina for CNN-based molecular docking. Builds a combinatorial
  chemical space from a scaffold ligand, then iteratively selects and evaluates
  candidates using Bayesian surrogate models.
tags: [fegrow, active-learning, gnina, docking, lead-optimization, chemical-space, drug-discovery]
---

# Workflow 009: FEGrow Active Learning

Active learning-guided lead optimization using [FEGrow](https://github.com/cole-group/FEgrow) and [gnina](https://github.com/gnina/gnina). This workflow takes a ligand scaffold and a protein target, builds a combinatorial chemical space by attaching linkers and R-groups, then uses active learning with gnina docking to efficiently explore that space and identify top-binding compounds.

## Overview

FEGrow is a molecular builder for functional group elaboration of lead compounds (Bieniek et al., 2022). Given a scaffold molecule with a defined attachment point, it enumerates candidate molecules by combining linker fragments and R-group substituents from built-in chemical libraries. Rather than docking every candidate, this workflow uses active learning to prioritize evaluation: a surrogate model (Gaussian Process, Random Forest, or Linear) is trained on initial docking results and iteratively selects the most promising candidates for evaluation using gnina, a CNN-based molecular docking tool that extends AutoDock Vina with convolutional neural network scoring (McNutt et al., 2021).

## When to use this workflow

Use this workflow when you have a known ligand (hit or lead compound) bound to a protein target and want to explore chemical modifications to improve binding affinity. You need an SDF file of the ligand and a PDB file of the protein. The workflow is designed for lead optimization — systematic elaboration of a starting molecule rather than de novo design or large-library virtual screening.

For de novo protein design, use workflow-013 (BoltzGen). For virtual screening of existing compound libraries, use workflow-004 (AutoDock Vina) or workflow-002/003 (Smina). For binding pocket identification before docking, use workflow-007 (Pocketeer).

## Architecture and data flow

```text
[01: Ligand Upload] ──> [02: Scaffold Creation] ──┐
                                                   ├──> [05: Chemical Space] ──> [06: Active Learning] ──> [07: Report]
[03: Protein Upload] ──> [04: Protein Prep] ──────┘
```

Nodes 01 and 03 run in parallel. Node 02 depends on 01; Node 04 depends on 03. Node 05 requires both 02 and 04. Nodes 06 and 07 run sequentially after 05.

## Input requirements

- **Ligand SDF file:** A 3D structure of the starting ligand in SDF format. Must contain explicit hydrogens and valid atom coordinates. The default test ligand is provided at `.chiral/test_inputs/ligand.sdf`.
- **Protein PDB file:** The target protein structure in PDB format. The default test protein is provided at `.chiral/test_inputs/protein.pdb`.
- **Attachment point:** The atom index on the scaffold where new functional groups will be attached (set via the `attachment_id` parameter).

## Workflow nodes

### Node 01: Ligand Upload

**Goal:** Validate the input ligand and extract its SMILES representation.

**Process:** Loads the SDF file with RDKit, validates that the molecule can be parsed, assigns atom map numbers for tracking, and writes both a canonical SMILES file and a cleaned SDF. Generates an HTML visualization of the ligand with atom indices labeled to help identify the attachment point.

**Outputs:**
- `*.smi` -- ligand SMILES
- `*.sdf` -- validated ligand SDF
- `ligand_viz.html` -- interactive ligand visualization with atom indices

### Node 02: Scaffold Creation

**Goal:** Create a FEGrow scaffold molecule with a defined attachment point.

**Process:** Loads the validated ligand from Node 01 and creates an FEGrow scaffold by setting the atom at `attachment_id` to atomic number 0 (wildcard), marking it as the site for fragment attachment. The scaffold is serialized as a pickle file for downstream use. Generates an HTML visualization showing the scaffold with the attachment point highlighted.

**Scientific notes:** The attachment point defines where new chemical groups are grown onto the scaffold. Choosing the right attachment point is critical -- it should be at a position where modifications can access the binding pocket without disrupting key protein-ligand interactions. The atom index visualization from Node 01 helps identify suitable positions.

**Outputs:**
- `scaffold.pkl` -- serialized FEGrow scaffold object
- `scaffold_viz.html` -- scaffold visualization with attachment point
- `*.smi` -- scaffold SMILES (passed through)

### Node 03: Protein Upload

**Goal:** Validate the input protein structure.

**Process:** Loads the PDB file using ProDy, validates the structure, and copies it to outputs.

**Outputs:**
- `*.pdb` -- validated protein PDB file

### Node 04: Protein Preparation

**Goal:** Clean the protein structure for docking.

**Process:** Removes water molecules, nucleic acids, and hetero atoms (ligands, cofactors) from the protein structure. Runs `fegrow.fix_receptor()` to add missing hydrogens and optimize side-chain conformations for docking.

**Scientific notes:** Removing co-crystallized ligands and water molecules ensures the docking algorithm evaluates binding in the apo pocket. The receptor fix step repairs common PDB issues (missing atoms, incorrect protonation) that would otherwise cause docking failures.

**Outputs:**
- `rec_final.pdb` -- cleaned and prepared protein structure

### Node 05: Chemical Space Creation

**Goal:** Build the combinatorial chemical space from the scaffold.

**Process:** Loads the scaffold from Node 02 and the prepared receptor from Node 04. Iterates over `num_linkers` linker fragments and `num_rgroups` R-group substituents from FEGrow's built-in chemical libraries, combining each linker-R-group pair with the scaffold to enumerate candidate molecules. Uses a Dask LocalCluster for parallel construction. The resulting chemical space is serialized as a `ChemSpace` object.

**Scientific notes:** The chemical space size is at most `num_linkers x num_rgroups` molecules (some combinations may fail valence or geometry checks). With the defaults (10 x 10), the space contains up to 100 candidates. Increasing these values explores more chemistry but proportionally increases the evaluation cost in Node 06.

**Outputs:**
- `chemspace.pkl` -- serialized ChemSpace object with all valid candidates
- `rec_final.pdb` -- receptor (passed through for Node 06)
- `chemspace_viz.html` -- visualization of the chemical space

### Node 06: Active Learning Docking

**Goal:** Efficiently identify top-binding compounds through iterative surrogate-guided docking.

**Process:** Runs active learning over the chemical space:
1. **Initial sampling:** Randomly selects `initial_molecules` candidates and evaluates them with gnina docking (5 conformations per molecule)
2. **Surrogate training:** Fits the selected surrogate model (Gaussian Process by default) on the docking results
3. **Iterative cycles:** For `num_cycles` iterations, uses the query strategy (UCB by default) to select the next `molecules_per_cycle` most promising candidates, evaluates them with gnina (10 conformations per molecule in cycles), and updates the surrogate model
4. **Result compilation:** Writes an SDF file with all evaluated molecules and their docking scores, plus per-iteration CSV logs

**Scientific notes:** Active learning reduces the number of expensive docking evaluations needed to find top binders. The surrogate model predicts binding affinity from molecular features, and the acquisition function balances exploitation (selecting predicted best) with exploration (selecting uncertain). UCB (Upper Confidence Bound) adds a multiple of the prediction uncertainty to the predicted value. Greedy selects purely by predicted score. PI (Probability of Improvement) and EI (Expected Improvement) use the probability distribution to balance risk. Gaussian Process surrogates provide calibrated uncertainty estimates; Random Forest and Linear models are faster but less principled in their uncertainty quantification.

**Outputs:**
- `chemspace_evaluated.sdf` -- all evaluated molecules with docking scores
- `iteration_*.csv` -- per-cycle evaluation logs

### Node 07: Report Generation

**Goal:** Rank evaluated compounds and generate a final report.

**Process:** Loads the evaluated chemical space SDF, converts pK scores to binding energy in kcal/mol using the thermodynamic relationship (deltaG = -1.36 x pK at 298 K), ranks compounds by binding energy (most negative = strongest binder), and extracts the top `top_n` compounds. Generates an SDF file of top compounds, a CSV report with SMILES, scores, and classifications, and a text summary.

**Scientific notes:** The pK-to-deltaG conversion derives from deltaG = -RT ln(K) = -2.303 RT log10(K). At 298 K, RT = 0.593 kcal/mol, giving deltaG = -1.364 x pK kcal/mol. Compounds are classified by binding energy: excellent (< -8 kcal/mol, nanomolar binders), good (-8 to -5, low micromolar), fair (-5 to 0, weak binding), poor (>= 0, no meaningful binding).

**Outputs:**
- `top_compounds.sdf` -- SDF file of the top-ranked compounds
- `top_compounds_report.csv` -- CSV with SMILES, binding energy, and classification
- `summary.txt` -- text summary of the active learning run
- `report_viz.html` -- interactive report visualization

## Parameters

### attachment_id

- **Type:** integer
- **Default:** `27`
- **Description:** Atom index on the scaffold ligand where new functional groups are attached. Use the atom index visualization from Node 01 to identify the correct position.
- **Guidance:** This is molecule-specific. The default works for the test ligand. For your own ligands, run Node 01 first, inspect `ligand_viz.html` to find the atom index at the desired growth position.

### num_linkers

- **Type:** integer
- **Default:** `10`
- **Description:** Number of linker fragments to sample from FEGrow's built-in library. Linkers connect the scaffold to R-groups.
- **Guidance:** Higher values explore more structural diversity but increase chemical space size proportionally. For thorough exploration, set to 20-50.

### num_rgroups

- **Type:** integer
- **Default:** `10`
- **Description:** Number of R-group substituents to sample from FEGrow's built-in library. Combined with each linker to form candidate molecules.
- **Guidance:** Higher values explore more chemical diversity. Chemical space size = num_linkers x num_rgroups.

### initial_molecules

- **Type:** integer
- **Default:** `10`
- **Description:** Number of molecules randomly selected for the initial docking round before active learning begins.
- **Guidance:** Should be large enough to give the surrogate model a reasonable training set. Too few (< 5) leads to poor initial model quality; too many wastes the active learning budget on random exploration.

### num_cycles

- **Type:** integer
- **Default:** `3`
- **Description:** Number of active learning iterations after the initial random round. Each cycle trains the surrogate model and selects the next batch of candidates.
- **Guidance:** More cycles allow the surrogate to refine its predictions. Diminishing returns after the model has explored most of the promising chemical space.

### molecules_per_cycle

- **Type:** integer
- **Default:** `50`
- **Description:** Number of molecules selected and evaluated per active learning cycle.
- **Guidance:** Larger batches cover more space per cycle but reduce the number of model updates. For large chemical spaces, increase this value.

### model_type

- **Type:** enum
- **Default:** `gaussian_process`
- **Node:** 06
- **Description:** Surrogate model used to predict binding affinity between docking evaluations.

| Value | Description |
|-------|-------------|
| `gaussian_process` (default) | Gaussian Process regression with calibrated uncertainty estimates. Best for small-to-medium chemical spaces |
| `random_forest` | Ensemble of decision trees. Faster than GP for large datasets but less principled uncertainty |
| `linear` | Linear regression. Fastest but assumes linear structure-activity relationship |

### query_type

- **Type:** enum
- **Default:** `UCB`
- **Node:** 06
- **Description:** Acquisition function that decides which molecules to evaluate next.

| Value | Description |
|-------|-------------|
| `UCB` (default) | Upper Confidence Bound. Balances predicted score with uncertainty. Good general choice |
| `Greedy` | Selects purely by predicted score. Exploits current knowledge, may miss better regions |
| `PI` | Probability of Improvement. Selects molecules most likely to beat the current best |
| `EI` | Expected Improvement. Balances probability and magnitude of improvement |

### top_n

- **Type:** integer
- **Default:** `10`
- **Node:** 07
- **Description:** Number of top-ranked compounds to include in the final report.

## Outputs and interpretation

### Binding energy classification

Compounds in the final report are classified by predicted binding energy (deltaG in kcal/mol):

| Classification | Binding energy | Interpretation |
|---------------|----------------|----------------|
| Excellent | < -8 kcal/mol | Nanomolar-range binding. Strong drug candidates |
| Good | -8 to -5 kcal/mol | Low micromolar binding. Viable leads for further optimization |
| Fair | -5 to 0 kcal/mol | Weak binding. May need significant modification |
| Poor | >= 0 kcal/mol | No meaningful binding predicted |

### top_compounds_report.csv

The primary results file. Contains SMILES, predicted binding energy (kcal/mol), gnina docking scores, and binding energy classification for the top-ranked compounds.

### top_compounds.sdf

3D structures of the top-ranked compounds in SDF format. Can be loaded into molecular viewers (PyMOL, ChimeraX) or used as input for downstream workflows.

### iteration_*.csv

Per-cycle logs showing which molecules were selected and their docking scores. Useful for analyzing the active learning convergence -- scores should generally improve across iterations.

## Quick start

### Running with Docker

```bash
docker pull ghcr.io/chiral-data/fegrow:2026_01_10_v4
```

### Running on Silva

1. Select "FEGrow Active Learning Workflow" from the workflow list
2. Upload your ligand SDF and protein PDB files
3. Run Node 01 first to identify atom indices, then set `attachment_id`
4. Adjust chemical space size (`num_linkers`, `num_rgroups`) and active learning parameters as needed
5. Click Run

### Test vs production settings

| Setting | Test (default) | Production |
|---------|---------------|------------|
| `num_linkers` | `10` | `20`--`50` |
| `num_rgroups` | `10` | `20`--`50` |
| `initial_molecules` | `10` | `20`--`30` |
| `num_cycles` | `3` | `5`--`10` |
| `molecules_per_cycle` | `50` | `50`--`100` |

The default test run uses the included ligand and protein files. A successful run produces a report with the top 10 compounds ranked by binding energy, along with per-iteration CSV logs showing active learning convergence.

## Troubleshooting

### Scaffold creation fails

If Node 02 fails, the `attachment_id` likely points to an invalid atom (e.g., a hydrogen or an atom at a ring junction where substitution violates valence rules). Inspect `ligand_viz.html` from Node 01 and choose a different attachment point.

### All compounds score poorly

If all top compounds have binding energy >= 0, the scaffold may not be well-positioned in the binding pocket, or the attachment point faces away from the pocket. Verify that the input ligand represents a reasonable binding pose and that the attachment point is oriented toward the binding site.

## References

- Bieniek, M. K. et al. "An open-source molecular builder and free energy preparation workflow." *Communications Chemistry* 5:136, 2022. DOI: https://doi.org/10.1038/s42004-022-00754-9
- McNutt, A. T. et al. "GNINA 1.0: molecular docking with deep learning." *Journal of Cheminformatics* 13:43, 2021. DOI: https://doi.org/10.1186/s13321-021-00522-2
- Graff, D. E. et al. "Accelerating high-throughput virtual screening through molecular pool-based active learning." *Chemical Science* 12:7866--7881, 2021. DOI: https://doi.org/10.1039/D0SC06805E
- [FEGrow GitHub](https://github.com/cole-group/FEgrow)
- [gnina GitHub](https://github.com/gnina/gnina)
