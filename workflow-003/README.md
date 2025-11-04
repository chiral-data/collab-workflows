# Workflow-002: In-Silico Virtual Screening

## Overview

This workflow automates **virtual screening** of ligand variants against a protein target using **Smina**, an enhanced fork of AutoDock Vina.  
It integrates molecular structure extraction, ligand modification, receptor preparation, docking, and result ranking into a reproducible computational pipeline.

The target system modeled here is based on **PDB ID: 4OHU**, which includes **chain A** and the **NAD cofactor**, and the ligand of interest **2TK**.

---

## Pipeline Structure

### Step 1: Ligand Extraction (`extract_2TK_from_pdb.py`)
Extracts the 2TK ligand coordinates from the downloaded PDB file.

**Process:**
- Reads the `4OHU.pdb` structure.
- Isolates the `2TK` residue from all chains.
- Exports the ligand as `2TK.sdf` (preserving its 3D coordinates).

**Output:**
```
2TK.sdf
```

**Usage Example:**
```bash
python3 extract_2TK_from_pdb.py
```

---

### Step 2: Ligand Modification (`generate_ligand_mods.py`)
Generates a **set of chemical variants** based on the extracted 2TK molecule.

**Process:**
1. Defines a base SMILES representation of 2TK.  
2. Applies systematic functional group substitutions:
   - Hydroxyl ↔ Amine/Thiol/Halogen
   - Carboxyl ↔ Amide/Ester
   - Halogen ↔ Hydroxyl/Amine
   - Nitro ↔ Amine  
3. Generates 2D preview (`variants.svg`) and 3D SDF structures (`ligand_library/*.sdf`).

**Output:**
```
ligand_library/
├── original.sdf
├── hydroxyl_to_amine.sdf
├── hydroxyl_to_thiol.sdf
├── ...
variants.svg
```

**Usage Example:**
```bash
python3 generate_ligand_mods.py
```

---

### Step 3: Protein Preparation (`protein_preparation.py`)
Prepares the **target protein** for docking.

**Process:**
1. **Downloads** the PDB file (`4OHU.pdb`) from RCSB.
2. **Extracts chain A** and the **NAD cofactor**, excluding 2TK.
3. **Computes binding site center** using 2TK’s coordinates (for grid setup).
4. **Generates docking configuration** (`config.txt`) with center and grid size.
5. **Fixes structural issues** using **PDBFixer**, adds hydrogens.
6. **Reattaches NAD** if removed by PDBFixer.
7. **Optionally assigns AMBER charges** using PDB2PQR.
8. Copies final receptor to the `results/` directory.

**Key Outputs:**
```
4OHU_A_NAD_fixed_with_NAD.pdb
4OHU_A_NAD_fixed_with_NAD.pqr (optional)
config.txt
results/4OHU_A_NAD_fixed_with_NAD.pdb
```

**Usage Example:**
```bash
python3 protein_preparation.py
```

---

### Step 4: In-Silico Screening (`in_silico_screening.py`)
Performs docking of each ligand against the prepared receptor using **Smina**.

**Process:**
1. Verifies Smina installation.
2. Converts receptor to `.pdbqt` format via `prepare_receptor4.py` (MGLTools).
3. Iterates through each `.sdf` ligand in `ligand_library/`.
4. Runs Smina docking with parameters from `config.txt`.
5. Extracts predicted binding affinities from docking logs.

**Configuration:**
- Scoring: Vina
- Number of modes: 1
- Exhaustiveness: 8 (from config)
- Energy range: 4 kcal/mol

**Outputs:**
```
docking_results/
├── variant_name_docked.sdf
├── variant_name_docking.log
```

**Usage Example:**
```bash
python3 in_silico_screening.py
```

---

### Step 5: Report Generation (`report.py`)
Aggregates docking results and generates a ranked summary of ligand performance.

**Process:**
1. Parses each `*_docking.log` file to extract the top (mode 1) affinity value.
2. Sorts ligands by binding energy (ascending order = stronger binding).
3. Generates a human-readable report in `results/docking_ranking.txt`.
4. Copies the top-ranked docked structure to `results/`.

**Outputs:**
```
results/
├── docking_ranking.txt
├── <best_ligand>_docked.sdf
├── 4OHU_A_NAD_fixed_with_NAD.pdb
```

**Usage Example:**
```bash
python3 report.py
```

---

## Automated Execution Scripts

The workflow is organized into three shell wrappers:

| Stage | Script | Description |
|--------|---------|-------------|
| **1. Pre-run** | `pre_run.sh` | Runs protein preparation and ligand generation. |
| **2. Run** | `run.sh` | Performs virtual screening with Smina. |
| **3. Post-run** | `post_run.sh` | Generates reports and rankings. |

Example pipeline execution:
```bash
bash pre_run.sh
bash run.sh
bash post_run.sh
```

---

## Expected Runtime

| Stage | Duration (Approx.) |
|--------|--------------------|
| Ligand Extraction + Preparation | ~2–4 min |
| Ligand Modification | ~3–5 min |
| Docking (per ligand) | ~2–5 min |
| Full Library (10–20 ligands) | ~1–2 hrs |
| Report Generation | <1 min |

---

## Output Structure

```
workflow-002/
├── 4OHU.pdb
├── ligand_library/
│   ├── *.sdf
│   └── variants.svg
├── docking_results/
│   ├── *_docked.sdf
│   ├── *_docking.log
├── results/
│   ├── 4OHU_A_NAD_fixed_with_NAD.pdb
│   ├── docking_ranking.txt
│   ├── best_ligand_docked.sdf
└── config.txt
```

---

## Docking Result Ranking (Strongest Binding First)

```
Rank 1: Compound halogen_to_amine, Binding Energy: -10.10 kcal/mol
Rank 2: Compound halogen_to_hydroxyl, Binding Energy: -10.10 kcal/mol
Rank 3: Compound 2TK, Binding Energy: -9.90 kcal/mol
Rank 4: Compound original, Binding Energy: -9.80 kcal/mol
Rank 5: Compound hydroxyl_to_amine, Binding Energy: -9.60 kcal/mol
Rank 6: Compound hydroxyl_to_halogen, Binding Energy: -9.10 kcal/mol
Rank 7: Compound hydroxyl_to_thiol, Binding Energy: -8.00 kcal/mol
```

**Interpretation:**
- Ligands with more negative binding energies exhibit stronger predicted affinity.  
- `halogen_to_amine` and `halogen_to_hydroxyl` show the best potential as drug candidates with ≈ -10.1 kcal/mol binding energies.

---

## Dependencies

Installed inside the Docker environment:
- **Python 3.x**
- **RDKit**
- **Smina**
- **BioPython**
- **OpenMM** / **PDBFixer**
- **MGLTools (prepare_receptor4.py)**
- **PDB2PQR**
- **NumPy**

---

## References
- **Smina:** https://sourceforge.net/projects/smina/  
- **AutoDock Vina:** Trott & Olson, *J. Comput. Chem.* 31, 455–461 (2010).  
- **PDBFixer:** https://github.com/openmm/pdbfixer  
- **PDB2PQR:** Dolinsky et al., *Nucleic Acids Res.*, 32, W665–W667 (2004).  
