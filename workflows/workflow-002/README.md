# Workflow-002: In-Silico Virtual Screening

## Overview

This workflow automates **virtual screening** of ligand compounds against a protein target using **SMINA**, an enhanced fork of AutoDock Vina.  
It integrates molecular structure extraction, ligand preparation, receptor preparation, docking, and result ranking into a reproducible computational pipeline managed by **Silva**.

The target system modeled here is based on **PDB ID: 5Y7J** (default), with the ligand of interest **8OL** (default).

---

## Pipeline Structure

The workflow is organized into **10 main nodes**, each containing one or more Python scripts executed sequentially:

### 1. Protein Input (`protein_input`)

Prepares the target protein structure for docking.

**Scripts:**
- **`download_pdb.py`**: Downloads the PDB file from RCSB PDB database
- **`protein_input.py`**: Verifies the downloaded protein input file

**Process:**
1. Downloads the PDB structure file based on `pdb_id` parameter
2. Verifies the downloaded file exists and is valid

**Input:**
- `pdb_id` (from `global_params.json` or environment variable `PARAM_PDB_ID`)

**Output:**
```
outputs/
└── {pdb_id}.pdb
```

**Dependencies:** None (entry point)

---

### 2. Ligand Download (`ligand_download`)

Downloads the ligand library ZIP file.

**Scripts:**
- **`download_ligands.py`**: Downloads ligand library ZIP file

**Process:**
1. Downloads a ligand library ZIP file from a specified URL
2. Copies the downloaded file to `ligand_unpack/input/` for the next step

**Input:**
- None (implicit)

**Output:**
```
outputs/
└── ligands.zip
```

**Dependencies:** None

---

### 3. Ligand Unpack (`ligand_unpack`)

Unpacks the ligand library.

**Scripts:**
- **`unpack_ligands.py`**: Unpacks the ligand library

**Process:**
1. Unpacks the library ZIP file to extract SDF files
2. Organizes ligands into `constructed_library/` directory

**Input:**
- `ligands.zip` from `ligand_download`

**Output:**
```
outputs/
└── constructed_library/
    └── *.sdf
```

**Dependencies:** `ligand_download`

---

### 4. Ligand Selection (`ligand_selection`)

Selects specific ligands from the library.

**Scripts:**
- **`ligand_selection.py`**: Selects specific ligands from the library

**Process:**
1. Reads ligands from `ligand_unpack/outputs/constructed_library/`
2. Selects specific compounds based on criteria
3. Saves selected ligands to `selected_compounds/` directory

**Input:**
- `*.sdf` from `ligand_unpack/outputs/constructed_library/`

**Output:**
```
outputs/
└── selected_compounds/
    └── *.sdf
```

**Dependencies:** `ligand_unpack`

---

### 5. True Ligand Addition (`true_ligand_addition`)

Downloads the true ligand from PDB and adds it to the workflow.

**Scripts:**
- **`true_ligand_addition.py`**: Downloads the true ligand SDF file from RCSB PDB database

**Process:**
1. Downloads the true ligand SDF file from RCSB PDB database using `ligand_name` from `global_params.json`
2. Saves the ligand as `true_ligand.sdf` in `outputs/`

**Input:**
- `pdb_id` and `ligand_name` from `global_params.json`

**Output:**
```
outputs/
└── true_ligand.sdf
```

**Dependencies:** `protein_input`, `ligand_selection`

---

### 6. Library Construction (`library_construction`)

Prepares all ligands (selected library ligands + true ligand) for docking.

**Scripts:**
- **`prepare_ligands.py`**: Prepares all ligands - adds hydrogens, assigns charges, generates 3D structures, minimizes energy using RDKit UFF (with timeout handling: 60s per ligand, 30s per step)
- **`ligand_view.py`**: Generates visualization of ligands

**Process:**
1. Collects ligands from two sources:
   - Selected compounds from `ligand_selection/outputs/selected_compounds/`
   - True ligand from `true_ligand_addition/outputs/true_ligand.sdf`
2. Prepares all ligands using Open Babel and RDKit:
   - Adds hydrogens with pH 7.4 consideration (proper protonation state)
   - Assigns Gasteiger partial charges
   - Generates 3D structures using RDKit
   - Minimizes energy using RDKit UFF force field
   - **Timeout handling**: Each ligand processing is limited to 60 seconds, with 30 seconds per step. Ligands exceeding the timeout are skipped to prevent workflow hanging.
3. Generates visualization files

**Input:**
- `*.sdf` from `ligand_selection/outputs/selected_compounds/`
- `true_ligand.sdf` from `true_ligand_addition/outputs/`

**Output:**
```
outputs/
├── selected_compounds/
│   ├── *.sdf (prepared ligands)
│   └── true_ligand.sdf
├── ligand.csv
└── variants.svg
```

**Dependencies:** `true_ligand_addition`

---

### 7. Docking Setup (`docking_setup`)

Prepares the docking configuration.

**Scripts:**
- **`ligand_loading.py`**: Extracts ligand information from PDB file and generates `ligand_list.txt`. Uses `ligand_name` from `global_params.json` to select the ligand (if multiple found, uses the first one)
- **`extract_chains.py`**: Extracts protein chains (default: all chains) and the real ligand from PDB
- **`ligand_center_identification.py`**: Generates `config.txt` for SMINA docking based on ligand center coordinates. Uses `ligand_name` from `global_params.json` to select the ligand (if multiple found, uses the first one)

**Process:**
1. Loads ligand information from the PDB file
2. Extracts protein chains (default: all chains) and the real ligand. The ligand is selected based on `ligand_name` from `global_params.json` (if multiple ligands with the same name are found, the first one is used)
3. Calculates the binding site center using the selected real ligand's coordinates from `global_params.json`
4. Generates docking configuration file (`config.txt`) with center coordinates and grid size

**Input:**
- `*.pdb` from `protein_input/outputs/`

**Output:**
```
outputs/
├── {pdb_id}_chain.pdb
├── real_ligand.sdf
├── ligand_list.txt
└── config.txt
```

**Dependencies:** `protein_input`

**Configuration File (`config.txt`):**
```
center_x = <x_coordinate>  # Calculated from ligand_name in global_params.json
center_y = <y_coordinate>  # Calculated from ligand_name in global_params.json
center_z = <z_coordinate>  # Calculated from ligand_name in global_params.json
size_x = 15.0
size_y = 15.0
size_z = 15.0
exhaustiveness = 8
num_modes = 5
energy_range = 4.0
```

---

### 8. Protein Extraction (`protein_extraction`)

Cleans up the protein structure using PDBFixer.

**Scripts:**
- **`protein_extraction.py`**: Cleans up the protein structure using PDBFixer (default: all chains)

**Process:**
1. Reads PDB file from `protein_input/outputs/`
2. Cleans the protein structure using PDBFixer (adds missing residues, adds hydrogens, fixes structural issues)
3. By default, all chains are included in the cleaned structure

**Input:**
- `*.pdb` from `protein_input/outputs/`

**Output:**
```
outputs/
└── {pdb_id}_clean.pdb (or {pdb_id}_chain_{chain_id}_clean.pdb if specific chain selected)
```

**Dependencies:** `protein_input`

---

### 9. Docking (`docking`)

Performs molecular docking of each ligand against the prepared receptor using **SMINA**.

**Scripts:**
- **`smina_screening.py`**: Executes SMINA docking for each ligand in `selected_compounds/`

**Process:**
1. Verifies SMINA installation
2. Locates the receptor PDB file (cleaned protein structure from `protein_extraction`)
3. Reads docking configuration from `config.txt` (from `docking_setup`)
4. Iterates through each ligand SDF file in `selected_compounds/` (from `library_construction`)
5. Runs SMINA docking with parameters from `config.txt`
6. Extracts predicted binding affinities from docking logs

**Configuration:**
- Scoring: Vina
- Number of modes: 1 (configurable via `config.txt`)
- Exhaustiveness: 8 (configurable via `config.txt`)
- Energy range: 4.0 kcal/mol (configurable via `config.txt`)
- Grid center and size: From `config.txt`

**Input:**
- `*.pdb` from `protein_extraction/outputs/` (cleaned protein structure)
- `selected_compounds/` from `library_construction/outputs/selected_compounds/` (ligand SDF files)
- `config.txt` from `docking_setup/outputs/` (docking configuration)

**Output:**
```
outputs/
└── docking_results/
    ├── {ligand_name}_docked.sdf
    └── {ligand_name}_log.txt
```

**Dependencies:** `library_construction`, `docking_setup`, `protein_extraction`

---

### 10. Report (`report`)

Aggregates docking results and generates a ranked summary of ligand performance.

**Scripts:**
- **`reporting.py`**: Generates docking result ranking and copies top compounds

**Process:**
1. Parses each `*_log.txt` file to extract the top (mode 1) binding affinity value
2. Sorts ligands by binding energy (ascending order = stronger binding)
3. Generates a human-readable report in `results/docking_ranking.txt`
4. Copies the top 3 ranked docked ligand structures and receptor PDB to `results/`

**Input:**
- `docking_results/` from `docking/outputs/docking_results/` (docking result files)

**Output:**
```
outputs/
└── results/
├── docking_ranking.txt
├── top1_{ligand_name}_docked.sdf
├── top2_{ligand_name}_docked.sdf
└── top3_{ligand_name}_docked.sdf
```

**Dependencies:** `docking_setup`, `docking`

---

## Execution Flow

The workflow is executed using **Silva**, which automatically manages dependencies and execution order:

```
protein_input ──┬──> docking_setup ────────┐
                │                          │
                └──> protein_extraction ───┼──> docking ──> report
                │                          │
                └──> true_ligand_addition ─┘
                                             │
ligand_download ──> ligand_unpack ──> ligand_selection ──> library_construction ──┘
```

**Execution Order:**
1. **`protein_input`** and **`ligand_download`** run first (no dependencies, can run in parallel)
2. **`ligand_unpack`** runs after `ligand_download` completes
3. **`ligand_selection`** runs after `ligand_unpack` completes
4. **`true_ligand_addition`** runs after both `protein_input` and `ligand_selection` complete
5. **`docking_setup`** and **`protein_extraction`** run after `protein_input` completes (can run in parallel)
6. **`library_construction`** runs after `true_ligand_addition` completes
7. **`docking`** runs after `library_construction`, `docking_setup`, and `protein_extraction` all complete
8. **`report`** runs last, after `docking_setup` and `docking` complete

---

## Configuration

### Global Parameters (`global_params.json`)

Located at `workflows/workflow-002/global_params.json`:

```json
{
  "pdb_id": "5Y7J",
  "ligand_name": "8OL"
}
```

**Parameters:**
- `pdb_id`: The PDB ID of the target protein structure (default: "5Y7J")
- `ligand_name`: The three-letter residue name of the ligand to extract (default: "8OL")

### Environment Variables

Parameters can also be set via environment variables (takes precedence over `global_params.json`):
- `PARAM_PDB_ID` or `PDB_ID`: PDB ID
- `PARAM_LIGAND_NAME` or `LIGAND_NAME`: Ligand name (used for downloading true ligand and selecting ligand for config generation)
- `CHAIN_ID`: Chain ID(s) to extract, comma-separated for multiple chains (e.g., "A" or "A,B"). Default: all chains (if not specified)

---

## Expected Runtime

| Stage | Duration (Approx.) |
|-------|--------------------|
| Protein Input | ~1-2 min |
| Ligand Download | ~1-2 min |
| Ligand Unpack | ~1-2 min |
| Ligand Selection | <1 min |
| True Ligand Addition | ~1 min |
| Library Construction | ~5-10 min (depends on number of ligands) |
| Docking Setup | ~2-3 min |
| Protein Extraction | ~1-2 min |
| Docking (per ligand) | ~2-5 min |
| Full Library (20-30 ligands) | ~10 min |
| Report Generation | <1 min |

**Total Runtime:** ~20-30 min for a typical screening of 20-30 ligands

---

## Output Structure

```
workflow-002/
├── global_params.json
├── protein_input/
│   └── outputs/
│       └── {pdb_id}.pdb
├── ligand_download/
│   └── outputs/
│       └── ligands.zip
├── ligand_unpack/
│   └── outputs/
│       └── constructed_library/
│           └── *.sdf
├── ligand_selection/
│   └── outputs/
│       └── selected_compounds/
│           └── *.sdf
├── true_ligand_addition/
│   └── outputs/
│       └── true_ligand.sdf
├── library_construction/
│   └── outputs/
│       ├── selected_compounds/
│       │   ├── *.sdf (prepared ligands)
│       │   └── true_ligand.sdf
│       ├── ligand.csv
│       └── variants.svg
├── docking_setup/
│   └── outputs/
│       ├── {pdb_id}_chain.pdb
│       ├── real_ligand.sdf
│       ├── ligand_list.txt
│       └── config.txt
├── protein_extraction/
│   └── outputs/
│       └── {pdb_id}_clean.pdb
├── docking/
│   └── outputs/
│       └── docking_results/
│           ├── {ligand_name}_docked.sdf
│           └── {ligand_name}_log.txt
└── report/
    └── outputs/
        └── results/
            ├── docking_ranking.txt
            ├── top1_{ligand_name}_docked.sdf
            ├── top2_{ligand_name}_docked.sdf
            └── top3_{ligand_name}_docked.sdf
```

---

## Docking Result Ranking

The final report (`outputs/results/docking_ranking.txt`) contains a ranked list of ligands sorted by binding energy (stronger binding = more negative energy):

```
Docking Result Ranking (sorted by binding strength)
==================================================
Rank 1: Compound {ligand_name}, Binding energy: -10.10 kcal/mol
Rank 2: Compound {ligand_name}, Binding energy: -9.90 kcal/mol
Rank 3: Compound {ligand_name}, Binding energy: -9.80 kcal/mol
...
```

**Interpretation:**
- Ligands with more negative binding energies exhibit stronger predicted affinity
- The top-ranked compounds are the most promising drug candidates
- Typical binding energies range from -5 to -12 kcal/mol

---

## Dependencies

All dependencies are installed inside the Docker container (`chiral.sakuracr.jp/smina:2025_11_06`):

- **Python 3.12**
- **RDKit** - Cheminformatics library for ligand manipulation
- **SMINA** - Enhanced AutoDock Vina for molecular docking
- **BioPython** - Biological computation library for PDB parsing
- **OpenMM / PDBFixer** - Protein structure cleaning and preparation
- **Open Babel** - Ligand preparation (hydrogen addition with pH 7.4 consideration, charge assignment)
- **NumPy** - Numerical computing

---

## Running the Workflow

### Using Silva (Recommended)

The workflow is designed to run in **Silva**, which automatically:
- Manages Docker containers
- Resolves node dependencies
- Mounts input/output files between nodes
- Executes nodes in the correct order

1. Configure `global_params.json` with your target PDB ID and ligand name
2. Run the workflow in Silva TUI
3. Monitor execution progress
4. Review results in `report/outputs/results/`

### Manual Execution (Testing)

For testing individual nodes, you can run them manually:

```bash
# Activate conda environment
conda activate in_silico

# Run individual nodes
cd protein_input
bash run.sh

cd ../ligand_download
bash run.sh

cd ../ligand_unpack
bash run.sh

cd ../ligand_selection
bash run.sh

cd ../true_ligand_addition
bash run.sh

cd ../library_construction
bash run.sh

cd ../docking_setup
bash run.sh

cd ../protein_extraction
bash run.sh

cd ../docking
bash run.sh

cd ../report
bash run.sh
```

---

## Notes

1. **File Mounting**: Silva mounts outputs from upstream nodes to the current node's working directory. Scripts check both the root directory and `outputs/` subdirectory to handle different mounting patterns.

2. **Ligand Preparation**: The `prepare_ligands.py` script prepares all ligands (library + true_ligand) using Open Babel and RDKit. It adds hydrogens with pH 7.4 consideration for proper protonation state, assigns Gasteiger charges, generates 3D structures, and minimizes energy. **Timeout handling**: Each ligand processing is limited to 60 seconds (30 seconds per step) to prevent workflow hanging. Ligands exceeding the timeout are automatically skipped. Failed preparations are logged but don't stop the workflow.
   Energy minimization is performed using RDKit UFF.

3. **True Ligand Download**: The `true_ligand_addition.py` script downloads the true ligand SDF file directly from RCSB PDB database using the `ligand_name` specified in `global_params.json`. This eliminates the dependency on `docking_setup` for the true ligand file.

4. **Ligand Selection**: When multiple ligands with the same name are found in the PDB file, the workflow automatically selects the first one. The ligand selection is based on `ligand_name` from `global_params.json`, which is used both for downloading the true ligand and for generating the docking configuration.

5. **Chain Selection**: By default, all protein chains are included in the docking preparation. To use specific chains, set the `CHAIN_ID` environment variable (e.g., `export CHAIN_ID=A` or `export CHAIN_ID=A,B`).

6. **Protein Cleaning**: PDBFixer may remove cofactors (like NAD) during cleaning. The workflow handles this by preserving important cofactors when possible.

7. **Docking Configuration**: The `config.txt` file is generated automatically based on the real ligand's coordinates from `global_params.json`. The ligand is selected based on `ligand_name`, and if multiple ligands with the same name are found, the first one is used. You can manually edit `config.txt` to adjust the docking search space.

8. **SMINA Path**: In Docker containers, SMINA is located at `/usr/local/bin/smina`. For local testing on macOS, a local `smina.osx.12` binary is used if available.

---

## References

- **SMINA:** https://sourceforge.net/projects/smina/  
- **AutoDock Vina:** Trott & Olson, *J. Comput. Chem.* 31, 455–461 (2010).  
- **PDBFixer:** https://github.com/openmm/pdbfixer  
- **Open Babel:** http://openbabel.org/  
- **RDKit:** https://www.rdkit.org/  
- **BioPython:** https://biopython.org/
