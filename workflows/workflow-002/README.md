# Workflow-002: In-Silico Virtual Screening

## Overview

This workflow performs **in-silico virtual screening** using molecular docking to identify potential drug candidates from a compound library. It uses **Smina** (a fork of AutoDock Vina with additional features) to dock compounds against a target protein and ranks them by binding affinity.

## Target Information

- **Protein**: PDB ID `5Y7J` (chains A and B)
- **Reference Ligand**: `8OL`
- **Compound Library**: Zenodo dataset with constructed drug-like molecules

## Workflow Pipeline

### 1. Data Download (`download.py`)
Downloads the compound library from Zenodo:
- **Source**: https://zenodo.org/records/17374422
- **Format**: ZIP archive containing SDF files
- **Output**: `constructed_library/` directory with ligand structures

### 2. Protein Preparation (`protein_preparation.py`)

Prepares the target protein for molecular docking:

1. **Download PDB**: Fetches structure from RCSB PDB
2. **Chain Extraction**: Extracts chains A and B along with the reference ligand
3. **Binding Site Analysis**:
   - Calculates ligand center coordinates
   - Determines optimal grid box dimensions
   - Generates `config.txt` for docking parameters
4. **Structure Cleanup**:
   - Uses PDBFixer to add missing residues and hydrogens
   - Removes non-standard residues and heterogens
5. **Charge Assignment**: Applies AMBER force field charges using PDB2PQR
6. **Output**:
   - `5Y7J_AB_chains_fixed.pdb` (prepared receptor)
   - `5Y7J_AB_chains_fixed.pqr` (receptor with charges)
   - `config.txt` (docking configuration)

### 3. Virtual Screening (`in_silico_screening.py`)

Performs high-throughput molecular docking:

- **Method**: Smina with Vina scoring function
- **Receptor**: Prepared protein structure
- **Ligands**: Compounds from the downloaded library
- **Configuration**:
  - Grid center: Based on reference ligand position
  - Grid size: 15 × 15 × 15 Å
  - Exhaustiveness: 8
  - Number of modes: 5
  - Energy range: 4 kcal/mol

**Test Mode**: By default, processes only 11 compounds (files matching `clean_drug108*.sdf`) for testing. To screen the full library, modify line 70 in `in_silico_screening.py`.

**Output**:
- `docking_results/` directory containing:
  - `*_docked.sdf`: Docking poses for each ligand
  - `*_log.txt`: Binding affinity scores and statistics

### 4. Report Generation (`report.py`)

Analyzes docking results and generates rankings:

1. **Result Parsing**: Extracts binding affinities from log files
2. **Ranking**: Sorts compounds by binding energy (lower = stronger binding)
3. **Output Files**:
   - `results/docking_ranking.txt`: Complete ranking of all screened compounds
   - `results/*_docked.sdf`: Top-ranked compound structure
   - `results/5Y7J_AB_chains_fixed.pdb`: Prepared receptor structure

## Execution

### Prerequisites

The workflow runs in a Docker container with all dependencies pre-installed:
- Python 3.x
- Smina (molecular docking software)
- BioPython, OpenMM, PDBFixer
- PDB2PQR
- RDKit (implicit dependency)

### Running the Workflow

The workflow is divided into three stages:

```bash
# Stage 1: Download and prepare protein
bash pre_run.sh

# Stage 2: Run virtual screening
bash run.sh

# Stage 3: Generate results report
bash post_run.sh
```

### Expected Runtime

- **Download + Preparation**: ~2-5 minutes
- **Screening (11 compounds)**: ~5-10 minutes
- **Screening (full library)**: Several hours depending on library size
- **Report Generation**: <1 minute

## Output Files

```
workflow-002/job-1/
├── constructed_library/        # Downloaded compound library
├── 5Y7J_AB_chains_fixed.pdb   # Prepared protein structure
├── 5Y7J_AB_chains_fixed.pqr   # Protein with AMBER charges
├── config.txt                  # Docking configuration
├── docking_results/            # Individual docking results
│   ├── *_docked.sdf           # Docked poses
│   └── *_log.txt              # Docking logs with scores
└── results/                    # Final results
    ├── docking_ranking.txt    # Ranked compounds by affinity
    ├── *_docked.sdf          # Top compound structure
    └── 5Y7J_AB_chains_fixed.pdb
```

## Interpreting Results

### Binding Affinity Scores

- Measured in **kcal/mol**
- **Lower (more negative) = stronger binding**
- Typical range: -5 to -12 kcal/mol for drug-like molecules
- Values < -8 kcal/mol generally indicate promising candidates

### Ranking File

The `docking_ranking.txt` file lists all screened compounds sorted by binding strength:

```
1位: 化合物 clean_drug5371, 結合エネルギー: -10.25 kcal/mol
2位: 化合物 clean_drug2143, 結合エネルギー: -9.87 kcal/mol
...
```

## Customization

### Screening the Full Library

Edit `in_silico_screening.py` line 70-72:

```python
# Current (test mode):
ligands = glob.glob("./constructed_library/clean_drug108*.sdf")

# Change to (full screening):
ligands = glob.glob("./constructed_library/clean_drug*.sdf")
```

### Adjusting Docking Parameters

Modify `protein_preparation.py` lines 169-174 or edit `config.txt` directly:

```
center_x = <x>           # Binding site coordinates
center_y = <y>
center_z = <z>
size_x   = 15            # Search space size (Å)
size_y   = 15
size_z   = 15
exhaustiveness = 8       # Search effort (higher = slower but more thorough)
num_modes = 5            # Number of poses to generate
energy_range = 4         # Energy window for poses (kcal/mol)
```

### Using a Different Target

1. Modify the PDB ID in `protein_preparation.py` line 53
2. Update the ligand name in line 71 if needed
3. Adjust chain selection in lines 80-87 if required

## References

- **Smina**: https://sourceforge.net/projects/smina/
- **AutoDock Vina**: Trott, O. & Olson, A. J. (2010). AutoDock Vina: improving the speed and accuracy of docking. *J. Comput. Chem.* 31, 455-461.
- **PDBFixer**: https://github.com/openmm/pdbfixer
- **PDB2PQR**: Dolinsky et al. (2004). *Nucleic Acids Res.* 32, W665-W667.

## License

See the main repository LICENSE file.
