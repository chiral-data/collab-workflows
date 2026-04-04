# Workflow-002: In-Silico Virtual Screening

## Overview

This workflow performs **in-silico virtual screening** using molecular docking to identify potential drug candidates from a compound library. It uses **Smina** (a fork of AutoDock Vina with additional features) to dock compounds against a target protein and ranks them by binding affinity.

## Workflow Structure

```
workflow-002/
├── .chiral/workflow.toml      # Workflow configuration and dependencies
├── global_params.json         # Global parameters
├── README.md
├── 01-download/               # Stage 1: Download compound library
│   ├── .chiral/job.toml
│   ├── run.sh
│   └── download.py
├── 02-prepare/                # Stage 2: Prepare target protein
│   ├── .chiral/job.toml
│   ├── run.sh
│   └── protein_preparation.py
├── 03-screening/              # Stage 3: Virtual screening
│   ├── .chiral/job.toml
│   ├── run.sh
│   └── in_silico_screening.py
└── 04-report/                 # Stage 4: Generate results report
    ├── .chiral/job.toml
    ├── run.sh
    └── report.py
```

## Parameters

### Global Parameters (workflow.toml)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pdb_id` | string | "5Y7J" | PDB ID of target protein |
| `ligand_name` | string | "8OL" | Reference ligand residue name |

### Stage 3: Screening Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `test_mode` | boolean | true | Run with 11 compounds (true) or full library (false) |
| `exhaustiveness` | integer | 8 | Search exhaustiveness (higher = slower but more thorough) |

### Stage 4: Report Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `top_n` | integer | 10 | Number of top compounds to display |

## Workflow Pipeline

### 1. Download Library (`01-download`)

Downloads the compound library from Zenodo:
- **Source**: https://zenodo.org/records/17374422
- **Format**: ZIP archive containing SDF files
- **Output**: `constructed_library/` directory with ligand structures

### 2. Prepare Protein (`02-prepare`)

Prepares the target protein for molecular docking:
1. Downloads PDB structure from RCSB
2. Extracts chains A and B with reference ligand
3. Calculates binding site coordinates
4. Generates `config.txt` for docking parameters
5. Fixes structure with PDBFixer
6. Adds AMBER charges with PDB2PQR

**Outputs**:
- `{pdb_id}_AB_chains_fixed.pdb` (prepared receptor)
- `{pdb_id}_AB_chains_fixed.pqr` (with charges)
- `config.txt` (docking configuration)

### 3. Virtual Screening (`03-screening`)

Performs high-throughput molecular docking:
- **Method**: Smina with Vina scoring function
- **Configuration**: 15x15x15 A grid, configurable exhaustiveness

**Outputs**:
- `docking_results/*_docked.sdf` (docked poses)
- `docking_results/*_log.txt` (scores)

### 4. Generate Report (`04-report`)

Analyzes results and generates rankings:
1. Parses binding affinities from log files
2. Ranks compounds by binding energy
3. Copies top compound structure

**Outputs**:
- `results/docking_ranking.txt`
- `results/*_docked.sdf` (top compound)

## Dependencies

```
01-download → 02-prepare → 03-screening → 04-report
```

## Container

All stages use: `ghcr.io/chiral-data/smina:2025_10_17_v2`

## Interpreting Results

### Binding Affinity Scores

- Measured in **kcal/mol**
- **Lower (more negative) = stronger binding**
- Typical range: -5 to -12 kcal/mol
- Values < -8 kcal/mol indicate promising candidates

## References

- **Smina**: https://sourceforge.net/projects/smina/
- **AutoDock Vina**: Trott, O. & Olson, A. J. (2010). J. Comput. Chem. 31, 455-461.
- **PDBFixer**: https://github.com/openmm/pdbfixer
- **PDB2PQR**: Dolinsky et al. (2004). Nucleic Acids Res. 32, W665-W667.
