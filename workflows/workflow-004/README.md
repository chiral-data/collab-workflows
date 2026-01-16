# 🧬 AutoDock Vina Virtual Screening Workflow

**Protein–Ligand Docking with Pocket Prediction and Interactive Visualization**  

A fully modular, production-ready workflow for protein–ligand virtual screening using AutoDock Vina with integrated P2Rank pocket prediction, executed on the Silva workflow platform.

---

## 📌 Overview

This workflow provides an **end-to-end, reproducible, and flexible pipeline** for computational drug discovery:

- **Six containerized nodes** with independent execution and beautiful HTML reports
- **Integrated pocket prediction** using P2Rank for intelligent grid box placement
- Automated **ligand ranking** by binding affinity with interaction analysis
- **Silva-compatible** with seamless data propagation between nodes
- **Interactive 3D visualization** using NGL.js molecular viewer

**Key Features:**

✨ **Modular Design**: Nodes run independently or sequentially  
🐳 **Containerized**: Each node runs in isolated environment  
🎨 **Rich Visualization**: Interactive HTML reports with 3D molecular viewers  
📊 **Comprehensive Analytics**: Real-time interaction analysis and pose comparison  
🔄 **Reproducible**: Full input/output tracking and metadata generation  
⚡ **GPU-Ready**: AutoDock Vina GPU acceleration support  

---

## 🧩 Workflow Structure

```
vina_workflow_final/
├── 01_Receptor_Acquisition/
│   ├── .chiral/
│   │   └── job.toml
│   ├── run.sh
│   ├── download_receptor_from_pdb.py
│   └── generate_receptor_report.py
│
├── 02_Ligand_Collection/
│   ├── .chiral/
│   │   └── job.toml
│   ├── run.sh
│   ├── download_ligands_from_pubchem.py
│   └── generate_ligand_report.py
│
├── 03_Receptor_Preparation/
│   ├── .chiral/
│   │   └── job.toml
│   ├── run.sh
│   ├── prepare_receptor_for_docking.py
│   └── generate_refinement_report.py
│
├── 04_Ligand_Preparation/
│   ├── .chiral/
│   │   └── job.toml
│   ├── run.sh
│   ├── prepare_ligands_for_docking.py
│   └── generate_ligand_refinement_report.py
│
├── 05_Pocket_Discovery/
│   ├── .chiral/
│   │   └── job.toml
│   ├── run.sh
│   ├── predict_binding_pockets.py
│   └── generate_pocket_analysis.py
│
├── 06_VS_Analytics/
│   ├── .chiral/
│   │   └── job.toml
│   ├── run.sh
│   ├── run_autodock_vina.py
│   └── generate_scientific_dashboard.py
│
├── global_params.json
└── README.md
```

Each node generates:
- **Processed data files** (PDB, PDBQT, SDF)
- **Metadata JSON** with comprehensive run information
- **Interactive HTML report** with 3D visualization

---

## 🔗 Workflow Dependency Diagram

```
┌─────────────────────┐
│ 01_Receptor         │
│    Acquisition      │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│ 03_Receptor         │
│   Preparation       │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐       ┌─────────────────────┐
│ 05_Pocket           │       │ 02_Ligand           │
│   Discovery         │       │   Collection        │
└──────────┬──────────┘       └──────────┬──────────┘
           │                             │
           │                             ▼
           │                  ┌─────────────────────┐
           │                  │ 04_Ligand           │
           │                  │   Preparation       │
           │                  └──────────┬──────────┘
           │                             │
           └──────────┬──────────────────┘
                      ▼
           ┌─────────────────────┐
           │ 06_VS_Analytics     │
           │  (Vina + Dashboard) │
           └─────────────────────┘
```

---

## 🔬 Node Descriptions

### ① Receptor Acquisition – `01_Receptor_Acquisition`

**Purpose:** Download and validate receptor structures from RCSB PDB

**Input Parameters:**
- `PARAM_PDB_ID`: PDB identifier (e.g., `5kir`, `1HSG`)

**Outputs:**
```
outputs/
├── {pdb_id}.pdb              # Downloaded structure
├── receptor_metadata.json     # Structure statistics
└── report.html                # Interactive 3D viewer
```

**Features:**
- ✅ Automatic download from RCSB PDB
- ✅ Structure validation and atom counting
- ✅ Interactive 3D visualization with multiple representations
- ✅ Chain and residue analysis

**Example:**
```bash
export PARAM_PDB_ID="5kir"
./run.sh
```

---

### ② Ligand Collection – `02_Ligand_Collection`

**Purpose:** Download ligand structures from PubChem database

**Input Parameters:**
- `PARAM_LIGAND_IDS`: JSON array of CIDs (e.g., `["3672", "2662"]`)

**Outputs:**
```
outputs/
├── {cid}.sdf                  # 3D conformers
├── ligands_metadata.json      # Download statistics
└── report.html                # Multi-ligand viewer
```

**Features:**
- ✅ Batch download from PubChem
- ✅ 3D conformer retrieval
- ✅ Individual 3D viewers for each ligand
- ✅ Download status tracking
- ✅ File size and format validation

**Example:**
```bash
export PARAM_LIGAND_IDS='["3672", "2662", "5280343"]'
./run.sh
```

---

### ③ Receptor Preparation – `03_Receptor_Preparation`

**Purpose:** Prepare receptor for docking (cleaning, protonation, PDBQT conversion)

**Input:** `{pdb_id}.pdb` from Node 1

**Outputs:**
```
outputs/
├── protein_fixed.pdb          # Cleaned PDB
├── receptor.pdbqt             # Docking-ready format
├── refined_receptor_metadata.json
└── report.html                # Preparation pipeline visualization
```

**Processing Steps:**
1. **PDBFixer**: Remove heterogens, add missing atoms
2. **Protonation**: Add hydrogens at pH 7.0
3. **OpenBabel**: Convert to PDBQT with charges
4. **Validation**: Atom count verification

**Features:**
- ✅ Automatic residue fixing
- ✅ Hydrogen addition with pH control
- ✅ PDBQT format with Gasteiger charges
- ✅ Step-by-step visualization
- ✅ Before/after statistics

---

### ④ Ligand Preparation – `04_Ligand_Preparation`

**Purpose:** Convert ligands to PDBQT format with torsion tree generation

**Input:** `*.sdf` files from Node 2

**Outputs:**
```
outputs/
├── {ligand}.pdbqt             # Docking-ready ligands
├── refined_ligands_metadata.json
└── report.html                # Multi-ligand preparation report
```

**Processing:**
- **OpenBabel** conversion with `-xh` flag (preserve hydrogens)
- **Torsion detection** for flexible docking
- **Charge assignment** using Gasteiger method

**Features:**
- ✅ Batch processing
- ✅ Rotatable bond identification
- ✅ Atom count preservation
- ✅ Individual 3D PDBQT visualization
- ✅ Conversion statistics per ligand

---

### ⑤ Pocket Discovery – `05_Pocket_Discovery`

**Purpose:** Predict binding pockets and define docking grid box

**Input:** `receptor.pdb` or `protein_fixed.pdb`

**Outputs:**
```
outputs/
├── pockets.pdb                # Predicted pocket locations
├── protein.pdb                # Copy for visualization
├── grid_config.json           # Grid box parameters
├── pocket_discovery_metadata.json
└── report.html                # Interactive pocket visualization
```

**P2Rank Analysis:**
- **Pocket ranking** by druggability score
- **Residue identification** per pocket
- **Center coordinates** calculation
- **Automatic grid box** placement on top pocket

**Grid Box Configuration:**
```json
{
  "center_x": 12.5,
  "center_y": -8.3,
  "center_z": 15.7,
  "size_x": 20.0,
  "size_y": 20.0,
  "size_z": 20.0
}
```

**Features:**
- ✅ ML-based pocket prediction (P2Rank)
- ✅ Interactive 3D pocket spheres
- ✅ Automatic grid box visualization
- ✅ Sortable pocket ranking table
- ✅ Click-to-focus pocket navigation

---

### ⑥ VS Analytics – `06_VS_Analytics`

**Purpose:** High-throughput docking with AutoDock Vina and comprehensive analysis

**Inputs:**
- `receptor.pdbqt` from Node 3
- `*.pdbqt` ligands from Node 4
- `grid_config.json` from Node 5

**Outputs:**
```
outputs/
├── receptor.pdbqt             # Copy for reference
├── {ligand}_docked.pdbqt      # Multi-pose docking results
├── {ligand}_pose{N}.pdbqt     # Individual poses
├── {ligand}_docking.log       # Vina output logs
├── virtual_screening_metadata.json
└── report.html                # Unified scientific dashboard
```

**Docking Parameters:**

| Parameter       | Environment Variable    | Default | Description                |
|-----------------|-------------------------|---------|----------------------------|
| Exhaustiveness  | `PARAM_EXHAUSTIVENESS`  | 32      | Search thoroughness        |
| Number of Modes | `PARAM_NUM_MODES`       | 10      | Max poses per ligand       |
| Energy Range    | `PARAM_ENERGY_RANGE`    | 5.0     | Energy window (kcal/mol)   |

**Dashboard Features:**

🎯 **3-Panel Layout:**
- **Left Panel**: Ranked ligand list with affinity badges
- **Center Panel**: Interactive 3D viewer with receptor + ligand
- **Right Panel**: Pose selector + interaction analysis

📊 **Real-Time Analysis:**
- Binding affinity ranking
- Distance-based interaction detection (< 4.0 Å)
- Residue contact mapping
- Pose-by-pose comparison

🔬 **Interactive Controls:**
- Surface representation toggle
- Grid box visualization
- Pose switching
- Auto-centering on binding site

**Example Usage:**
```bash
export PARAM_EXHAUSTIVENESS=32
export PARAM_NUM_MODES=10
export PARAM_ENERGY_RANGE=5.0
./run.sh
```

**Output Ranking:**
```
Rank #1: 3672 (-9.2 kcal/mol)
Rank #2: 2662 (-7.8 kcal/mol)
```

---

## 🧪 Complete Example Workflow

**Target:** HIV-1 Protease (PDB: 1HSG)  
**Ligands:** Known inhibitors from PubChem

```bash
# Node 1: Download receptor
cd 01_Receptor_Acquisition
export PARAM_PDB_ID="1hsg"
./run.sh
cd ..

# Node 2: Download ligands
cd 02_Ligand_Collection
export PARAM_LIGAND_IDS='["392622", "60823", "65028"]'
./run.sh
cd ..

# Node 3: Prepare receptor
cd 03_Receptor_Preparation
./run.sh
cd ..

# Node 4: Prepare ligands
cd 04_Ligand_Preparation
./run.sh
cd ..

# Node 5: Find binding pockets
cd 05_Pocket_Discovery
./run.sh
cd ..

# Node 6: Virtual screening
cd 06_VS_Analytics
export PARAM_EXHAUSTIVENESS=32
export PARAM_NUM_MODES=10
./run.sh
cd ..
```

**Results:** Open `06_VS_Analytics/outputs/report.html` for the complete interactive dashboard.

---

## ⚙️ Installation & Dependencies

### Container Requirements

Each node uses specific Docker containers with pre-installed dependencies:

**Node 1-2:** `python:3.9-slim`
- BioPython, requests, urllib

**Node 3-4:** `continuumio/miniconda3`
- OpenBabel, PDBFixer, OpenMM

**Node 5:** `openjdk:11`
- P2Rank 2.5.1

**Node 6:** `python:3.9-slim`
- AutoDock Vina 1.2.5

### GPU Acceleration (Optional)

For Vina GPU support:
```dockerfile
FROM nvidia/cuda:11.8.0-base-ubuntu22.04
RUN apt-get update && apt-get install -y vina-gpu
```

---

## 📊 Output Formats

### Metadata JSON Structure

**Node 1: `receptor_metadata.json`**
```json
{
  "pdb_id": "5KIR",
  "pdb_file": "5kir.pdb",
  "num_atoms": 2453,
  "num_residues": 318,
  "num_chains": 2,
  "download_timestamp": "2025-01-15T10:30:00"
}
```

**Node 5: `pocket_discovery_metadata.json`**
```json
{
  "pockets": [
    {
      "pocket_name": "pocket1",
      "rank": 1,
      "score": 18.5,
      "probability": 0.95,
      "center_x": 12.5,
      "center_y": -8.3,
      "center_z": 15.7
    }
  ],
  "grid_config": { ... }
}
```

**Node 6: `virtual_screening_metadata.json`**
```json
{
  "summary": {
    "total_compounds": 3,
    "successful_dockings": 3,
    "best_affinity": -9.2
  },
  "ligand_results": [
    {
      "ligand": "3672",
      "rank": 1,
      "best_affinity": -9.2,
      "poses": [
        {"pose": 1, "affinity": -9.2, "file": "3672_pose1.pdbqt"},
        {"pose": 2, "affinity": -8.8, "file": "3672_pose2.pdbqt"}
      ]
    }
  ]
}
```

---

## 🎨 Visualization Features

### Interactive HTML Reports

Each node generates a modern, responsive HTML report with:

**Design Elements:**
- 🌈 Gradient backgrounds (purple/blue theme)
- 📱 Responsive grid layouts
- ✨ Smooth hover animations
- 🎯 Card-based information display

**3D Visualization (NGL.js):**
- Cartoon, surface, licorice, ball-and-stick representations
- Color by chain, element, or residue
- Mouse controls: rotate, zoom, pan
- Automatic centering and reset view

**Node-Specific Features:**

| Node | Visualization                          |
|------|----------------------------------------|
| 1    | Single protein structure viewer       |
| 2    | Grid of ligand structures             |
| 3    | Before/after comparison               |
| 4    | Multi-ligand preparation status       |
| 5    | Pocket spheres + grid box overlay     |
| 6    | Unified dashboard with interaction map|

---

## 📚 References

### Software
- **AutoDock Vina** – https://github.com/ccsb-scripps/AutoDock-Vina  
- **P2Rank** – https://github.com/rdk/p2rank  
- **Open Babel** – https://github.com/openbabel/openbabel  
- **PDBFixer** – https://github.com/openmm/pdbfixer  
- **OpenMM** – https://github.com/openmm/openmm  
- **BioPython** – https://biopython.org/  
- **NGL Viewer** – https://github.com/nglviewer/ngl  

### Documentation
- **Silva Workflow Platform** – https://github.com/chiral-data/silva  
- **AutoDock Vina Manual** – https://autodock-vina.readthedocs.io/  
- **P2Rank Documentation** – https://p2rank.cz/