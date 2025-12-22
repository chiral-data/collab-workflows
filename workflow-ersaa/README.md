# 🧬 AutoDock Vina Virtual Screening Workflow

**Protein–Ligand Docking with Pocket Prediction and Automated Ranking**  

A fully modular workflow for protein–ligand virtual screening using AutoDock Vina, with optional binding pocket prediction via P2Rank, executed on the Silva workflow platform.

---

## 📌 Overview

This workflow provides an **end-to-end, reproducible, and flexible pipeline** for virtual screening:

- **Eight containerized nodes**, each with independent inputs and outputs  
- **Optional pocket prediction** using P2Rank  
- Automated **ligand ranking** by binding affinity  
- **Silva-compatible**, allowing seamless orchestration of nodes  

**Key Features:**

- Modular design: Nodes can run independently or sequentially  
- Clean environments: Each node runs in its own container  
- Reproducibility: Inputs, outputs, and parameters are fully tracked    

---

## 🧩 Workflow Layout

```
vina_workflow_final/
├── 01_ProteinInput
├── 02_LigandInput
├── 03_ProteinPreparation
├── 04_LigandPreparation
├── 05_PocketPrediction
├── 06_PocketSelection
├── 07_Docking
├── 08_Reporting
├── global_params.json
└── README.md
```

Each node is structured as follows:

```
<node_name>/
├── .chiral/
│   ├── job.toml      # Node inputs, outputs, and dependencies
│   └── node.json     # Node description
├── run.sh           # Entry point executed by Silva
├── *.py             # Node logic
└── params.json      # Optional node-specific parameters
```

---

## 🔗 Workflow Dependency Diagram

```
01_ProteinInput
        │
        ▼
03_ProteinPreparation ───────┐
        ▼                    │
05_PocketPrediction          │
        ▼                    │
06_PocketSelection           │
        ▼                    │
07_Docking ──────────────────┼──▶ 08_Reporting
        ▲                    │
04_LigandPreparation         │
        ▲                    │
02_LigandInput ──────────────┘
```

---

## 🔬 Node Descriptions

### ① Protein Input – `01_ProteinInput`
**Purpose:** Download receptor structures from RCSB PDB  
**Inputs:** Protein PDB ID (e.g., `5KIR`)  
**Outputs:** `receptor.pdb`

---

### ② Ligand Input – `02_LigandInput`
**Purpose:** Download ligand structures from PubChem  
**Inputs:** List of PubChem CIDs (e.g., `2662`, `3672`)  
**Outputs:**
```
ligands/
├── 2662.sdf
└── 3672.sdf
```

---

### ③ Protein Preparation – `03_ProteinPreparation`
**Purpose:** Clean, protonate, and convert receptor to PDBQT using MGLTools  
**Inputs:** `receptor.pdb`  
**Outputs:** `receptor.pdbqt`

---

### ④ Ligand Preparation – `04_LigandPreparation`
**Purpose:** Convert ligands to PDBQT format for docking  
**Inputs:** Ligand structure files (`.sdf`)  
**Outputs:**
```
ligands_prepared/
├── 2662.pdbqt
└── 3672.pdbqt
```

---

### ⑤ Pocket Prediction – `05_PocketPrediction`
**Purpose:** Predict potential binding pockets using P2Rank  
**Inputs:** `receptor.pdb`  
**Outputs:** Pocket prediction files (`JSON/CSV`)

---

### ⑥ Pocket Selection – `06_PocketSelection`
**Purpose:** Convert P2Rank output into AutoDock Vina grid box definitions  
**Inputs:** Pocket prediction files  
**Outputs:** Grid box center coordinates and dimensions

---

### ⑦ Docking – `07_Docking`
**Purpose:** Perform virtual screening using AutoDock Vina  
**Inputs:**
- `receptor.pdbqt`  
- `ligands in pdbqt format`  
- Grid box definition  

**Outputs:**
- Docked poses (`*.pdbqt`)  
- Vina log files (`*.log`)

**Docking Parameters (`params.json`):**
| Parameter       | Type    | Default | Description                     |
|-----------------|---------|---------|---------------------------------|
| exhaustiveness   | integer | 8       | Search thoroughness             |
| num_modes        | integer | 9       | Maximum poses per ligand        |
| energy_range     | integer | 4       | Energy window (kcal/mol)       |

---

### ⑧ Reporting – `08_Reporting`
**Purpose:** Parse docking logs and generate a ranked Excel report  
**Inputs:** Vina log files (`*.log`)  
**Outputs:** `binding_affinities.xlsx`

**Example Report:**

| Ligand | Affinity (kcal/mol) |
|--------|--------------------|
| 2662   | -8.4               |
| 3672   | -7.9               |

Lower (more negative) values indicate stronger predicted binding.

---

## 🧪 Tested Example

- **Protein:** `5KIR`  
- **Ligands:** `2662`, `3672`  
- **Docking Engine:** AutoDock Vina  
- **Pocket Prediction:** P2Rank  

After running the workflow, the **final ranked Excel file** (`binding_affinities.xlsx`) lists ligand affinities.

---

## ⚙️ Global Parameters

**`global_params.json`** contains workflow-wide configuration values shared across nodes.

---

## 📚 References

- **AutoDock Vina** – https://github.com/ccsb-scripps/AutoDock-Vina  
- **P2Rank** – https://github.com/rdk/p2rank  
- **Open Babel** – https://github.com/openbabel/openbabel  
- **MGLTools** – https://ccsb.scripps.edu/mgltools/  
- **Pandas** – https://pandas.pydata.org/  
- **OpenPyXL** – https://openpyxl.readthedocs.io/  
- **Silva Workflow Platform** – https://github.com/chiral-data/silva

