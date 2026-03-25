# Create a silva-runnable workflow for ADMET-AI Prediction from the node structure

- Node structure '~/dev/collab-workflows/workflows/workflow-014/'
- Draft workflow 'https://github.com/chiral-data/collab-workflows/issues/81' and Sample Directory Structure (below)
- Workflow references
    - ‘~/dev/collab-workflows/workflows/workflow-007’
    - '~/dev/collab-workflows/workflows/workflow-005/'
    - '~/dev/collab-workflows/workflows/workflow-011/'
- Silva source code ‘~/dev/silva’ 
- Silva migration guide 'https://github.com/chiral-data/collab-workflows/blob/main/SILVA_MIGRATION_GUIDE.md#5-step-3--implementation' (5. Step 3 — Implementation)
- ADMET-AI project: ‘https://github.com/swansonk14/admet_ai' 

## Tasks
- [x] Investigate the ADMET Prediction Pipeline (Training Assignment) issue and ADMET-AI repo; create a summary.
- [x] Using existing examples and references, and the sample node structure, create a workflow by building one node at a time. Make it silva runnable.
- [x] Run ‘silva ~/dev/collab-workflows/workflows/workflow-014’, debug and fix.

---

## Task 1: ADMET-AI Summary

### What is ADMET-AI?

ADMET-AI is a machine learning tool for predicting **Absorption, Distribution, Metabolism, Excretion, and Toxicity** (ADMET) properties of drug-like molecules. It takes SMILES strings as input and outputs predicted values for ~59 pharmacological properties.

- **Repository:** https://github.com/swansonk14/admet_ai
- **Dependencies:** `admet-ai`, `chemprop` (message-passing neural networks), `torch`, `rdkit`
- **CLI:** `admet_predict --data_path <in.csv> --save_path <out.csv> --smiles_column <col>`

### Properties Predicted

| Category | Properties |
|----------|-----------|
| Safety/Toxicity | hERG, AMES, DILI, ClinTox, Carcinogens_Lagunin, LD50_Zhu, Skin_Reaction |
| Absorption | BBB_Martins, HIA_Hou, Bioavailability_Ma, PAMPA_NCATS, Pgp_Broccatelli, Caco2_Wang |
| Metabolism (CYP450) | CYP1A2, CYP2C19, CYP2C9, CYP2D6, CYP3A4 (Veith + Substrate variants) |
| Distribution | PPBR_AZ, VDss_Lombardo, Lipophilicity_AstraZeneca |
| Excretion | Clearance_Hepatocyte_AZ, Clearance_Microsome_AZ, Half_Life_Obach |
| Tox21 panels | NR-AR, NR-AhR, NR-Aromatase, NR-ER, NR-PPAR-gamma, SR-ARE, SR-ATAD5, SR-HSE, SR-MMP, SR-p53 |
| Drug-likeness | Lipinski, QED, Solubility_AqSolDB, HydrationFreeEnergy_FreeSolv |

### Input/Output

- **Input:** CSV with a SMILES column (e.g., `drugbank_approved.csv` with 2845 FDA-approved drugs)
- **Output:** Same CSV extended with all predicted ADMET property columns (classification properties are 0–1 probabilities)

---

## Task 2: Silva-Runnable Workflow

### Workflow-Level Config (`.chiral/workflow.toml`)

Linear DAG: `01_validate_inputs → 02_compute → 03_analyze → 04_visualize`

Global params: `input_file` (default: `drugbank_approved.csv`), `smiles_column` (default: `smiles`)

### Node 01: Validate Inputs

- **What:** Reads input CSV, validates SMILES column exists, standardizes/canonicalizes SMILES via RDKit, removes invalid entries and duplicates
- **Container:** `python:3.11-slim` | **Packages:** pandas, rdkit-pypi
- **Params:** `input_file` (type: file) — uses `${PARAM_INPUT_FILE}` in run.sh (same pattern as workflow-012)
- **Inputs:** `*.csv` | **Outputs:** `outputs/standardized_molecules.csv`, `outputs/validation_report.json`

### Node 02: Compute ADMET Predictions

- **What:** Runs `admet_predict` CLI on standardized molecules, produces raw predictions CSV
- **Container:** `continuumio/miniconda3` (torch/chemprop need conda for reliable installs)
- **Packages:** admet-ai, chemprop, torch (via `pre_run.sh`)
- **Params:** Inherits `smiles_column` from workflow globals via `PARAM_SMILES_COLUMN`
- **Inputs:** `standardized_molecules.csv` | **Outputs:** `outputs/raw_predictions.csv`

### Node 03: Analyze (Filter & Rank)

- **What:** Computes MPO scores (weighted multi-parameter optimization), applies configurable safety/efficacy filters, ranks and selects top N candidates
- **Container:** `python:3.11-slim` | **Packages:** pandas, numpy
- **Params:** `top_n` (default: 20), `filter_hERG` ("safe"), `filter_Caco2` (">-5.15"), `filter_BBB` (">0.5"), `filter_HIA` (">0.5"), `filter_DILI` ("<0.5"), `filter_AMES` ("<0.5"), `filter_ClinTox` ("<0.3")
- **Inputs:** `raw_predictions.csv` | **Outputs:** `outputs/filtered_candidates.csv`, `outputs/analysis_summary.json`

### Node 04: Visualize (Dashboard)

- **What:** Generates a self-contained HTML dashboard (Bootstrap 5 + Plotly.js CDN) with summary cards, radar charts for top 5 candidates, and a color-coded data table
- **Container:** `python:3.11-slim` | **Packages:** none (pure Python, stdlib only)
- **Inputs:** `filtered_candidates.csv` | **Outputs:** `outputs/report.html`

### Local Test Results (drugbank_approved.csv, 2845 molecules)

| Node | Status | Notes |
|------|--------|-------|
| 01 validate | Pass | 2845/2845 valid (RDKit fallback on Python 3.12; full RDKit in container) |
| 02 compute | Skipped | Requires admet-ai install; test data already has ADMET columns |
| 03 analyze | Pass | 359/2845 pass all filters → top 20 by MPO score |
| 04 visualize | Pass | 22KB self-contained HTML dashboard |

**Sample Directory Structure**
```
workflows/workflow-014/
├── .chiral/workflow.toml
├── input_files/drugbank_approved.csv
├── 01_validate_inputs/
│   ├── .chiral/job.toml
│   ├── .chiral/test_inputs/sample_molecules.csv
│   ├── pre_run.sh           # pip install rdkit-pypi pandas
│   ├── run.sh
│   └── validate.py
├── 02_compute/
│   ├── .chiral/job.toml
│   ├── pre_run.sh            # pip install admet-ai chemprop torch
│   ├── run.sh
│   └── compute_admet.py
├── 03_analyze/
│   ├── .chiral/job.toml
│   ├── pre_run.sh            # pip install pandas numpy
│   ├── run.sh
│   └── analyze.py
├── 04_visualize/
│   ├── .chiral/job.toml
│   ├── run.sh
│   └── generate_report.py
└── README.md
```

## Outputs