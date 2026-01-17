# 🧪 QSAR Modeling & Prediction Workflow

**Machine Learning-Based Drug Discovery with Premium Interactive Visualization**

A fully modular, production-ready pipeline for Quantitative Structure-Activity Relationship (QSAR) analysis. This workflow automates data ingestion, feature engineering, deep learning model training, and binding affinity prediction, all wrapped in a standardized containerized environment.

---

## 📌 Overview

This workflow provides an **end-to-end, reproducible, and flexible pipeline** for computational chemistry and drug design:

- **Four Standardized Nodes** with independent execution and state-of-the-art reporting.
- **Advanced Machine Learning** using TensorFlow/Keras Neural Networks.
- **Premium "008-Style" Visualizations** for deep data insights.
- **Dual 2D/3D Visualization** for molecular structure analysis.
- **Reproducible Science** with comprehensive JSON metadata tracking.

**Key Features:**

✨ **Modular Design**: Run the full pipeline or individual stages independently.  
🐳 **Containerized**: Fully Dockerized environment (RDKit, TensorFlow, Scikit-learn).  
🎨 **Rich Visualization**: Interactive Heatmaps, PCA Scatter Plots, and 3D Molecular Viewers.  
📊 **Deep Analytics**: Automatic calculation of R², RMSE, and Applicability Domain (AD).  
🔍 **Interactive dashboards**: Client-side rendered HTML reports with Plotly.js and NGL.js.

---

## 🧩 Workflow Structure

```text
QSAR-workflow-script/
├── 01_Data_Preparation/
│   ├── .chiral/
│   │   └── job.toml              # Node Configuration
│   ├── run_data_prep.sh          # Entrypoint
│   ├── run_data_prep.py          # Logic: Descriptors & MACCS Keys
│   └── generate_data_prep_report.py # Viz: Heatmaps & Data Grid
│
├── 02_Feature_Engineering/
│   ├── .chiral/
│   │   └── job.toml
│   ├── run_feature_eng.sh
│   ├── run_feature_eng.py        # Logic: Scaling, PCA, Split
│   └── generate_feature_eng_report.py # Viz: PCA Outliers
│
├── 03_Model_Training/
│   ├── .chiral/
│   │   └── job.toml
│   ├── run_training.sh
│   ├── run_training.py           # Logic: Neural Network Training
│   └── generate_training_report.py # Viz: Loss Curves & Parity Plots
│
├── 04_Prediction/
│   ├── .chiral/
│   │   └── job.toml
│   ├── run_prediction.sh
│   ├── run_prediction.py         # Logic: Inference & AD Check
│   └── generate_prediction_report.py # Viz: Dual 2D/3D Cards
│
├── Dockerfile                    # Container Definition
├── workflow.toml                 # Master DAG Definition
├── run_workflow.sh               # Master Execution Script
└── build_image.ps1               # Build Helper
```

Each node generates:
- **Processed data artifacts** (.csv, .npz, .h5, .pkl)
- **Metadata JSON** (`outputs/data.json`)
- **Interactive HTML report** (`report.html`)

---

## 🔗 Workflow Dependency Diagram

```text
┌─────────────────────┐
│ 01_Data_Preparation │
│   (Descriptors)     │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│ 02_Feature_Engineering
│   (PCA & Scaling)   │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│ 03_Model_Training   │
│   (Deep Learning)   │
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│ 04_Prediction       │
│   (Inference & AD)  │
└─────────────────────┘
```

---

## 🔬 Node Descriptions

### ① Data Preparation – `01_Data_Preparation`

**Purpose:** Ingest raw SMILES data, compute physicochemical descriptors, and generate MACCS keys.

**Input:** `data/raw_data.csv` (SMILES, Activity)

**Outputs:**
```text
outputs/
├── data.json                  # Metadata, Correlation Matrix, Sample Images
└── report.html                # Interactive Heatmap & Data Grid
```

**Features:**
- ✅ **RDKit Integration**: Computes MolWt, LogP, TPSA, H-Donors, H-Acceptors.
- ✅ **MACCS Keys**: Generates 167-bit structural fingerprints.
- ✅ **Correlation Heatmap**: "008-Style" interactive matrix with threshold filtering.
- ✅ **Data Preview**: Vina-style grid of training molecule structures.

---

### ② Feature Engineering – `02_Feature_Engineering`

**Purpose:** Clean data, normalize features (StandardScaler), and perform dimensionality reduction (PCA).

**Input:** `descriptors.csv` from Node 1

**Outputs:**
```text
outputs/
├── processed_data.npz         # Train/Test arrays
├── scaler.pkl                 # Saved StandardScaler
├── ad_stats.json              # Applicability Domain Min/Max
├── data.json                  # PCA Coordinates
└── report.html                # Interactive PCA Plot
```

**Features:**
- ✅ **Outlier Detection**: Identifies outliers via Isolation Forest or Variance logic.
- ✅ **Normalization**: Z-score standardization for Neural Network stability.
- ✅ **PCA Visualization**: Interactive 2D scatter plot showing data distribution.
- ✅ **Applicability Domain**: Calculates feature ranges for future AD checks.

---

### ③ Model Training – `03_Model_Training`

**Purpose:** Train a Deep Neural Network (DNN) to predict binding affinity (pIC50).

**Inputs:** `processed_data.npz` from Node 2

**Outputs:**
```text
outputs/
├── model.h5                   # Trained Keras Model
├── data.json                  # Training History & Metrics
└── report.html                # Learning Curves
```

**Model Architecture:**
- **Input Layer**: Matches feature dimension.
- **Hidden Layers**: Dense (64, 32) with ReLU activation and Dropout (0.2).
- **Output Layer**: Linear (Regression).
- **Optimizer**: Adam (lr=0.001).

**Features:**
- ✅ **Live Metrics**: Tracks MSE, MAE, and R² per epoch.
- ✅ **Parity Plot**: Interactive Actual vs. Predicted scatter plot.
- ✅ **Loss Curves**: Zoomable training/validation loss history.

---

### ④ Prediction – `04_Prediction`

**Purpose:** Predict affinity for new molecules and visualize results in 2D/3D.

**Inputs:** 
- `model.h5` and `scaler.pkl` from Node 3.
- `ad_stats.json` from Node 2.
- New SMILES strings.

**Outputs:**
```text
outputs/
├── predictions.csv            # Results Table
├── data.json                  # Predictions + 3D PDB Blocks + Images
└── report.html                # Dual 2D/3D Dashboard
```

**Features:**
- ✅ **Dual Visualization**: Hybrid cards with **toggle switch** for 2D Image / 3D Model.
- ✅ **NGL Viewer**: Fully interactive 3D molecule inspection.
- ✅ **Applicability Domain**: Flags molecules outside of training distribution.
- ✅ **Switchable Layout**: Toggle between "Grid View" and "Table View".

---

## 🎨 Visualization Features

### Interactive HTML Reports ("Premium Style")

Each node generates a responsive, high-aesthetic HTML dashboard:

**Design Elements:**
- 🌈 **Modern Typography**: Inter font family.
- 📱 **Responsive Grid**: Adapts to any screen size.
- ✨ **Interactive Charts**: Plotly.js for zooming, panning, and tooltips.
- 🧪 **Chemical Intelligence**: RDKit and NGL.js integration.

| Node | Visualization | Key Interaction |
|------|---------------|-----------------|
| **01** | **Correlation Heatmap** | Slider to filter weak correlations; Search to highlight descriptors. |
| **02** | **PCA Scatter** | Zoom into clusters; Hover to see sample indices. |
| **03** | **Learning Curves** | Toggle train/val traces; Inspect specific epochs. |
| **04** | **Structure Grid** | **2D/3D Toggle**: Switch between static image and rotating 3D model. |

---

## ⚙️ Installation & Usage

### 1. Build the Environment
This pipeline uses a unified Docker image (`qsar-workflow`).
```powershell
./build_image.ps1
```

### 2. Run the Workflow
Execute the master script to run the full pipeline:
```bash
./run_workflow.sh
```

### 3. Run Individual Nodes
You can debug or re-run specific stages:
```bash
cd 02_Feature_Engineering
bash run_feature_eng.sh
```

---

## 📊 Output Formats

### 4. Prediction Metadata (`outputs/data.json`)
```json
{
  "smiles": ["CC(=O)Oc1ccccc1C(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"],
  "predictions": [4.5, 6.2],
  "ad_results": ["IN", "OUT"],
  "images": ["data:image/png;base64...", "..."],
  "structures_3d": ["ATOM      1  C   UNL     1...", "..."]
}
```

---

## 📚 References

- **RDKit**: Open-Source Chemoinformatics. https://www.rdkit.org
- **TensorFlow**: End-to-End Machine Learning Platform. https://www.tensorflow.org
- **Plotly**: Interactive Graphing Library. https://plotly.com
- **NGL Viewer**: WebGL Protein/Molecule Viewer. http://nglviewer.org
- **Scikit-learn**: Machine Learning in Python. https://scikit-learn.org

---
*Created for the Chiral Blueprint Migration Project.*
