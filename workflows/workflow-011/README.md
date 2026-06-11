---
doc_id: workflow-011
domain: admet
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  QSAR modeling pipeline using RDKit descriptors and MACCS keys with a Keras
  deep neural network. Performs SMILES augmentation, feature selection, model
  training with regularization, overfitting analysis, and batch prediction
  with applicability domain assessment.
tags: [qsar, deep-learning, rdkit, maccs, docking-score, drug-discovery, keras, active-learning]
---

# Workflow 011: QSAR Modeling Pipeline

Quantitative structure-activity relationship (QSAR) modeling using RDKit molecular descriptors, MACCS structural keys, and a Keras deep neural network. This workflow ingests SMILES with docking scores, trains a regularized regression model with SMILES augmentation, analyzes overfitting diagnostics, and predicts binding scores for new molecules with applicability domain assessment.

## Overview

The pipeline builds a QSAR model that predicts molecular docking scores from 2D molecular structure. Features are computed as a combination of RDKit 2D descriptors (~200 physicochemical properties) and 167-bit MACCS structural keys. The dataset is augmented via SMILES enumeration (Bjerrum, 2017) -- generating multiple valid SMILES representations per molecule to expand the training set. After feature selection with Random Forest importance, a 3-layer dense neural network with batch normalization, dropout, and L2 regularization is trained. The model is evaluated for overfitting via train/test R-squared gap analysis, and new molecules can be scored with applicability domain flags based on descriptor-space distance from the training set (Roy et al., 2015).

## When to use this workflow

Use this workflow when you have a set of molecules with known docking scores (or other continuous activity values) in CSV format and want to build a predictive model for scoring new compounds. The input requires a CSV with `smiles` and `DockingScore` columns.

For ADMET property prediction using pre-trained models (no training data needed), use workflow-014 (ADMET-AI). For physics-based molecular docking rather than ML-predicted scores, use workflow-004 (AutoDock Vina) or workflow-002/003 (Smina). For active learning-guided lead optimization with docking, use workflow-009 (FEGrow).

## Architecture and data flow

```text
                                ┌──> [3: Analyze Overfitting] ──> overfitting report
[1: Data Prep] ──> [2: Train] ─┤
                                └──> [4: Predict from CSV] ──> predictions.csv
```

Node 1 runs first. Node 2 depends on Node 1. Nodes 3 and 4 run independently after Node 2: Node 3 requires both Node 1 and Node 2 outputs; Node 4 requires only Node 2 outputs.

## Input requirements

- **Training data CSV:** A CSV file with at least two columns: `smiles` (SMILES strings) and `DockingScore` (continuous activity values). Place in `input_files/`. The default test file is `SpikeRBD_DD.csv` (SARS-CoV-2 Spike RBD docking data).
- **Prediction input (Node 4):** A plain text file with one SMILES per line, placed in the Node 4 inputs directory as `molecules_to_predict.txt`.

## Workflow nodes

### Node 1: Data Preparation

**Goal:** Compute molecular descriptors, select features, split data, and remove outliers.

**Process:**
1. Loads training CSV and augments SMILES using random enumeration (3 variants per molecule)
2. Computes RDKit 2D descriptors and 167-bit MACCS keys for each SMILES variant
3. Trains a Random Forest regressor and selects features above median importance
4. Splits data 70/30 train/test with StandardScaler normalization
5. Removes outliers using Isolation Forest (contamination=0.1) on PCA-projected descriptor space

**Scientific notes:** SMILES enumeration generates equivalent SMILES strings by randomizing atom traversal order, effectively augmenting the dataset without creating new chemistry (Bjerrum, 2017). Feature selection via Random Forest importance reduces the combined descriptor+fingerprint space (typically 350+ features) to the most predictive subset, improving model generalization and training speed. Isolation Forest identifies molecules whose descriptor profiles are anomalous in PCA space, removing potential data quality issues.

**Outputs:**
- `processed_data.pkl` -- pickled data bundle (X_train, X_test, y_train, y_test, scaler, selector, feature names)
- `data.json` -- visualization data (feature importances, PCA outlier map, distributions, correlation matrix)
- `report.html` -- interactive data preparation dashboard

### Node 2: Model Training

**Goal:** Train a regularized deep neural network to predict docking scores.

**Process:** Builds a 3-layer Keras Sequential model (192 -> 96 -> 48 units) with ReLU activation, batch normalization, dropout (0.4/0.3/0.2), and L2 regularization (0.002/0.002/0.001). Applies label smoothing (epsilon=0.05, shrinking targets toward the global mean) and Gaussian noise injection (sigma=0.02) to training features. Trains for up to 200 epochs with early stopping (patience=40 on validation loss) and learning rate reduction (factor=0.5, patience=15). Uses Adam optimizer with gradient clipping (clipnorm=1.0).

**Scientific notes:** The combination of dropout, batch normalization, and L2 regularization provides multiple layers of overfitting protection for the typically small datasets in QSAR. Label smoothing (originally a classification technique from Szegedy et al., 2016) is applied here as regression target shrinkage toward the mean, providing mild regularization. Feature noise injection during training further encourages the model to learn robust representations rather than memorizing descriptor noise.

**Outputs:**
- `hybrid_model.keras` -- trained Keras model
- `training_history.pkl` -- per-epoch loss, R-squared, learning rate
- `hybrid_ultimate_pipeline.pkl` -- preprocessing pipeline (scaler, selector, feature names)
- `data.json` -- training dashboard data (architecture, learning curves, predictions)
- `report.html` -- interactive training dashboard

### Node 3: Analyze Overfitting

**Goal:** Diagnose model overfitting by comparing train and test performance.

**Process:** Loads the trained model and preprocessed data, generates predictions on both train and test sets, and computes R-squared, RMSE, and MAE for each. Calculates the R-squared gap (|train_R2 - test_R2| x 100) and classifies overfitting severity. Reports total model parameters and data-to-parameter ratio.

**Scientific notes:** The R-squared gap between training and test sets is the primary overfitting diagnostic. Classification thresholds: no overfitting (gap < 5%), minimal (5-10%), moderate (10-20%), severe (> 20%). A high data-to-parameter ratio (> 1.0) generally indicates sufficient training data relative to model complexity.

**Outputs:**
- `data.json` -- metrics, residuals, training history for dashboard
- `report.html` -- interactive overfitting analysis dashboard

### Node 4: Predict from CSV

**Goal:** Generate docking score predictions for new molecules with applicability domain assessment.

**Process:** Reads SMILES from a text file, validates each with RDKit, computes the same descriptor+MACCS feature set as Node 1, applies the saved feature selector and scaler, and runs inference through the trained model. Assesses applicability domain using the standardization approach: computes the mean absolute z-score of scaled features for each molecule and flags compounds exceeding a threshold of 2.0 standard deviations as outside the model's domain. Generates molecule structure images and an interactive dashboard.

**Scientific notes:** The applicability domain (AD) check follows the standardization approach of Roy et al. (2015). Molecules with descriptor values far from the training set distribution (mean |z-score| > 2.0) are flagged as "OUT" because the model has not seen similar chemistry during training and predictions may be unreliable. Predictions for "OUT" compounds should be treated with caution.

**Outputs:**
- `predictions.csv` -- SMILES, predicted docking score, AD status, AD score
- `data.json` -- prediction dashboard data with molecule images
- `report.html` -- interactive prediction dashboard

## Parameters

### input_csv

- **Type:** string
- **Default:** (empty -- uses `SpikeRBD_DD.csv` from input_files)
- **Node:** 1
- **Description:** Name of the input CSV file in the input_files directory. Must contain `smiles` and `DockingScore` columns.

### epochs

- **Type:** string
- **Default:** `"200"`
- **Node:** 3
- **Description:** Number of training epochs (used in overfitting analysis context). The actual training in Node 2 uses early stopping with patience=40, so training typically stops before 200 epochs.

### batch_size

- **Type:** string
- **Default:** `"400"`
- **Node:** 3
- **Description:** Batch size for model training reference. Node 2 uses a fixed batch size of 256.

## Outputs and interpretation

### Model performance metrics

| Metric | Description |
|--------|-------------|
| Train R-squared | Coefficient of determination on training set. Values > 0.8 indicate good fit |
| Test R-squared | Coefficient of determination on held-out test set. Primary quality metric |
| RMSE | Root mean squared error in docking score units |
| R-squared gap | |Train_R2 - Test_R2| x 100. Should be < 10% for a well-generalized model |

### Overfitting verdict

| Verdict | R-squared gap | Interpretation |
|---------|--------------|----------------|
| No overfitting | < 5% | Model generalizes well |
| Minimal overfit | 5--10% | Acceptable for most applications |
| Moderate overfit | 10--20% | Consider more training data or stronger regularization |
| Severe overfit | > 20% | Model memorizes training data; predictions unreliable |

### Applicability domain status

| Status | Meaning |
|--------|---------|
| IN | Molecule's descriptors are within 2 standard deviations of training data. Prediction is reliable |
| OUT | Molecule falls outside the training descriptor space. Prediction may be unreliable |

### predictions.csv

Key columns: `smiles`, `Predicted_DockingScore` (predicted binding score), `AD_Status` (IN/OUT), `AD_Score` (mean absolute z-score; lower is more reliable).

## Quick start

### Running with Docker

```bash
docker pull ghcr.io/chiral-data/qsar-hybrid-model:v3
```

### Running on Silva

1. Select "QSAR Modeling Pipeline" from the workflow list
2. Upload a CSV with `smiles` and `DockingScore` columns
3. Run the pipeline (Nodes 1-2 train the model; Node 3 analyzes overfitting; Node 4 predicts new molecules)
4. Check the overfitting analysis dashboard before trusting predictions

### Test vs production settings

| Setting | Test (default) | Production |
|---------|---------------|------------|
| Training data | `SpikeRBD_DD.csv` (included) | Your own docking/activity data |
| SMILES augmentation | 3 variants | 3--10 variants (edit source) |
| Chemical space size | ~100 molecules (test) | Thousands of candidates |

A successful test run with the default SARS-CoV-2 dataset trains a model, produces overfitting diagnostics, and generates prediction-ready outputs.

## Troubleshooting

### Low test R-squared

If the model achieves high training R-squared but low test R-squared (gap > 20%), the dataset may be too small or chemically diverse for the model architecture. Try collecting more training data with broader chemical coverage, or use workflow-005 (a simpler QSAR pipeline) which may generalize better on small datasets.

### All predictions flagged as OUT of domain

If most new molecules are flagged as outside the applicability domain, the training set chemistry is too different from the prediction set. The model was trained on a specific chemical series and cannot reliably extrapolate to structurally dissimilar compounds.

## References

- Bjerrum, E. J. "SMILES Enumeration as Data Augmentation for Neural Network Modeling of Molecules." *arXiv:1703.07076*, 2017. DOI: https://doi.org/10.48550/arXiv.1703.07076
- Roy, K., Kar, S. & Ambure, P. "On a simple approach for determining applicability domain of QSAR models." *Chemometrics and Intelligent Laboratory Systems* 145:22--29, 2015. DOI: https://doi.org/10.1016/j.chemolab.2015.04.013
- [RDKit documentation](https://www.rdkit.org/docs/)
- [MACCS keys (RDKit)](https://www.rdkit.org/docs/source/rdkit.Chem.MACCSkeys.html)
- [scikit-learn Isolation Forest](https://scikit-learn.org/stable/modules/generated/sklearn.ensemble.IsolationForest.html)
