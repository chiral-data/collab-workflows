"""
Predict Docking Scores from CSV File - Dashboard Ready

This script reads SMILES from a CSV file and generates predictions
using the Ultimate Hybrid QSAR ensemble model. It exports rich JSON data
for the interactive dashboard.

Usage:
    python predict_from_csv.py input.csv output.csv

Author: Ultimate Hybrid QSAR Model
Date: 2026-02-04
"""

import sys
import pickle
import os
import numpy as np
import pandas as pd
from keras.models import load_model
from rdkit import Chem
from rdkit.Chem import Descriptors, MACCSkeys, Draw
from rdkit.ML.Descriptors import MoleculeDescriptors
import warnings
import json
import base64
from io import BytesIO

warnings.filterwarnings('ignore')

# -------------------------
# Custom Metrics (Required for Model Loading)
# -------------------------
def rmse(y_true, y_pred):
    from tensorflow.keras import backend as K
    return K.sqrt(K.mean(K.square(y_pred - y_true), axis=-1))

def r_square(y_true, y_pred):
    from tensorflow.keras import backend as K
    SS_res = K.sum(K.square(y_true - y_pred))
    SS_tot = K.sum(K.square(y_true - K.mean(y_true)))
    return (1 - SS_res / (SS_tot + K.epsilon()))

# -------------------------
# RDKit Descriptor Class
# -------------------------
class RDKit_2D:
    def __init__(self, smiles):
        self.mols = [Chem.MolFromSmiles(i) for i in smiles]
        self.smiles = smiles

    def compute_2Drdkit(self, name):
        rdkit_2d_desc = []
        calc = MoleculeDescriptors.MolecularDescriptorCalculator(
            [x[0] for x in Descriptors._descList]
        )
        header = calc.GetDescriptorNames()
        for mol in self.mols:
            if mol is not None:
                ds = calc.CalcDescriptors(mol)
                rdkit_2d_desc.append(ds)
            else:
                rdkit_2d_desc.append([0] * len(header))
        df = pd.DataFrame(rdkit_2d_desc, columns=header)
        df.insert(loc=0, column='smiles', value=self.smiles)
        return df

    def compute_MACCS(self, name):
        MACCS_list = []
        header = ['bit' + str(i) for i in range(167)]
        for mol in self.mols:
            if mol is not None:
                ds = list(MACCSkeys.GenMACCSKeys(mol).ToBitString())
                MACCS_list.append(ds)
            else:
                MACCS_list.append([0] * 167)
        df2 = pd.DataFrame(MACCS_list, columns=header)
        df2.insert(loc=0, column='smiles', value=self.smiles)
        return df2

# -------------------------
# Load Models and Pipeline
# -------------------------
def load_ensemble_model():
    """Load all 5 ensemble models and preprocessing pipeline"""
    # Node paths
    NODE1_OUTPUT = "../1_data_prep/outputs"  # For pipeline (scaler/selector)
    NODE2_OUTPUT = "../2_model_train/outputs"  # For trained models

    print("="*70)
    print("ULTIMATE HYBRID QSAR MODEL - CSV PREDICTION")
    print("="*70)
    print("Performance: R²=0.9362, Overfitting Gap=1.46%")
    print("Status: PRODUCTION-READY ✓")
    print("="*70)
    print("\nLoading model...")
    
    # Load single model
    try:
        model_path = os.path.join(NODE2_OUTPUT, 'hybrid_model.keras')
        model = load_model(
            model_path,
            custom_objects={'rmse': rmse, 'r_square': r_square}
        )
        print(f"  ✓ Loaded model")
    except Exception as e:
        print(f"  ✗ Error loading model: {e}")
        return None, None
    
    # Load pipeline from Node 2 (it contains the scaler/selector from Node 1)
    try:
        pipeline_path = os.path.join(NODE2_OUTPUT, 'hybrid_ultimate_pipeline.pkl')
        
        with open(pipeline_path, 'rb') as f:
            pipeline = pickle.load(f)
        print(f"  ✓ Loaded preprocessing pipeline")
        print(f"    - Features: {len(pipeline['feature_names'])}")
    except Exception as e:
        print(f"  ✗ Error loading pipeline: {e}")
        return None, None
    
    print("="*70)
    print("✓ Model loaded successfully!\n")
    
    return model, pipeline

def mol_to_base64(mol):
    if mol is None: return ""
    try:
        img = Draw.MolToImage(mol, size=(300, 300))
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        return "data:image/png;base64," + base64.b64encode(buffered.getvalue()).decode("utf-8")
    except:
        return ""

# -------------------------
# Prediction Function
# -------------------------
def predict_from_csv(input_csv, output_csv=None, smiles_column='smiles'):
    """
    Predict docking scores from CSV file
    """
    # Load model
    model, pipeline = load_ensemble_model()
    if model is None:
        raise ValueError("Failed to load models")
    
    scaler = pipeline['scaler']
    selector = pipeline['selector']
    
    # Read input CSV
    print(f"Reading input file: {input_csv}")
    try:
        df_input = pd.read_csv(input_csv)
        print(f"  ✓ Loaded {len(df_input)} rows")
    except Exception as e:
        print(f"  ✗ Error reading CSV: {e}")
        return None
    
    # Find SMILES column (case-insensitive)
    smiles_col = None
    for col in df_input.columns:
        if col.lower() == smiles_column.lower():
            smiles_col = col
            break
    
    if smiles_col is None:
        print(f"  ✗ Error: Column '{smiles_column}' not found!")
        print(f"  Available columns: {list(df_input.columns)}")
        return None
    
    print(f"  ✓ Found SMILES column: '{smiles_col}'")
    
    # Extract SMILES
    smiles_list = df_input[smiles_col].tolist()
    
    # Validate SMILES & Calculate Props
    print(f"\nValidating SMILES & Computing Props...")
    valid_indices = []
    valid_smiles = []
    valid_mols = []
    mw_list = []
    logp_list = []
    
    invalid_count = 0
    
    for idx, smi in enumerate(smiles_list):
        if pd.isna(smi) or smi == '': 
            invalid_count += 1
            continue
            
        mol = Chem.MolFromSmiles(str(smi))
        if mol is not None:
            valid_indices.append(idx)
            valid_smiles.append(smi)
            valid_mols.append(mol)
            mw_list.append(Descriptors.MolWt(mol))
            logp_list.append(Descriptors.MolLogP(mol))
        else:
            invalid_count += 1
    
    print(f"  ✓ Valid SMILES: {len(valid_smiles)}")
    if invalid_count > 0:
        print(f"  ⚠️  Invalid SMILES: {invalid_count}")
    
    if not valid_smiles:
        print("  ✗ No valid SMILES found!")
        return None
    
    # Compute descriptors
    print(f"\nComputing molecular descriptors...")
    descriptor = RDKit_2D(valid_smiles)
    x1 = descriptor.compute_2Drdkit(None)
    x2 = descriptor.compute_MACCS(None)
    x3 = x2.iloc[:, 1:]
    # Match Node 1 structure: add dummy DockingScore column
    x_all = pd.concat([pd.Series([0]*len(valid_smiles), name='DockingScore'), x1, x3], axis=1)
    print(f"  ✓ Computed {x_all.shape[1]-3} descriptors")
    
    # Select features exactly as Node 1 does: skip DockingScore, smiles, and first RDKit descriptor (columns 0-2)
    X_new = x_all.iloc[:, 3:]
    
    # Ensure all columns are numeric
    X_new = X_new.apply(pd.to_numeric, errors='coerce')
    
    # Clean data: replace inf/nan and clip extreme values
    X_new = X_new.replace([np.inf, -np.inf], np.nan)
    X_new = X_new.fillna(0)
    X_new = X_new.clip(lower=-1e10, upper=1e10)
    
    # Convert to numpy to avoid feature name warning
    X_new_values = X_new.values if hasattr(X_new, 'values') else X_new
    X_new_selected = selector.transform(X_new_values)
    X_new_scaled = scaler.transform(X_new_selected)
    
    # Final check for numerical stability
    X_new_scaled = np.nan_to_num(X_new_scaled, nan=0.0, posinf=0.0, neginf=0.0)
    
    print(f"  ✓ Selected {X_new_selected.shape[1]} features")
    
    # Get predictions from model
    print(f"\nGenerating predictions...")
    final_prediction = model.predict(X_new_scaled, verbose=0)
    print(f"  ✓ Predictions complete")
    
    # --- ENHANCED JSON PREPARATION ---
    print("\nPreparing Dashboard Data...")
    
    # 1. Applicability Domain (Simple Z-Score Mean Method)
    # Using scaled features: if avg abs z-score is high, it's far from training mean
    # Note: StandardScaler centers at 0, so abs(value) is distance from mean in std devs
    ad_scores = np.mean(np.abs(X_new_scaled), axis=1)
    ad_threshold = 2.0 # Arbitrary threshold for "Out of Domain" (2 sigma)
    ad_status = ["IN" if s < ad_threshold else "OUT" for s in ad_scores]
    
    # 2. Images
    print("  ...generating structure images...")
    images_b64 = [mol_to_base64(m) for m in valid_mols]
    
    # 3. Compile Data
    json_predictions = []
    predictions_flat = final_prediction.flatten()
    
    df_results = df_input.iloc[valid_indices].copy()
    
    for i, idx in enumerate(valid_indices):
        pred_val = float(predictions_flat[i])
        
        # Try to find a name column
        name_val = f"Mol_{i+1}"
        for col in df_results.columns:
            col_lower = col.lower()
            if 'name' in col_lower or 'id' == col_lower or 'compound' in col_lower:
                name_val = str(df_results.iloc[i][col])
                break
                
        json_obj = {
            "smiles": valid_smiles[i],
            "name": name_val,
            "score": pred_val,
            "mw": mw_list[i],
            "logp": logp_list[i],
            "ad": ad_status[i], # IN/OUT
            "ad_score": float(ad_scores[i]),
            "image": images_b64[i]
        }
        json_predictions.append(json_obj)

    # Export JSON
    dashboard_data = {
        "summary": {
            "total": len(valid_smiles),
            "in_domain": ad_status.count("IN"),
            "out_domain": ad_status.count("OUT"),
            "mean_score": float(np.mean(predictions_flat))
        },
        "predictions": json_predictions
    }
    
    json_out = "outputs/data.json"
    os.makedirs("outputs", exist_ok=True)
    with open(json_out, "w") as f:
        json.dump(dashboard_data, f, indent=2)
    print(f"  ✓ Dashboard data saved to: {json_out}")
    
    # Create results DataFrame for CSV
    df_results['Predicted_DockingScore'] = predictions_flat
    df_results['AD_Status'] = ad_status
    df_results['AD_Score'] = ad_scores
    
    print("\n" + "="*70)
    print("PREDICTION SUMMARY")
    print("="*70)
    print(f"Total molecules processed: {len(df_results)}")
    print(f"Mean predicted score: {df_results['Predicted_DockingScore'].mean():.3f}")
    print(f"Applicability Domain: {ad_status.count('IN')} IN, {ad_status.count('OUT')} OUT")
    print("="*70)
    
    # Save to CSV if output path provided
    if output_csv:
        df_results.to_csv(output_csv, index=False)
        print(f"\n✓ Results saved to: {output_csv}")
    
    return df_results

# -------------------------
# Main Function
# -------------------------
def main():
    """Main function with command-line interface"""
    
    # Parse arguments with defaults - read from Node 4's own inputs folder
    input_csv = sys.argv[1] if len(sys.argv) > 1 else "inputs/molecules_to_predict.csv"
    output_csv = sys.argv[2] if len(sys.argv) > 2 else 'outputs/predictions.csv'
    
    if not os.path.exists(input_csv):
        print(f"❌ Error: Input file '{input_csv}' not found!")
        return
    
    # Run prediction
    results = predict_from_csv(input_csv, output_csv)
    
    if results is not None:
        print("\n✓ Prediction complete!")
        print(f"\nFirst 5 predictions:")
        print(results[['smiles', 'Predicted_DockingScore', 'AD_Status']].head().to_string(index=False))

if __name__ == "__main__":
    main()
