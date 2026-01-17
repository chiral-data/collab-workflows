import pandas as pd
import numpy as np
import os
import json
import pickle
import base64
from io import BytesIO
from tensorflow.keras.models import load_model
from tensorflow.keras import backend as K
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, MACCSkeys, Draw
from rdkit.ML.Descriptors import MoleculeDescriptors

RDLogger.DisableLog("rdApp.*")

# Custom metrics for model loading
def rmse(y_true, y_pred): return K.sqrt(K.mean(K.square(y_pred - y_true), axis=-1))
def r_square(y_true, y_pred): return 0.0 # Placeholder

def compute_descriptors(smiles_list):
    mols = [Chem.MolFromSmiles(s) for s in smiles_list]
    
    # RDKit 2D
    calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
    desc_list = [calc.CalcDescriptors(m) if m else [0]*len(calc.GetDescriptorNames()) for m in mols]
    df_rdkit = pd.DataFrame(desc_list, columns=calc.GetDescriptorNames())
    
    # MACCS
    maccs_list = [list(MACCSkeys.GenMACCSKeys(m).ToBitString()) if m else [0]*167 for m in mols]
    df_maccs = pd.DataFrame(maccs_list, columns=[f'bit{i}' for i in range(167)])
    
    # Combine
    final = pd.concat([df_rdkit, df_maccs], axis=1)
    final = final.apply(pd.to_numeric, errors='coerce').fillna(0)
    return final

def check_ad(features, ad_stats_path):
    if not os.path.exists(ad_stats_path): return ["Unknown"] * len(features)
    with open(ad_stats_path, 'r') as f: stats = json.load(f)
    
    min_v = np.array(stats['min'])
    max_v = np.array(stats['max'])
    results = []
    
    # Only check columns that match the training data
    # (Assuming features match index-wise for simplicity in this snippet)
    vals = features.values
    for row in vals:
        is_in = np.all((row >= min_v) & (row <= max_v))
        results.append("IN" if is_in else "OUT")
    return results

def run():
    print("--- Node 4: Prediction ---")
    os.makedirs("outputs", exist_ok=True)
    
    # Input SMILES
    test_smiles = [
        "CCO", 
        "CC(=O)OC1=CC=CC=C1C(=O)O", 
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    ]
    
    # Paths - Silva copies dependency outputs to inputs/ folder
    model_path = "inputs/model.h5"
    scaler_path = "inputs/scaler.pkl"
    ad_path = "inputs/ad_stats.json"

    if not os.path.exists(model_path):
        print("Model not found.")
        return

    # 1. Compute Descriptors
    print("Computing descriptors...")
    features = compute_descriptors(test_smiles)
    
    # 2. Scale
    with open(scaler_path, "rb") as f:
        scaler = pickle.load(f)
    features_scaled = scaler.transform(features)
    
    # 3. Predict
    print("Loading model and predicting...")
    model = load_model(model_path, custom_objects={'rmse': rmse, 'r_square': r_square})
    preds = model.predict(features_scaled).flatten()
    
    # 4. Applicability Domain
    ad_status = check_ad(features, ad_path)
    
    # 5. Generate Images for Report
    images = []
    for s in test_smiles:
        mol = Chem.MolFromSmiles(s)
        img_str = ""
        if mol:
            buf = BytesIO()
            Draw.MolToImage(mol, size=(200,200)).save(buf, format="PNG")
            img_str = base64.b64encode(buf.getvalue()).decode("utf-8")
        images.append(f"data:image/png;base64,{img_str}")

    # Output JSON
    results = []
    for i, s in enumerate(test_smiles):
        results.append({
            "smiles": s,
            "prediction": float(preds[i]),
            "ad_status": ad_status[i],
            "image": images[i]
        })
        
    with open("outputs/data.json", "w") as f:
        json.dump(results, f)
    
    # Save CSV
    pd.DataFrame(results).to_csv("outputs/predictions.csv", index=False)
    print("Predictions saved.")

if __name__ == "__main__":
    run()