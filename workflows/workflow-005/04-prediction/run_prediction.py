import pandas as pd
import numpy as np
import os
import json
import ast
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

def find_file(filename, search_path):
    for root, dirs, files in os.walk(search_path):
        if filename in files:
            return os.path.join(root, filename)
    return None

import re

def get_param(param_name, default_if_missing=None):
    # 1. Try Environment Variable
    val = os.environ.get(param_name.upper()) or os.environ.get(f"PARAMS_{param_name.upper()}")
    if val is not None:
        return val

    # 2. Try Parsing .chiral/job.toml
    try:
        toml_path = os.path.join(".chiral", "job.toml")
        if os.path.exists(toml_path):
            with open(toml_path, "r") as f:
                content = f.read()
            # Regex to find [params.param_name] ... default = "value"
            # Handling multiline strings or simple values
            # Find the section
            section_match = re.search(f'\\[params\\.{param_name}\\]', content, re.IGNORECASE)
            if section_match:
                # Search for default = ... after the section
                start_index = section_match.end()
                # Look for default value until next section or EOF
                # Simple extraction for quoted string or number
                default_match = re.search(r'default\s*=\s*(.+)', content[start_index:])
                if default_match:
                    raw_val = default_match.group(1).strip()
                    # Handle quotes
                    if raw_val.startswith('"') and raw_val.endswith('"'):
                         return raw_val[1:-1].replace('\\"', '"')
                    elif raw_val.startswith("'") and raw_val.endswith("'"):
                         return raw_val[1:-1]
                    else:
                         return raw_val # number or bool
    except Exception as e:
        print(f"Warning: Failed to parse job.toml for {param_name}: {e}")
        
    if default_if_missing is not None:
        return default_if_missing
        
    raise ValueError(f"Parameter {param_name} missing (Env Var and job.toml default)")

def run():
    with open("outputs/debug_env.txt", "w") as f:
        for k, v in os.environ.items():
            f.write(f"{k}={v}\n")
    
    # Check if we should predict on all data from input CSV
    predict_all = os.environ.get("PREDICT_ALL_DATA", "true").lower() == "true"
    
    test_smiles = []
    
    if predict_all:
        # Load all SMILES from original input CSV
        print("=" * 60)
        print("PREDICT_ALL_DATA=true: Loading all compounds from input dataset")
        print("=" * 60)
        
        try:
            import glob
            csv_files = glob.glob("../01-data-preparation/inputs/*.csv")
            if csv_files:
                import pandas as pd
                original_data = pd.read_csv(csv_files[0])
                # Find SMILES column (case-insensitive)
                smiles_col = next((col for col in original_data.columns if col.lower() == 'smiles'), None)
                
                if smiles_col:
                    test_smiles = original_data[smiles_col].dropna().tolist()
                    print(f"Loaded {len(test_smiles)} compounds from {csv_files[0]}")
                else:
                    print("ERROR: No 'smiles' column found in input CSV")
                    return
            else:
                print("ERROR: No CSV file found in ../01-data-preparation/inputs/")
                return
        except Exception as e:
            print(f"ERROR loading SMILES from input CSV: {e}")
            return
    else:
        # Use provided TEST_SMILES parameter
        test_smiles_raw = os.environ.get("TEST_SMILES", "").strip()
        
        if not test_smiles_raw:
            print("=" * 60)
            print("ERROR: No test SMILES provided!")
            print("=" * 60)
            print("")
            print("Set PREDICT_ALL_DATA=true to predict on all dataset compounds,")
            print("or provide TEST_SMILES as a JSON array, for example:")
            print('  TEST_SMILES=\'["CCO", "c1ccccc1", "CC(C)O"]\'')
            print("=" * 60)
            return
        
        # Parse JSON array of SMILES
        try:
            test_smiles = json.loads(test_smiles_raw)
            if not isinstance(test_smiles, list) or len(test_smiles) == 0:
                print("ERROR: TEST_SMILES must be a non-empty JSON array")
                return
        except json.JSONDecodeError as e:
            print(f"ERROR: Invalid JSON in TEST_SMILES: {e}")
            print(f"Received: {test_smiles_raw}")
            # Fallback: try comma-separated
            test_smiles = [s.strip() for s in test_smiles_raw.replace(',', '\n').split('\n') if s.strip()]
    
    print(f"\nPredicting on {len(test_smiles)} compounds...")

    
    # Silva shared workspace: read directly from previous nodes' outputs
    model_path = os.path.join("..", "03-model-training", "outputs", "model.h5")
    if not os.path.exists(model_path): model_path = find_file("model.h5", "..")
    
    scaler_path = os.path.join("..", "02-feature-engineering", "outputs", "scaler.pkl")
    if not os.path.exists(scaler_path): scaler_path = find_file("scaler.pkl", "..")
    
    ad_path = os.path.join("..", "02-feature-engineering", "outputs", "ad_stats.json")
    if not os.path.exists(ad_path): ad_path = find_file("ad_stats.json", "..")
    
    if not model_path:
        print("Model not found.")
        for root, dirs, files in os.walk(".."):
            if "model.h5" in files:
                print(f"  Found in: {root}")
        return
        
    print(f"Found inputs: {model_path}, {scaler_path}, {ad_path}")

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

    # Try to load drug names from original input CSV (best effort)
    drug_names = {}
    try:
        import glob
        csv_files = glob.glob("../01-data-preparation/inputs/*.csv")
        if csv_files:
            import pandas as pd
            original_data = pd.read_csv(csv_files[0])
            # Find SMILES and Drug Name columns (case-insensitive)
            smiles_col = next((col for col in original_data.columns if col.lower() == 'smiles'), None)
            drug_col = next((col for col in original_data.columns if col.lower().replace(' ', '').replace('_', '') == 'drugname'), None)
            
            if smiles_col and drug_col:
                for _, row in original_data.iterrows():
                    drug_names[row[smiles_col]] = row[drug_col]
    except Exception as e:
        print(f"Could not load drug names: {e}")
    
    # Output JSON
    results = []
    for i, s in enumerate(test_smiles):
        results.append({
            "drug_name": drug_names.get(s, f"Compound {i+1}"),
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