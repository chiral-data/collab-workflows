import pandas as pd
import numpy as np
import os
import json
import pickle
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.ensemble import IsolationForest

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
            section_match = re.search(f'\\[params\\.{param_name}\\]', content, re.IGNORECASE)
            if section_match:
                start_index = section_match.end()
                default_match = re.search(r'default\s*=\s*(.+)', content[start_index:])
                if default_match:
                    raw_val = default_match.group(1).strip()
                    if raw_val.startswith('"') and raw_val.endswith('"'):
                         return raw_val[1:-1].replace('\\"', '"')
                    elif raw_val.startswith("'") and raw_val.endswith("'"):
                         return raw_val[1:-1]
                    else:
                         return raw_val
    except Exception as e:
        print(f"Warning: Failed to parse job.toml for {param_name}: {e}")
        
    if default_if_missing is not None:
        return default_if_missing
        
    raise ValueError(f"Parameter {param_name} missing (Env Var and job.toml default)")

def run():
    print("--- Node 2: Feature Engineering ---")
    os.makedirs("outputs", exist_ok=True)

    # Silva shared workspace: read directly from previous node's outputs
    input_path = os.path.join("..", "01-data-preparation", "outputs", "descriptors.csv")
    
    if not os.path.exists(input_path):
        # Fallback to find in current directory or inputs
        print(f"Path {input_path} not found, searching for descriptors.csv...")
        input_path = find_file("descriptors.csv", "..")

    if not input_path:
        print(f"Error: descriptors.csv not found.")
        # Debug listing
        print("DEBUG: Listing parent directory:")
        for root, dirs, files in os.walk(".."):
            if "descriptors.csv" in files:
                print(f"  Found in: {root}")
        return

    print(f"Found input: {input_path}")
    df = pd.read_csv(input_path)
    
    # 1. Cleaning
    # Column structure: Drug_Name (if present), DockingScore, smiles, then numeric descriptors
    # Extract labels first
    score_col = next((col for col in df.columns if col.lower() == 'dockingscore'), None)
    if not score_col:
        print("Error: No DockingScore column found")
        return
    
    labels = df[score_col]
    
    # Get only numeric columns (skip Drug_Name, DockingScore, smiles)
    # This automatically handles whether Drug_Name exists or not
    features = df.select_dtypes(include=[np.number])
    
    # Drop DockingScore if it's in the numeric columns
    if score_col in features.columns:
        features = features.drop(columns=[score_col])
    
    # Force numeric and clean
    features = features.apply(pd.to_numeric, errors='coerce')
    features.replace([np.inf, -np.inf], np.nan, inplace=True)
    features.dropna(inplace=True)
    features = features.clip(lower=-1e10, upper=1e10)
    
    # Align labels
    labels = labels.loc[features.index]
    
    X = features.values.astype(np.float64)
    Y = labels.values.astype(np.float64)
    
    # Save stats for Applicability Domain (AD)
    ad_stats = {
        "min": features.min().tolist(),
        "max": features.max().tolist(),
        "columns": features.columns.tolist()
    }
    with open("outputs/ad_stats.json", "w") as f:
        json.dump(ad_stats, f)

    # 2. Split & Scale
    test_size = float(get_param("test_size"))
    X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=test_size, random_state=0)
    
    scaler = StandardScaler()
    X_train_s = scaler.fit_transform(X_train)
    X_test_s = scaler.transform(X_test)
    
    # Save Scaler
    with open("outputs/scaler.pkl", "wb") as f:
        pickle.dump(scaler, f)

    # 3. Outlier Removal (PCA + IsolationForest)
    # CORRECT LOGIC: Fit on Train, Transform Test
    pca = PCA(n_components=2)
    X_train_pca = pca.fit_transform(X_train_s)
    X_test_pca = pca.transform(X_test_s)
    
    contamination = float(get_param("contamination"))
    iso = IsolationForest(contamination=contamination, random_state=3)
    # Fit and predict on train
    train_outliers = iso.fit_predict(X_train_pca) 
    # Predict on test
    test_outliers = iso.predict(X_test_pca)
    
    # Filter Train Data
    mask_train = train_outliers != -1
    X_train_final = X_train_s[mask_train]
    y_train_final = y_train[mask_train]
    
    # Filter Test Data (Optional, but done in original script)
    mask_test = test_outliers != -1
    X_test_final = X_test_s[mask_test]
    y_test_final = y_test[mask_test]
    
    # Save Processed Data
    np.savez("outputs/processed_data.npz", 
             X_train=X_train_final, y_train=y_train_final, 
             X_test=X_test_final, y_test=y_test_final)
    print("Saved processed_data.npz")

    # JSON for Visualization
    json_data = {
        "train_samples": len(X_train_final),
        "test_samples": len(X_test_final),
        "removed_train": int(np.sum(train_outliers == -1)),
        "pca_data": {
            "train_x": X_train_pca[:, 0].tolist(),
            "train_y": X_train_pca[:, 1].tolist(),
            "train_label": ["Inlier" if x == 1 else "Outlier" for x in train_outliers],
            "test_x": X_test_pca[:, 0].tolist(),
            "test_y": X_test_pca[:, 1].tolist(),
        }
    }
    with open("outputs/data.json", "w") as f:
        json.dump(json_data, f)
    print("Saved outputs/data.json")

if __name__ == "__main__":
    run()