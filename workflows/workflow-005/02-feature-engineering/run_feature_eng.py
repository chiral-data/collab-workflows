import pandas as pd
import numpy as np
import os
import json
import pickle
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.ensemble import IsolationForest

def run():
    print("--- Node 2: Feature Engineering ---")
    os.makedirs("outputs", exist_ok=True)

    input_path = "../01_Data_Preparation/outputs/descriptors.csv"
    if not os.path.exists(input_path):
        print(f"Error: {input_path} not found.")
        return

    df = pd.read_csv(input_path)
    
    # 1. Cleaning
    # Assuming col 0 is DockingScore, col 1 is smiles, rest are features
    labels = df['DockingScore']
    features = df.iloc[:, 2:] # Skip DockingScore and smiles
    
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
    X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.3, random_state=0)
    
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
    
    iso = IsolationForest(contamination=0.1, random_state=3)
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