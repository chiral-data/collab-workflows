"""
QSAR Workflow Node 1: Data Preparation
Responsibilities:
- Data Loading & SMILES Augmentation
- Descriptor Calculation (RDKit + MACCS)
- Variance Threshold / Feature Cleaning
- Feature Selection (Random Forest)
- Train/Test Split & Scaling
- Outlier Removal
- Saving Processed Data for Node 2
"""
import os
import sys
import pickle
import json
import base64
import io
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, MACCSkeys, Draw
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.ensemble import IsolationForest, RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectFromModel
from sklearn.decomposition import PCA
from numpy.random import seed

# Disable RDKit logging
RDLogger.DisableLog("rdApp.*")

# Reproducibility
seed(42)
np.random.seed(42)

# Ensure output directory exists (Local to this node)
OUTPUT_DIR = "outputs"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# -------------------------
# Utilities
# -------------------------
def mol_to_base64_img(mol):
    try:
        img = Draw.MolToImage(mol, size=(200, 200))
        buffer = io.BytesIO()
        img.save(buffer, format="PNG")
        return "data:image/png;base64," + base64.b64encode(buffer.getvalue()).decode("utf-8")
    except:
        return ""

def enumerate_smiles(smiles, n_variants=3):
    """Generate multiple SMILES variants for augmentation"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: return [smiles]
    
    variants = set([smiles])
    for _ in range(n_variants * 3):
        try:
            random_smiles = Chem.MolToSmiles(mol, doRandom=True)
            variants.add(random_smiles)
            if len(variants) >= n_variants + 1: break
        except: continue
    return list(variants)[:n_variants + 1]

def augment_molecular_data(smiles_list, labels, n_variants=3):
    """Augment dataset by generating SMILES variants"""
    augmented_smiles = []
    augmented_labels = []
    
    # Visualization Data: Capture one example of augmentation
    viz_example = {}
    
    print(f"Generating {n_variants} variants per molecule...")
    for i, (smi, label) in enumerate(zip(smiles_list, labels)):
        if i % 100 == 0: print(f"  Progress: {i}/{len(smiles_list)}")
        variants = enumerate_smiles(smi, n_variants)
        augmented_smiles.extend(variants)
        augmented_labels.extend([label] * len(variants))
        
        # Capture the first successful augmentation for visualization
        if i == 0 and len(variants) > 1:
            mol_orig = Chem.MolFromSmiles(variants[0])
            viz_example['original_smiles'] = variants[0]
            viz_example['original_img'] = mol_to_base64_img(mol_orig)
            viz_example['variants'] = []
            for v_smi in variants[1:]:
                v_mol = Chem.MolFromSmiles(v_smi)
                viz_example['variants'].append({
                    'smiles': v_smi,
                    'img': mol_to_base64_img(v_mol)
                })

    print(f"✓ Augmented dataset: {len(smiles_list)} → {len(augmented_smiles)} samples")
    return augmented_smiles, augmented_labels, viz_example

class RDKit_2D:
    def __init__(self, smiles):
        self.mols = [Chem.MolFromSmiles(i) for i in smiles]
        self.smiles = smiles

    def compute_2Drdkit(self, name):
        rdkit_2d_desc = []
        calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
        header = calc.GetDescriptorNames()
        for i in range(len(self.mols)):
            if self.mols[i] is not None:
                ds = calc.CalcDescriptors(self.mols[i])
                rdkit_2d_desc.append(ds)
            else:
                rdkit_2d_desc.append([0] * len(header))
        df = pd.DataFrame(rdkit_2d_desc, columns=header)
        df.insert(loc=0, column='smiles', value=self.smiles)
        return df

    def compute_MACCS(self, name):
        MACCS_list = []
        header = ['bit' + str(i) for i in range(167)]
        for i in range(len(self.mols)):
            if self.mols[i] is not None:
                ds = list(MACCSkeys.GenMACCSKeys(self.mols[i]).ToBitString())
                MACCS_list.append(ds)
            else:
                MACCS_list.append([0] * 167)
        df2 = pd.DataFrame(MACCS_list, columns=header)
        df2.insert(loc=0, column='smiles', value=self.smiles)
        return df2

# -------------------------
# Main Execution
# -------------------------
if __name__ == "__main__":
    print("="*80)
    print("NODE 1: DATA PREPARATION (SCIENTIFIC VIZ ENHANCED)")
    print("="*80)

    stats = {}

    # 1. Load Data
    input_path = "inputs/SpikeRBD_DD.csv"
    if not os.path.exists(input_path):
        print(f"❌ Error: {input_path} not found!")
        sys.exit(1)

    print(f"Step 1: Loading data from {input_path}...")
    data = pd.read_csv(input_path)
    print(f"✓ Loaded {len(data)} molecules")
    stats['count_loaded'] = len(data)

    # 2. Augmentation
    print("\nStep 2: Data Augmentation (SMILES enumeration)...")
    augmented_smiles, augmented_labels, augmentation_viz = augment_molecular_data(
        data['smiles'].tolist(),
        data['DockingScore'].tolist(),
        n_variants=3
    )
    stats['count_augmented'] = len(augmented_smiles)

    # 3. Descriptors
    print("\nStep 3: Computing descriptors...")
    RDKit_descriptor = RDKit_2D(augmented_smiles)
    x1 = RDKit_descriptor.compute_2Drdkit(None)
    x2 = RDKit_descriptor.compute_MACCS(None)
    x3 = x2.iloc[:, 1:]
    x4 = pd.concat([pd.Series(augmented_labels, name='DockingScore'), x1, x3], axis=1)

    labels = x4['DockingScore']
    features = x4.iloc[:, 3:]
    features = features.apply(pd.to_numeric, errors='coerce')
    
    # Cleaning
    features.replace([np.inf, -np.inf], np.nan, inplace=True)
    features.fillna(0, inplace=True)
    features = features.clip(lower=-1e10, upper=1e10)
    
    X_full = features.astype(np.float64)
    Y_full = np.ravel(labels).astype(np.float64)
    Y_full = np.nan_to_num(Y_full, nan=0.0)
    SMILES_full = np.array(augmented_smiles) # Track SMILES for visualization
    
    stats['count_features_raw'] = X_full.shape[1]
    stats['desc_rdkit'] = x1.shape[1] - 1 # Exclude SMILES
    stats['desc_maccs'] = x3.shape[1]

    # 4. Feature Selection
    print("\nStep 4: Feature selection using Random Forest...")
    rf = RandomForestRegressor(n_estimators=100, random_state=42, n_jobs=-1, verbose=0)
    rf.fit(X_full, Y_full)
    
    selector = SelectFromModel(rf, threshold='median', prefit=True)
    X_selected = selector.transform(X_full.values if hasattr(X_full, 'values') else X_full)
    print(f"✓ Selected features: {X_selected.shape[1]} (reduced by {100*(1-X_selected.shape[1]/X_full.shape[1]):.1f}%)")
    
    stats['count_features_selected'] = X_selected.shape[1]

    # Capture Feature Importance for Viz
    importances = rf.feature_importances_
    indices = np.argsort(importances)[::-1]
    
    # Cumulative Importance Curve
    cumulative_importance = np.cumsum(importances[indices])
    
    top_n = 50 # Increased for Dashboard Slider
    top_features = []
    if hasattr(features, 'columns'):
        all_feat_names = features.columns
        for i in range(min(top_n, len(indices))):
            top_features.append({
                "name": all_feat_names[indices[i]],
                "importance": float(importances[indices[i]])
            })

    # 5. Split & Scale
    print("\nStep 5: Splitting & Scaling...")
    X_train, X_test, y_train, y_test, smiles_train, smiles_test = train_test_split(
        X_selected, Y_full, SMILES_full, test_size=0.30, shuffle=True, random_state=42
    )
    
    # Save raw split stats
    stats['count_train_initial'] = len(X_train)
    stats['count_test_initial'] = len(X_test)
    
    # Store labels for Histogram Overlay
    viz_train_labels = y_train.tolist()
    viz_test_labels = y_test.tolist()

    sc = StandardScaler()
    X_train = sc.fit_transform(X_train)
    X_test = sc.transform(X_test)
    X_train = np.nan_to_num(X_train, nan=0.0)
    X_test = np.nan_to_num(X_test, nan=0.0)

    # 6. Outlier Removal (with detailed Viz capture)
    print("\nStep 6: Removing outliers (Isolation Forest)...")
    iso = IsolationForest(contamination=0.1, n_estimators=100, random_state=42, verbose=0)
    
    # Compute PCA first for ALL training data (to see outliers in context)
    pca = PCA(n_components=2)
    X_train_pca = pca.fit_transform(X_train)
    pca_variance = pca.explained_variance_ratio_.tolist()
    
    # Train outliers
    yhat_train = iso.fit_predict(X_train_pca) # Predict using PCA components or full space? Usually full, but staying consistent with original approach which used PCA for fitting iso? 
    # WAIT: Original code used pca.fit_transform in the loop, implying it learned on full space but transformed to 2D for IsoForest??
    # "yhat_train = iso.fit_predict(pca.fit_transform(X_train))" -> Yes, it used 2D PCA for outlier detection.
    # We will stick to that logic to avoid changing model behavior.
    
    mask_train = yhat_train != -1
    
    # Separate Inliers and Outliers for Viz
    inliers_pca = X_train_pca[mask_train]
    outliers_pca = X_train_pca[~mask_train]
    inliers_y = y_train[mask_train]
    outliers_y = y_train[~mask_train]
    
    # Apply filter
    X_train, y_train, smiles_train = X_train[mask_train], y_train[mask_train], smiles_train[mask_train]
    
    # Test outliers (Repeat same logic)
    X_test_pca = pca.transform(X_test) # Use same PCA
    yhat_test = iso.fit_predict(X_test_pca)
    mask_test = yhat_test != -1
    X_test, y_test, smiles_test = X_test[mask_test], y_test[mask_test], smiles_test[mask_test]
    
    print(f"✓ Final Training set: {X_train.shape}")
    print(f"✓ Final Test set: {X_test.shape}")
    
    stats['count_train_final'] = len(X_train)
    stats['count_test_final'] = len(X_test)
    stats['outliers_removed_train'] = int(np.sum(~mask_train))

    # 7. Save Output (Processed Data)
    print(f"\nStep 7: Saving processed data for Node 2 to {OUTPUT_DIR}...")
    output_file = os.path.join(OUTPUT_DIR, "processed_data.pkl")
    
    data_bundle = {
        'X_train': X_train,
        'y_train': y_train,
        'X_test': X_test,
        'y_test': y_test,
        'smiles_train': smiles_train, # Added for Node 2 Visualization
        'smiles_test': smiles_test,   # Added for Node 2 Visualization
        'scaler': sc,
        'selector': selector,
        'feature_names': list(features.columns[selector.get_support()]) if hasattr(features, 'columns') else []
    }
    
    with open(output_file, 'wb') as f:
        pickle.dump(data_bundle, f)
    print(f"✓ Saved processed data to: {output_file}")

    # 8. Save JSON Data for Visualization
    print(f"\nStep 8: Generatng JSON Report Data...")
    json_output_file = os.path.join(OUTPUT_DIR, "data.json")

    # Correlation Matrix (Top 20 features)
    top_indices = indices[:20]
    X_top = X_full.iloc[:, top_indices] if hasattr(X_full, 'iloc') else pd.DataFrame(X_full).iloc[:, top_indices]
    corr_matrix = X_top.corr()

    json_data = {
        "summary": stats,
        "augmentation": augmentation_viz,
        "sankey_counts": {
            "loaded": stats['count_loaded'],
            "augmented": stats['count_augmented'],
            "selected": stats['count_augmented'], # same samples, fewer features
            "split_train": stats['count_train_initial'],
            "split_test": stats['count_test_initial'],
            "final_train": stats['count_train_final'],
            "final_test": stats['count_test_final']
        },
        "feature_selection": {
            "top_features": top_features,
            "cumulative_importance": cumulative_importance.tolist()[:100], # Top 100 for curve
            "n_features_total": int(X_full.shape[1])
        },
        "pca_analysis": {
            "variance_ratio": pca_variance,
            "inliers": {
                "x": inliers_pca[:, 0].tolist(),
                "y": inliers_pca[:, 1].tolist(),
                "c": inliers_y.tolist(),
                "s": smiles_train.tolist() # Add SMILES to PCA data
            },
            "outliers": {
                "x": outliers_pca[:, 0].tolist(),
                "y": outliers_pca[:, 1].tolist(),
                "c": outliers_y.tolist()
            }
        },
        "distributions": {
            "train_labels": viz_train_labels,
            "test_labels": viz_test_labels
        },
        "correlation": {
            "z": corr_matrix.values.tolist(),
            "x": list(corr_matrix.columns),
            "y": list(corr_matrix.index)
        }
    }
    
    with open(json_output_file, 'w') as f:
        json.dump(json_data, f, indent=4)
    print(f"✓ Saved JSON report data to: {json_output_file}")

    print("="*80)
    print("NODE 1 COMPLETION SUCCESSFUL")
    print("="*80)
