import pandas as pd
import numpy as np
import os
import glob
import json
import base64
from io import BytesIO
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, MACCSkeys, Draw
from rdkit.ML.Descriptors import MoleculeDescriptors

# Disable RDKit logging
RDLogger.DisableLog("rdApp.*")

class RDKit_2D:
    def __init__(self, smiles):
        self.smiles = smiles if isinstance(smiles, list) else smiles.tolist()
        self.mols = [Chem.MolFromSmiles(s) for s in self.smiles]

    def compute_2Drdkit(self):
        rdkit_2d_desc = []
        calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
        header = calc.GetDescriptorNames()
        for mol in self.mols:
            if mol:
                ds = calc.CalcDescriptors(mol)
            else:
                ds = [np.nan] * len(header)
            rdkit_2d_desc.append(ds)
        df = pd.DataFrame(rdkit_2d_desc, columns=header)
        df.insert(loc=0, column='smiles', value=self.smiles)
        return df

    def compute_MACCS(self):
        MACCS_list = []
        header = [f'bit{i}' for i in range(167)]
        for mol in self.mols:
            if mol:
                ds = list(MACCSkeys.GenMACCSKeys(mol).ToBitString())
            else:
                ds = [0] * 167
            MACCS_list.append(ds)
        df = pd.DataFrame(MACCS_list, columns=header)
        df.insert(loc=0, column='smiles', value=self.smiles)
        return df

def generate_images(df, n=24):
    """Generate base64 images for top N molecules for visualization"""
    images = []
    for _, row in df.head(n).iterrows():
        smi = row.get('smiles')
        if not smi: continue
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                img = Draw.MolToImage(mol, size=(200, 200))
                buffer = BytesIO()
                img.save(buffer, format="PNG")
                img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")
                images.append({"smiles": smi, "img": f"data:image/png;base64,{img_str}"})
        except:
            pass
    return images

def run():
    print("--- Node 1: Data Preparation ---")
    os.makedirs("outputs", exist_ok=True)

    # Find CSV file in inputs/ directory (copied from input_files/ by silva)
    csv_files = glob.glob("inputs/*.csv")
    if not csv_files:
        print("Error: No CSV file found in inputs/ directory.")
        return

    input_file = csv_files[0]  # Use the first CSV file found
    print(f"Using input file: {input_file}")

    data = pd.read_csv(input_file)
    print(f"Loaded {len(data)} rows.")
    
    # Compute Descriptors
    print("Computing descriptors...")
    desc_eng = RDKit_2D(data['smiles'])
    x1 = desc_eng.compute_2Drdkit()
    x2 = desc_eng.compute_MACCS()
    
    # Merge
    x3 = x2.iloc[:, 1:] # Drop duplicated smiles col
    # Keep DockingScore, smiles, then descriptors
    final_df = pd.concat([data['DockingScore'], x1, x3], axis=1)
    
    # Save CSV
    final_df.to_csv("outputs/descriptors.csv", index=False)
    print("Saved outputs/descriptors.csv")
    
    # Prepare JSON Data for Report
    # Calculate correlation matrix for a subset of features to avoid huge JSON
    numeric_df = x1.select_dtypes(include=[np.number]).iloc[:, :20] 
    corr_matrix = numeric_df.corr().round(2)
    
    json_data = {
        "total_samples": len(final_df),
        "feature_count": final_df.shape[1] - 2,
        "correlation": {
            "z": corr_matrix.values.tolist(),
            "x": corr_matrix.columns.tolist(),
            "y": corr_matrix.index.tolist()
        },
        "sample_images": generate_images(data)
    }
    
    with open("outputs/data.json", "w") as f:
        json.dump(json_data, f)
    print("Saved outputs/data.json")

if __name__ == "__main__":
    run()