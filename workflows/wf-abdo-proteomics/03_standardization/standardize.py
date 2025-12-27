"""Node 3: Statistical Standardization - Z-score scaling of amino acid concentrations"""
import json
import pandas as pd
from pathlib import Path
from sklearn.preprocessing import StandardScaler

# Setup paths
current_dir = Path(__file__).parent
root_dir = current_dir.parent
input_file = current_dir / "data_cleaned.pkl"
metadata_file = root_dir / "global_params.json"
output_file = current_dir / "data_standardized.pkl"

# Load metadata
print(f"Loading metadata from {metadata_file}...", flush=True)
with open(metadata_file, 'r') as f:
    metadata = json.load(f)
    if 'amino_acids_Conc' in metadata:
        AMINO_ACIDS = metadata['amino_acids_Conc']
    else:
        AMINO_ACIDS = metadata.get('amino_acids', [])

print(f"Loading cleaned data from {input_file}...", flush=True)
df = pd.read_pickle(input_file)

print("Standardizing amino acid concentrations (z-score)...", flush=True)
scaler = StandardScaler()
df[AMINO_ACIDS] = scaler.fit_transform(df[AMINO_ACIDS])

df.to_pickle(output_file)
print(f"Saved: {output_file.name}", flush=True)
