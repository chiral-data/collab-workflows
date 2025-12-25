"""Node 6: Data Transformation - Wide to Long format for visualization"""
import json
import pandas as pd
from pathlib import Path

# Setup paths
current_dir = Path(__file__).parent
root_dir = current_dir.parent
input_file = root_dir / "03_standardization" / "data_standardized.pkl"
metadata_file = root_dir / "global_params.json"
output_file = current_dir / "data_melted.pkl"

ID_VARS = ['status', 'ID', 'type', 'Duration', 'EDSS', 'age', 'sex', 'drug', 'place']

# Load metadata
print(f"Loading metadata from {metadata_file}...", flush=True)
with open(metadata_file, 'r') as f:
    metadata = json.load(f)
    if 'amino_acids_Conc' in metadata:
        AMINO_ACIDS = metadata['amino_acids_Conc']
    else:
        AMINO_ACIDS = metadata.get('amino_acids', [])

print(f"Loading standardized data from {input_file}...", flush=True)
df = pd.read_pickle(input_file)

print("Transforming from wide to long format...", flush=True)
melted = df.melt(
    id_vars=ID_VARS,
    value_vars=AMINO_ACIDS,
    var_name='AA',
    value_name='conc'
)

print(f"Melted data shape: {melted.shape}", flush=True)

melted.to_pickle(output_file)
print(f"Saved: {output_file.name}", flush=True)
