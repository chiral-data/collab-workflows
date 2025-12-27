"""Node 3: Statistical Standardization - Z-score scaling of amino acid concentrations"""
import os
import pandas as pd
from pathlib import Path
from sklearn.preprocessing import StandardScaler

# Setup paths
current_dir = Path(__file__).parent
input_file = current_dir / "data_cleaned.pkl"
output_file = current_dir / "data_standardized.pkl"

# Load parameters from environment variables
DEFAULT_AMINO_ACIDS = "SER,GLN,ARG,CIT,ASN,1MHIS,3MHIS,HYP,GLY,THR,ALA,GABA,SAR,BAIB,ABA,ORN,MET,PRO,LYS,ASP,HIS,VAL,TRP,AAA,LEU,PHE,ILE,C-C,TYR"
AMINO_ACIDS = os.environ.get("PARAM_AMINO_ACIDS", DEFAULT_AMINO_ACIDS).split(",")
print(f"Using {len(AMINO_ACIDS)} amino acids from environment", flush=True)

print(f"Loading cleaned data from {input_file}...", flush=True)
df = pd.read_pickle(input_file)

print("Standardizing amino acid concentrations (z-score)...", flush=True)
scaler = StandardScaler()
df[AMINO_ACIDS] = scaler.fit_transform(df[AMINO_ACIDS])

df.to_pickle(output_file)
print(f"Saved: {output_file.name}", flush=True)
