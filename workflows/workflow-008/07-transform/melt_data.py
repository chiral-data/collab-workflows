"""Node 7: Data Transformation - Wide to Long format for visualization"""
import os
import pandas as pd

# Setup paths
input_file = "data_standardized.pkl"
output_file = "data_melted.pkl"

ID_VARS = ['status', 'ID', 'type', 'Duration', 'EDSS', 'age', 'sex', 'drug', 'place']

# Load parameters from environment variables
DEFAULT_AMINO_ACIDS = "SER,GLN,ARG,CIT,ASN,1MHIS,3MHIS,HYP,GLY,THR,ALA,GABA,SAR,BAIB,ABA,ORN,MET,PRO,LYS,ASP,HIS,VAL,TRP,AAA,LEU,PHE,ILE,C-C,TYR"
AMINO_ACIDS = os.environ.get("PARAM_AMINO_ACIDS", DEFAULT_AMINO_ACIDS).split(",")
print(f"Using {len(AMINO_ACIDS)} amino acids from environment", flush=True)

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
print(f"Saved: {output_file}", flush=True)
