"""Node 6: Data Transformation - Wide to Long format for visualization"""
import pandas as pd

AMINO_ACIDS = [
    'SER_conc', 'GLN_conc', 'ARG_conc', 'CIT_conc', 'ASN_conc', '1MHIS_conc', '3MHIS_conc', 'HYP_conc', 'GLY_conc',
    'THR_conc', 'ALA_conc', 'GABA_conc', 'SAR_conc', 'BAIB_conc', 'ABA_conc', 'ORN_conc', 'MET_conc', 'PRO_conc',
    'LYS_conc', 'ASP_conc', 'HIS_conc', 'VAL_conc', 'TRP_conc', 'AAA_conc', 'LEU_conc', 'PHE_conc', 'ILE_conc',
    'C-C_conc', 'TYR_conc'
]

ID_VARS = ['status', 'ID', 'type', 'Duration', 'EDSS', 'age', 'sex', 'drug', 'place']

print("Loading standardized data...", flush=True)
df = pd.read_pickle("../03_standardization/data_standardized.pkl")

print("Transforming from wide to long format...", flush=True)
melted = df.melt(
    id_vars=ID_VARS,
    value_vars=AMINO_ACIDS,
    var_name='AA',
    value_name='conc'
)

print(f"Melted data shape: {melted.shape}", flush=True)

melted.to_pickle("data_melted.pkl")
print("Saved: data_melted.pkl", flush=True)
