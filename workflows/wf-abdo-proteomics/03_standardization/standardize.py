"""Node 3: Statistical Standardization - Z-score scaling of amino acid concentrations"""
import pandas as pd
from sklearn.preprocessing import StandardScaler

AMINO_ACIDS = [
    'SER_conc', 'GLN_conc', 'ARG_conc', 'CIT_conc', 'ASN_conc', '1MHIS_conc', '3MHIS_conc', 'HYP_conc', 'GLY_conc',
    'THR_conc', 'ALA_conc', 'GABA_conc', 'SAR_conc', 'BAIB_conc', 'ABA_conc', 'ORN_conc', 'MET_conc', 'PRO_conc',
    'LYS_conc', 'ASP_conc', 'HIS_conc', 'VAL_conc', 'TRP_conc', 'AAA_conc', 'LEU_conc', 'PHE_conc', 'ILE_conc',
    'C-C_conc', 'TYR_conc'
]

print("Loading cleaned data...", flush=True)
df = pd.read_pickle("../02_feature_engineering/data_cleaned.pkl")

print("Standardizing amino acid concentrations (z-score)...", flush=True)
scaler = StandardScaler()
df[AMINO_ACIDS] = scaler.fit_transform(df[AMINO_ACIDS])

df.to_pickle("data_standardized.pkl")
print("Saved: data_standardized.pkl", flush=True)
