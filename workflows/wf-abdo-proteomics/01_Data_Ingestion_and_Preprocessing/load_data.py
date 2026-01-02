"""Node 1: Data Ingestion & Preprocessing"""
import os
import pandas as pd
import numpy as np

# Configuration
DATA_FILE = os.environ.get("PARAM_DATA_FILE", "../01_Data_Ingestion_and_Preprocessing/database-multiple-sclerosis-myasthenia.csv")
OUTPUT_DIR = "output"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

print(">>> NODE 1: DATA INGESTION & PREPROCESSING...", flush=True)

# Load data
print(f"Loading data from: {DATA_FILE}", flush=True)
try:
    df = pd.read_csv(DATA_FILE, sep='\t')
except:
    df = pd.read_csv(DATA_FILE)

print(f"Loaded {len(df)} records with {len(df.columns)} columns", flush=True)

# Rename Columns
df = df.rename(columns={
    'postać': 'Type', 'status': 'Status', 'wiek': 'Age', 
    'Plec': 'Sex', 'Czas trwania': 'Duration'
})

# Amino Acids
aa_cols = ['SER_conc', 'GLN_conc', 'ARG_conc', 'CIT_conc', 'ASN_conc', 
           '1MHIS_conc', '3MHIS_conc', 'HYP_conc', 'GLY_conc', 'THR_conc', 
           'ALA_conc', 'GABA_conc', 'SAR_conc', 'BAIB_conc', 'ABA_conc', 
           'ORN_conc', 'MET_conc', 'PRO_conc', 'LYS_conc', 'ASP_conc', 
           'HIS_conc', 'VAL_conc', 'TRP_conc', 'AAA_conc', 'LEU_conc', 
           'PHE_conc', 'ILE_conc', 'C-C_conc', 'TYR_conc']

# Standardization
df_std = df.copy()
# Calculate Total AA from raw concentrations first
df_std['Total_AA'] = df[aa_cols].sum(axis=1)
# Then standardize the individual amino acid columns and Total_AA
cols_to_standardize = aa_cols + ['Total_AA']
df_std[cols_to_standardize] = (df_std[cols_to_standardize] - df_std[cols_to_standardize].mean()) / df_std[cols_to_standardize].std()

# Clean DMTs
dmt_map = {
    'betaferon': 'Interferon-beta', 'Copaxone': 'Glatiramer acetate',
    'Tecfidera': 'Dimethyl fumarate', 'mitoksantron': 'Mitoxantrone',
    'Tysabri': 'Natalizumab', 'natalizumab': 'Natalizumab',
    'fingolimod': 'Fingolimod'
}
df_std['DMT_Clean'] = df['Lek'].map(dmt_map)

# Define Cohorts
ms_types = ['RRMS', 'SPMS', 'PPMS']
mg_types = ['general', 'eye-type']

# Save processed data
df_std.to_pickle(os.path.join(OUTPUT_DIR, "data_standardized.pkl"))
print(f"Saved: {os.path.join(OUTPUT_DIR, 'data_standardized.pkl')}", flush=True)

# Save amino acid columns list
with open(os.path.join(OUTPUT_DIR, "aa_cols.txt"), 'w') as f:
    f.write('\n'.join(aa_cols))
print(f"Saved: {os.path.join(OUTPUT_DIR, 'aa_cols.txt')}", flush=True)

print(f"Processing complete. Standardized {len(df_std)} records.", flush=True)
