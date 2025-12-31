"""Node 1: Data Ingestion - Load TSV proteomics database"""
import os
import pandas as pd

data_file = os.environ.get("PARAM_DATA_FILE", "database-multiple-sclerosis-myasthenia.csv")

print(f"Loading data from: {data_file}", flush=True)
df = pd.read_csv(data_file, sep='\t', header=0)
print(f"Loaded {len(df)} records with {len(df.columns)} columns", flush=True)

df.to_pickle("data_raw.pkl")
print("Saved: data_raw.pkl", flush=True)
