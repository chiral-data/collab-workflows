"""Node 2: Feature Engineering - Clean and translate columns"""
import pandas as pd

print("Loading raw data...", flush=True)
df = pd.read_pickle("../01_data_ingestion/data_raw.pkl")

print("Renaming Polish columns to English...", flush=True)
df = df.rename(columns={
    'ID Pacjenta': 'ID',
    'postać': 'type',
    'Czas trwania': 'Duration',
    'wiek': 'age',
    'Plec': 'sex',
    'Lek': 'drug',
    'miejsce': 'place'
})

print("Merging disease sub-categories...", flush=True)
df['type'] = df['type'].replace({'general': 'GMG', 'eye-type': 'OMG'})

print("Adding status column...", flush=True)
df['status'] = df['type'].apply(lambda x: 'case' if x in ['PPMS', 'SPMS', 'RRMS'] else 'control')

print("Removing '_conc' suffix from column names...", flush=True)
df.columns = df.columns.str.replace('_conc', '')


df.to_pickle("data_cleaned.pkl")
print("Saved: data_cleaned.pkl", flush=True)
