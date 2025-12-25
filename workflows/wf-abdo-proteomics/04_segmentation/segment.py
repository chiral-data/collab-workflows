"""Node 4: Group Segmentation - Split data into Case and Control groups"""
import pandas as pd

print("Loading standardized data...", flush=True)
df = pd.read_pickle("../03_standardization/data_standardized.pkl")

print("Segmenting into Case and Control groups...", flush=True)
case = df[df['status'] == 'case'].copy()
control = df[df['status'] == 'control'].copy()

print(f"Cases: {len(case)}, Controls: {len(control)}", flush=True)

case.to_pickle("case.pkl")
control.to_pickle("control.pkl")
print("Saved: case.pkl, control.pkl", flush=True)
