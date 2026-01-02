"""Node 8: Global Metabolic Load Analysis - Generate Fig 8"""
import os
import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from html_generator import generate_metabolic_load_html

# Configuration
INPUT_FILE = "../01_Data_Ingestion_and_Preprocessing/output/data_standardized.pkl"
OUTPUT_DIR = "output"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'sans-serif'

print(">>> NODE 8: GLOBAL METABOLIC LOAD ANALYSIS...", flush=True)

# Load data
df = pd.read_pickle(INPUT_FILE)

# Define masks
ms_types = ['RRMS', 'SPMS', 'PPMS']
mg_types = ['general', 'eye-type']
masks = {
    'MS': df['Type'].isin(ms_types),
    'MG': df['Type'].isin(mg_types)
}

df_mg_ms = df[masks['MS'] | masks['MG']].copy()

# ==========================================
# FIG 8: Total AA by Type
# ==========================================
print("Generating Fig 8: Total AA by Type...", flush=True)

df_mg_ms['Plot_Type'] = np.where(df_mg_ms['Type'].isin(['general', 'eye-type']), 'GMG', df_mg_ms['Type'])

plt.figure(figsize=(8,6))

type_order = ['PPMS', 'SPMS', 'RRMS', 'GMG']
type_colors = {'PPMS': '#1f77b4', 'SPMS': '#ff7f0e', 'RRMS': '#2ca02c', 'GMG': '#d62728'}

sns.boxplot(data=df_mg_ms, x='Plot_Type', y='Total_AA', order=type_order,
            palette=type_colors, showfliers=True, width=0.5, linewidth=1.5,
            flierprops={'marker': 'o', 'markerfacecolor': 'white', 'markeredgecolor': 'black', 'markersize': 5})

plt.ylabel('Concentration [nmol/ml]')
plt.xlabel('')
plt.yticks([-4, -2, 0, 2, 4, 6, 8, 10, 12])
plt.title('Total Amino Acid Concentration by Disease Type', fontsize=12, fontweight='bold')

plt.savefig(os.path.join(OUTPUT_DIR, 'Fig8_TotalAA_Type.png'))
plt.close()
print(f"Saved: Fig8_TotalAA_Type.png", flush=True)


# ==========================================
# GENERATE JSON/HTML DATA
# ==========================================
print("Generating JSON and HTML...", flush=True)

# Create JSON data structure (simplified - enhance as needed)
json_data = {'metadata': {'title': 'Global Metabolic Load'}}

json_path = os.path.join(OUTPUT_DIR, 'metabolic_load_data.json')
with open(json_path, 'w', encoding='utf-8') as f:
    json.dump(json_data, f, indent=2)
print(f"Saved: {json_path}", flush=True)

html_content = generate_metabolic_load_html(json_filename='metabolic_load_data.json')

html_path = os.path.join(OUTPUT_DIR, 'metabolic_load.html')
with open(html_path, 'w', encoding='utf-8') as f:
    f.write(html_content)
print(f"Saved: {html_path}", flush=True)

print("Node 8 completed successfully.", flush=True)
