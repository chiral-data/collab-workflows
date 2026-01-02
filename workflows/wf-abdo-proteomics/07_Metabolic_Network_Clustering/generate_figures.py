"""Node 7: Metabolic Network Clustering - Generate Fig 6 & 7"""
import os
import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from html_generator import generate_clustering_html

# Configuration
INPUT_FILE = "../01_Data_Ingestion_and_Preprocessing/output/data_standardized.pkl"
AA_COLS_FILE = "../01_Data_Ingestion_and_Preprocessing/output/aa_cols.txt"
OUTPUT_DIR = "output"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'sans-serif'

print(">>> NODE 7: METABOLIC NETWORK CLUSTERING...", flush=True)

# Load data
df = pd.read_pickle(INPUT_FILE)
with open(AA_COLS_FILE, 'r') as f:
    aa_cols = [line.strip() for line in f.readlines()]

# Define masks
ms_types = ['RRMS', 'SPMS', 'PPMS']
mg_types = ['general', 'eye-type']
masks = {
    'MS': df['Type'].isin(ms_types),
    'MG': df['Type'].isin(mg_types),
    'Control': df['Status'] == 'control'
}

df_mg_ms = df[masks['MS'] | masks['MG']].copy()

# Clean labels (remove '_conc' suffix)
aa_labels_clean = [col.replace('_conc', '') for col in aa_cols]

# ==========================================
# FIG 6: MS+MG Clustered Correlation
# ==========================================
print("Generating Fig 6: MS+MG Correlation Matrix...", flush=True)

corr_ms_mg = df_mg_ms[aa_cols].corr()
corr_ms_mg.index = aa_labels_clean
corr_ms_mg.columns = aa_labels_clean

g6 = sns.clustermap(corr_ms_mg, cmap='YlGnBu', figsize=(12, 12),
                    linewidths=0.5, linecolor='white',
                    dendrogram_ratio=(0.15, 0.15),
                    cbar_pos=(0.02, 0.83, 0.02, 0.15),
                    cbar_kws={'shrink': 1.0},
                    vmin=-0.2, vmax=1.0)
g6.ax_heatmap.set_title('')
g6.ax_heatmap.tick_params(axis='both', which='major', labelsize=8)
g6.figure.suptitle('MS+MG Correlation Matrix (Clustered)', fontsize=16, fontweight='bold', y=1.02)
g6.savefig(os.path.join(OUTPUT_DIR, 'Fig6_Corr_MS_MG.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: Fig6_Corr_MS_MG.png", flush=True)

# ==========================================
# FIG 7: Control Clustered Correlation
# ==========================================
print("Generating Fig 7: Control Correlation Matrix...", flush=True)

corr_ctrl = df[masks['Control']][aa_cols].corr()
corr_ctrl.index = aa_labels_clean
corr_ctrl.columns = aa_labels_clean

g7 = sns.clustermap(corr_ctrl, cmap='YlGnBu', figsize=(12, 12),
                    linewidths=0.5, linecolor='white',
                    dendrogram_ratio=(0.15, 0.15),
                    cbar_pos=(0.02, 0.83, 0.02, 0.15),
                    cbar_kws={'shrink': 1.0},
                    vmin=-0.2, vmax=1.0)
g7.ax_heatmap.set_title('')
g7.ax_heatmap.tick_params(axis='both', which='major', labelsize=8)
g7.figure.suptitle('Control Correlation Matrix (Clustered)', fontsize=16, fontweight='bold', y=1.02)
g7.savefig(os.path.join(OUTPUT_DIR, 'Fig7_Corr_Control.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: Fig7_Corr_Control.png", flush=True)


# ==========================================
# GENERATE JSON DATA
# ==========================================
print("Generating JSON data...", flush=True)

json_data = {
    'metadata': {'title': 'Metabolic Network Clustering'},
    'fig6': {
        'x': corr_ms_mg.columns.str.replace('_conc', '').tolist(),
        'y': corr_ms_mg.index.str.replace('_conc', '').tolist(),
        'z': corr_ms_mg.values.tolist(),
        'title': 'MS & MG Correlation Matrix',
        'colorscale': 'RdBu'
    },
    'fig7': {
        'x': corr_ctrl.columns.str.replace('_conc', '').tolist(),
        'y': corr_ctrl.index.str.replace('_conc', '').tolist(),
        'z': corr_ctrl.values.tolist(),
        'title': 'Control Correlation Matrix',
        'colorscale': 'RdBu'
    }
}

json_path = os.path.join(OUTPUT_DIR, 'clustering_data.json')
with open(json_path, 'w', encoding='utf-8') as f:
    json.dump(json_data, f, indent=2)
print(f"Saved: {json_path}", flush=True)

html_content = generate_clustering_html(json_filename='clustering_data.json')

html_path = os.path.join(OUTPUT_DIR, 'clustering.html')
with open(html_path, 'w', encoding='utf-8') as f:
    f.write(html_content)
print(f"Saved: {html_path}", flush=True)

print("Node 7 completed successfully.", flush=True)
