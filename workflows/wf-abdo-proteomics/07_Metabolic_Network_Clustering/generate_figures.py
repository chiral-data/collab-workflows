"""Node 7: Metabolic Network Clustering - Generate Fig 6 & 7"""
import os
import sys
import json
import base64
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram

from html_generator import generate_clustering_html

# Configuration
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_FILE = os.path.join(SCRIPT_DIR, "../01_Data_Ingestion_and_Preprocessing/output/data_standardized.pkl")
AA_COLS_FILE = os.path.join(SCRIPT_DIR, "../01_Data_Ingestion_and_Preprocessing/output/aa_cols.txt")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "output")

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

# Capture clustered data for JSON
dendro_ms_mg = {'row': None, 'col': None}
try:
    if hasattr(g6, 'data2d'):
        corr_ms_mg_clustered = g6.data2d
    else:
        # Reconstruct from indices
        row_idx = g6.dendrogram_row.reordered_ind
        col_idx = g6.dendrogram_col.reordered_ind
        corr_ms_mg_clustered = corr_ms_mg.iloc[row_idx, col_idx]
    
    # Capture dendrogram traces
    if hasattr(g6, 'dendrogram_row') and g6.dendrogram_row is not None:
        d = dendrogram(g6.dendrogram_row.linkage, no_plot=True)
        # Scale icoords to match heatmap indices (0..N)
        # Scipy dendrogram uses 10x scaling and centers at 5, 15, 25...
        icoords = [[(x - 5.0) / 10.0 for x in seg] for seg in d['icoord']]
        dendro_ms_mg['row'] = {'icoords': icoords, 'dcoords': d['dcoord']}
    if hasattr(g6, 'dendrogram_col') and g6.dendrogram_col is not None:
        d = dendrogram(g6.dendrogram_col.linkage, no_plot=True)
        # Scale icoords
        icoords = [[(x - 5.0) / 10.0 for x in seg] for seg in d['icoord']]
        dendro_ms_mg['col'] = {'icoords': icoords, 'dcoords': d['dcoord']}
        
except Exception as e:
    print(f"Warning: Could not capture clustering order for Fig 6: {e}")
    corr_ms_mg_clustered = corr_ms_mg

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

# Capture clustered data for JSON
dendro_ctrl = {'row': None, 'col': None}
try:
    if hasattr(g7, 'data2d'):
        corr_ctrl_clustered = g7.data2d
    else:
        # Reconstruct from indices
        row_idx = g7.dendrogram_row.reordered_ind
        col_idx = g7.dendrogram_col.reordered_ind
        corr_ctrl_clustered = corr_ctrl.iloc[row_idx, col_idx]

    # Capture dendrogram traces
    if hasattr(g7, 'dendrogram_row') and g7.dendrogram_row is not None:
        d = dendrogram(g7.dendrogram_row.linkage, no_plot=True)
        # Scale icoords
        icoords = [[(x - 5.0) / 10.0 for x in seg] for seg in d['icoord']]
        dendro_ctrl['row'] = {'icoords': icoords, 'dcoords': d['dcoord']}
    if hasattr(g7, 'dendrogram_col') and g7.dendrogram_col is not None:
        d = dendrogram(g7.dendrogram_col.linkage, no_plot=True)
        # Scale icoords
        icoords = [[(x - 5.0) / 10.0 for x in seg] for seg in d['icoord']]
        dendro_ctrl['col'] = {'icoords': icoords, 'dcoords': d['dcoord']}

except Exception as e:
    print(f"Warning: Could not capture clustering order for Fig 7: {e}")
    corr_ctrl_clustered = corr_ctrl

plt.close()
print(f"Saved: Fig7_Corr_Control.png", flush=True)


# ==========================================
# GENERATE JSON DATA
# ==========================================
print("Generating JSON data...", flush=True)

# Helper to encode images
def get_base64_image(filename):
    try:
        with open(os.path.join(OUTPUT_DIR, filename), "rb") as img_file:
            return f"data:image/png;base64,{base64.b64encode(img_file.read()).decode('utf-8')}"
    except Exception as e:
        print(f"Warning: Could not encode image {filename}: {e}")
        return None

json_data = {
    'metadata': {'title': 'Metabolic Network Clustering'},
    'fig6': {
        'x': corr_ms_mg_clustered.columns.str.replace('_conc', '').tolist(),
        'y': corr_ms_mg_clustered.index.str.replace('_conc', '').tolist(),
        'z': corr_ms_mg_clustered.values.tolist(),
        'dendrogram': dendro_ms_mg,
        'title': 'MS & MG Correlation Matrix (Clustered)',
        'colorscale': [
            [0.0, '#ffffd9'],
            [0.2, '#c7e9b4'],
            [0.4, '#41b6c4'],
            [0.6, '#1d91c0'],
            [0.8, '#225ea8'],
            [1.0, '#0c2c84']
        ],
        'static_image': get_base64_image('Fig6_Corr_MS_MG.png')
    },
    'fig7': {
        'x': corr_ctrl_clustered.columns.str.replace('_conc', '').tolist(),
        'y': corr_ctrl_clustered.index.str.replace('_conc', '').tolist(),
        'z': corr_ctrl_clustered.values.tolist(),
        'dendrogram': dendro_ctrl,
        'title': 'Control Correlation Matrix (Clustered)',
        'colorscale': [
            [0.0, '#ffffd9'],
            [0.2, '#c7e9b4'],
            [0.4, '#41b6c4'],
            [0.6, '#1d91c0'],
            [0.8, '#225ea8'],
            [1.0, '#0c2c84']
        ],
        'static_image': get_base64_image('Fig7_Corr_Control.png')
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
