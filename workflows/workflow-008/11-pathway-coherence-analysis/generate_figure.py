"""Node 11: Pathway Coherence Analysis - Generate Fig 12"""
import os
import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from html_generator import generate_coherence_html

# Configuration
INPUT_FILE = "inputs/data_standardized.pkl"
AA_COLS_FILE = "inputs/aa_cols.txt"
OUTPUT_DIR = "outputs"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'sans-serif'

print(">>> NODE 11: PATHWAY COHERENCE ANALYSIS...", flush=True)

# Load data
df = pd.read_pickle(INPUT_FILE)
with open(AA_COLS_FILE, 'r') as f:
    aa_cols = [line.strip() for line in f.readlines()]

# Define masks
mg_types = ['general', 'eye-type']
masks = {
    'RRMS': df['Type'] == 'RRMS',
    'MG': df['Type'].isin(mg_types)
}

# ==========================================
# FIG 12: Split Correlation (MG vs RRMS)
# ==========================================
print("Generating Fig 12: Split Correlation Matrix...", flush=True)

fig, ax = plt.subplots(1, 2, figsize=(20, 8))

aa_labels = [col.replace('_conc', '') for col in aa_cols]

# MG correlation matrix (LEFT)
corr_mg = df[masks['MG']][aa_cols].corr()
corr_mg.index = aa_labels
corr_mg.columns = aa_labels
mask_mg = np.triu(np.ones_like(corr_mg, dtype=bool))

sns.heatmap(corr_mg, ax=ax[0], cmap='YlGnBu', mask=mask_mg, annot=True, fmt='.2f', 
            annot_kws={'size': 6}, cbar=True, vmin=-0.2, vmax=1.0, square=True,
            cbar_kws={'shrink': 0.8, 'label': 'Correlation'})
ax[0].set_title('MG', fontsize=14, fontweight='bold')
ax[0].tick_params(labelsize=7)

# RRMS correlation matrix (RIGHT)
corr_rrms = df[masks['RRMS']][aa_cols].corr()
corr_rrms.index = aa_labels
corr_rrms.columns = aa_labels
mask_rrms = np.triu(np.ones_like(corr_rrms, dtype=bool))

sns.heatmap(corr_rrms, ax=ax[1], cmap='YlGnBu', mask=mask_rrms, annot=True, fmt='.2f',
            annot_kws={'size': 6}, cbar=True, vmin=-0.2, vmax=1.0, square=True,
            cbar_kws={'shrink': 0.8, 'label': 'Correlation'})
ax[1].set_title('MS-RRMS', fontsize=14, fontweight='bold')
ax[1].tick_params(labelsize=7)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig12_Split_Corr.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: Fig12_Split_Corr.png", flush=True)


# ==========================================
# GENERATE JSON/HTML DATA
# ==========================================
print("Generating JSON and HTML...", flush=True)

from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list

# ... (Imports preserved above)

# Calculate dictionaries for all cohorts
cohorts = {
    'Control': df['Status'] == 'control',
    'RRMS': df['Type'] == 'RRMS',
    'SPMS': df['Type'] == 'SPMS',
    'PPMS': df['Type'] == 'PPMS',
    'GMG': df['Type'] == 'general',
    'OMG': df['Type'] == 'eye-type'
}

correlation_data = {}

for name, mask in cohorts.items():
    subset = df[mask][aa_cols]
    if len(subset) < 2:
        continue
        
    # Calculate Correlation
    corr = subset.corr()
    # Ensure standard index/columns
    corr.index = aa_labels
    corr.columns = aa_labels
    
    # Perform Clustering using Seaborn Clustermap (to get order and linkage)
    # We use a temporary figure to extract data, then close it.
    try:
        # Clustermap (Standard call, no invalid args)
        g = sns.clustermap(corr, cmap='YlGnBu', 
                           method='average', metric='euclidean') 
        
        # Reorder correlation matrix
        reordered_ind_row = g.dendrogram_row.reordered_ind
        reordered_ind_col = g.dendrogram_col.reordered_ind
        
        corr_clustered = corr.iloc[reordered_ind_row, reordered_ind_col]
        
        # Extract Dendrogram Data
        dendro = {'row': None, 'col': None}
        
        # Row Dendrogram
        if hasattr(g, 'dendrogram_row') and g.dendrogram_row is not None:
             d = dendrogram(g.dendrogram_row.linkage, no_plot=True)
             icoords = [[(x - 5.0) / 10.0 for x in seg] for seg in d['icoord']]
             dendro['row'] = {'icoords': icoords, 'dcoords': d['dcoord']}
             
        # Col Dendrogram
        if hasattr(g, 'dendrogram_col') and g.dendrogram_col is not None:
             d = dendrogram(g.dendrogram_col.linkage, no_plot=True)
             icoords = [[(x - 5.0) / 10.0 for x in seg] for seg in d['icoord']]
             dendro['col'] = {'icoords': icoords, 'dcoords': d['dcoord']}

        plt.close(g.figure) # Close to save memory
        
    except Exception as e:
        print(f"Warning: Clustering failed for {name}, falling back to unclustered. Error: {e}")
        corr_clustered = corr
        dendro = None

    # Mask Upper Triangle AND Diagonal (Keep strictly lower)
    # np.triu(..., k=0) includes diagonal. We want to mask that.
    mask = np.triu(np.ones_like(corr_clustered, dtype=bool), k=0)
    corr_masked = corr_clustered.mask(mask)

    # Prepare Z-values (handle NaNs for JSON)
    z_vals = corr_masked.where(pd.notnull(corr_masked), None).values.tolist()
    
    correlation_data[name] = {
        'x': corr_clustered.columns.tolist(),
        'y': corr_clustered.index.tolist(),
        'z': z_vals,
        'dendrogram': dendro,
        'title': f'{name} Clustered Network'
    }

json_data = {
    'metadata': {'title': 'Pathway Coherence Comparator'},
    'cohorts': correlation_data
}

json_path = os.path.join(OUTPUT_DIR, 'coherence_data.json')
with open(json_path, 'w', encoding='utf-8') as f:
    json.dump(json_data, f, indent=2)
print(f"Saved: {json_path}", flush=True)

html_content = generate_coherence_html(json_filename='coherence_data.json')

html_path = os.path.join(OUTPUT_DIR, 'coherence.html')
with open(html_path, 'w', encoding='utf-8') as f:
    f.write(html_content)
print(f"Saved: {html_path}", flush=True)

print("Node 11 completed successfully.", flush=True)
