"""Node 3: MS Pathology Overview - Generate Fig 1A & 1B"""
import os
import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Import local HTML generator
from html_generator import generate_pathology_html

# Configuration
INPUT_FILE = "data_standardized.pkl"
AA_COLS_FILE = "aa_cols.txt"
OUTPUT_DIR = "outputs"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'sans-serif'

print(">>> NODE 3: MS PATHOLOGY OVERVIEW...", flush=True)

# Load data
df = pd.read_pickle(INPUT_FILE)
with open(AA_COLS_FILE, 'r') as f:
    aa_cols = [line.strip() for line in f.readlines()]

print(f"Loaded {len(df)} records", flush=True)

# Define masks
ms_types = ['RRMS', 'SPMS', 'PPMS']
masks = {
    'MS': df['Type'].isin(ms_types),
    'Control': df['Status'] == 'control'
}

# ==========================================
# FIG 1A: MS vs Control (29 Grid)
# ==========================================
print("Generating Fig 1A: MS vs Control...", flush=True)

df_p1a = df[masks['MS'] | masks['Control']].copy()
df_p1a['Group'] = np.where(df_p1a['Status']=='control', 'Control', 'MS')

df_melt = df_p1a.melt(id_vars='Group', value_vars=aa_cols, var_name='Amino Acid', value_name='LogC')
df_melt['Amino Acid'] = df_melt['Amino Acid'].str.replace('_conc', '', regex=False)

g = sns.FacetGrid(df_melt, col='Amino Acid', col_wrap=6, sharey=False, height=2.5, aspect=1,
                  gridspec_kws={'hspace': 0.5, 'wspace': 0.3})
g.map_dataframe(sns.boxplot, x='Group', y='LogC', palette={'MS':'steelblue', 'Control':'salmon'}, 
                hue='Group', palette={'MS':'steelblue', 'Control':'salmon'}, legend=False,
                showfliers=True, 
                flierprops={'marker': 'o', 'markerfacecolor': 'black', 'markersize': 3})

g.set_titles("{col_name}")
g.set_xlabels("")

y_ticks = [0, 10]
for i, ax in enumerate(g.axes.flat):
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)
    ax.set_yticks(y_ticks)
    if i % 6 != 0:
        ax.set_yticklabels([])
    ax.set_ylabel('')

g.figure.supylabel('Log C')
g.figure.suptitle('Fig 1A: MS vs Control', y=1.02, fontsize=16)
plt.tight_layout(rect=[0.03, 0, 1, 0.98])
g.savefig(os.path.join(OUTPUT_DIR, 'Fig1A_MS_vs_Control.png'))
plt.close()
print(f"Saved: Fig1A_MS_vs_Control.png", flush=True)

# ==========================================
# FIG 1B: MS Subtypes (29 Grid)
# ==========================================
print("Generating Fig 1B: MS Subtypes...", flush=True)

df_p1b = df[masks['MS']]
df_melt = df_p1b.melt(id_vars='Type', value_vars=aa_cols, var_name='Amino Acid', value_name='LogC')
df_melt['Amino Acid'] = df_melt['Amino Acid'].str.replace('_conc', '', regex=False)

g = sns.FacetGrid(df_melt, col='Amino Acid', col_wrap=5, sharey=False, height=2.5, aspect=1,
                  gridspec_kws={'hspace': 0.5, 'wspace': 0.3})
g.map_dataframe(sns.boxplot, x='Type', y='LogC', 
                order=['RRMS', 'SPMS', 'PPMS'],
                palette={'RRMS':'lightgreen', 'SPMS':'orange', 'PPMS':'purple'},
                showfliers=True,
                flierprops={'marker': 'o', 'markerfacecolor': 'black', 'markersize': 3})

g.set_titles("{col_name}")
g.set_xlabels("")

y_ticks = [-5, 0, 5]
for i, ax in enumerate(g.axes.flat):
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)
    ax.set_yticks(y_ticks)
    if i % 5 != 0:
        ax.set_yticklabels([])
    ax.set_ylabel('')

g.figure.supylabel('Log C')
g.figure.suptitle('Fig 1B: MS Subtypes', y=1.02, fontsize=16)
plt.tight_layout(rect=[0.03, 0, 1, 0.98])
g.savefig(os.path.join(OUTPUT_DIR, 'Fig1B_MS_Subtypes.png'))
plt.close()
print(f"Saved: Fig1B_MS_Subtypes.png", flush=True)

# ==========================================
# GENERATE JSON DATA
# ==========================================
print("Generating JSON data...", flush=True)

json_data = {
    'metadata': {
        'title': 'MS Pathology Overview',
        'amino_acids': [aa.replace('_conc', '') for aa in aa_cols]
    },
    'fig1a': {
        'subplots': []
    },
    'fig1b': {
        'subplots': []
    }
}

# Fig 1A data
for aa in aa_cols:
    aa_clean = aa.replace('_conc', '')
    ms_vals = df_p1a[df_p1a['Group']=='MS'][aa].dropna().tolist()
    ctrl_vals = df_p1a[df_p1a['Group']=='Control'][aa].dropna().tolist()
    
    json_data['fig1a']['subplots'].append({
        'title': aa_clean,
        'traces': [
            {'name': 'MS', 'y': ms_vals, 'color': '#4682b4'},
            {'name': 'Control', 'y': ctrl_vals, 'color': '#fa8072'}
        ]
    })

# Fig 1B data  
for aa in aa_cols:
    aa_clean = aa.replace('_conc', '')
    rrms_vals = df_p1b[df_p1b['Type']=='RRMS'][aa].dropna().tolist()
    spms_vals = df_p1b[df_p1b['Type']=='SPMS'][aa].dropna().tolist()
    ppms_vals = df_p1b[df_p1b['Type']=='PPMS'][aa].dropna().tolist()
    
    json_data['fig1b']['subplots'].append({
        'title': aa_clean,
        'traces': [
            {'name': 'RRMS', 'y': rrms_vals, 'color': '#90ee90'},
            {'name': 'SPMS', 'y': spms_vals, 'color': '#ffa500'},
            {'name': 'PPMS', 'y': ppms_vals, 'color': '#800080'}
        ]
    })

json_path = os.path.join(OUTPUT_DIR, 'pathology_data.json')
with open(json_path, 'w', encoding='utf-8') as f:
    json.dump(json_data, f, indent=2, ensure_ascii=False)
print(f"Saved: {json_path}", flush=True)

# ==========================================
# GENERATE HTML
# ==========================================
html_content = generate_pathology_html(json_filename="pathology_data.json")

html_path = os.path.join(OUTPUT_DIR, 'pathology.html')
with open(html_path, 'w', encoding='utf-8') as f:
    f.write(html_content)
print(f"Saved: {html_path}", flush=True)

print("Node 3 completed successfully.", flush=True)

