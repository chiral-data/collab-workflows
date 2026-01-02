"""Node 4: Clinical Confounders & Validation - Generate Fig 2A, 2B, 2C"""
import os
import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D

from html_generator import generate_confounders_html

# Configuration
INPUT_FILE = "../01_Data_Ingestion_and_Preprocessing/output/data_standardized.pkl"
AA_COLS_FILE = "../01_Data_Ingestion_and_Preprocessing/output/aa_cols.txt"
OUTPUT_DIR = "output"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'sans-serif'

print(">>> NODE 4: CLINICAL CONFOUNDERS & VALIDATION...", flush=True)

# Load data
df = pd.read_pickle(INPUT_FILE)
with open(AA_COLS_FILE, 'r') as f:
    aa_cols = [line.strip() for line in f.readlines()]

# Define masks
ms_types = ['RRMS', 'SPMS', 'PPMS']
masks = {
    'MS': df['Type'].isin(ms_types),
    'Control': df['Status'] == 'control'
}

df_p1a = df[masks['MS'] | masks['Control']].copy()
df_p1a['Group'] = np.where(df_p1a['Status']=='control', 'Control', 'MS')

# ==========================================
# FIG 2A: Age Grid
# ==========================================
print("Generating Fig 2A: Age vs AA...", flush=True)

df_melt_2a = df_p1a.melt(id_vars=['Age', 'Group'], value_vars=aa_cols, var_name='Amino Acid', value_name='LogC')
df_melt_2a['Amino Acid'] = df_melt_2a['Amino Acid'].str.replace('_conc', '', regex=False)

g = sns.lmplot(data=df_melt_2a, x='Age', y='LogC', col='Amino Acid', hue='Group', 
               palette={'MS':'blue', 'Control':'red'},
               col_wrap=6, height=2.5, aspect=1, 
               scatter_kws={'s': 10, 'alpha': 0.5}, line_kws={'lw': 2},
               legend=False)
               
g.set_titles("{col_name}")
g.set(xticks=[25, 50, 75], yticks=[0, 10])

for i, ax in enumerate(g.axes.flat):
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)

g.figure.suptitle('Fig 2A: Age vs AA', y=1.02, fontsize=16)

legend_elements = [Line2D([0], [0], marker='o', color='w', label='MS', markerfacecolor='blue', markersize=8),
                   Line2D([0], [0], marker='o', color='w', label='Control', markerfacecolor='red', markersize=8)]
g.figure.legend(handles=legend_elements, title='Group', loc='lower right', bbox_to_anchor=(0.98, 0.05))

plt.tight_layout(rect=[0, 0, 1, 0.98])
g.savefig(os.path.join(OUTPUT_DIR, 'Fig2A_Age_Grid.png'))
plt.close()
print(f"Saved: Fig2A_Age_Grid.png", flush=True)

# ==========================================
# FIG 2B: Duration Grid
# ==========================================
print("Generating Fig 2B: Duration vs AA...", flush=True)

df_p2b = df[masks['MS']]
df_melt = df_p2b.melt(id_vars='Duration', value_vars=aa_cols, var_name='Amino Acid', value_name='LogC')
df_melt['Amino Acid'] = df_melt['Amino Acid'].str.replace('_conc', '', regex=False)

g = sns.lmplot(data=df_melt, x='Duration', y='LogC', col='Amino Acid',
               col_wrap=6, height=2.5, aspect=1,
               scatter_kws={'s': 10, 'alpha': 0.5}, line_kws={'lw': 2, 'color': 'blue'})
               
g.set_titles("{col_name}")
g.set(xticks=[0, 20, 40], yticks=[-5, 0, 5])

for i, ax in enumerate(g.axes.flat):
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)

g.figure.suptitle('Fig 2B: Duration vs AA', y=1.02, fontsize=16)

plt.tight_layout(rect=[0, 0, 1, 0.98])
g.savefig(os.path.join(OUTPUT_DIR, 'Fig2B_Duration_Grid.png'))
plt.close()
print(f"Saved: Fig2B_Duration_Grid.png", flush=True)

# ==========================================
# FIG 2C: EDSS Grid
# ==========================================
print("Generating Fig 2C: EDSS vs AA...", flush=True)

df_p2c = df[masks['MS']]
df_melt = df_p2c.melt(id_vars='EDSS', value_vars=aa_cols, var_name='Amino Acid', value_name='LogC')
df_melt['Amino Acid'] = df_melt['Amino Acid'].str.replace('_conc', '', regex=False)

g = sns.lmplot(data=df_melt, x='EDSS', y='LogC', col='Amino Acid',
               col_wrap=6, height=2.5, aspect=1,
               scatter_kws={'s': 10, 'alpha': 0.5}, line_kws={'lw': 2, 'color': 'blue'})
               
g.set_titles("{col_name}")
g.set(xticks=[0, 5, 10], yticks=[-5, 0, 5])

for i, ax in enumerate(g.axes.flat):
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)

g.figure.suptitle('Fig 2C: EDSS vs AA', y=1.02, fontsize=16)

plt.tight_layout(rect=[0, 0, 1, 0.98])
g.savefig(os.path.join(OUTPUT_DIR, 'Fig2C_EDSS_Grid.png'))
plt.close()
print(f"Saved: Fig2C_EDSS_Grid.png", flush=True)

# ==========================================
# GENERATE JSON DATA
# ==========================================
print("Generating JSON data...", flush=True)

json_data = {
    'metadata': {
        'title': 'Clinical Confounders & Validation',
        'amino_acids': [aa.replace('_conc', '') for aa in aa_cols]
    },
    'fig2a': {'type': 'scatter', 'variable': 'Age', 'traces': []},
    'fig2b': {'type': 'scatter', 'variable': 'Duration', 'traces': []},
    'fig2c': {'type': 'scatter', 'variable': 'EDSS', 'traces': []}
}

# Collect data for each amino acid
for aa in aa_cols:
    aa_clean = aa.replace('_conc', '')
    
    # Fig 2A: Age
    ms_data = df_p1a[df_p1a['Group']=='MS'][[aa, 'Age']].dropna()
    ctrl_data = df_p1a[df_p1a['Group']=='Control'][[aa, 'Age']].dropna()
    json_data['fig2a']['traces'].append({
        'aa': aa_clean,
        'MS': {'x': ms_data['Age'].tolist(), 'y': ms_data[aa].tolist()},
        'Control': {'x': ctrl_data['Age'].tolist(), 'y': ctrl_data[aa].tolist()}
    })
    
    # Fig 2B: Duration (MS only)
    dur_data = df_p1a[df_p1a['Group']=='MS'][[aa, 'Duration']].dropna()
    json_data['fig2b']['traces'].append({
        'aa': aa_clean,
        'MS': {'x': dur_data['Duration'].tolist(), 'y': dur_data[aa].tolist()}
    })
    
    # Fig 2C: EDSS (MS only)
    edss_data = df_p1a[df_p1a['Group']=='MS'][[aa, 'EDSS']].dropna()
    json_data['fig2c']['traces'].append({
        'aa': aa_clean,
        'MS': {'x': edss_data['EDSS'].tolist(), 'y': edss_data[aa].tolist()}
    })

json_path = os.path.join(OUTPUT_DIR, 'confounders_data.json')
with open(json_path, 'w', encoding='utf-8') as f:
    json.dump(json_data, f, indent=2, ensure_ascii=False)
print(f"Saved: {json_path}", flush=True)

# ==========================================
# GENERATE HTML
# ==========================================
html_content = generate_confounders_html(json_filename='confounders_data.json')

html_path = os.path.join(OUTPUT_DIR, 'confounders.html')
with open(html_path, 'w', encoding='utf-8') as f:
    f.write(html_content)
print(f"Saved: {html_path}", flush=True)

print("Node 4 completed successfully.", flush=True)
