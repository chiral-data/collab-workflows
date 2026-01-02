"""Node 6: Differential Diagnosis Biomarkers - Generate Fig 4 & 5"""
import os
import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from matplotlib.patches import Patch

from html_generator import generate_biomarkers_html

# Configuration
INPUT_FILE = "../01_Data_Ingestion_and_Preprocessing/output/data_standardized.pkl"
OUTPUT_DIR = "output"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'sans-serif'

print(">>> NODE 6: DIFFERENTIAL DIAGNOSIS BIOMARKERS...", flush=True)

# Load data
df = pd.read_pickle(INPUT_FILE)

# Define masks
ms_types = ['RRMS', 'SPMS', 'PPMS']
mg_types = ['general', 'eye-type']
masks = {
    'MS': df['Type'].isin(ms_types),
    'MG': df['Type'].isin(mg_types)
}

# ==========================================
# FIG 4: Specific Amino Acids (CIT, GABA, AAA)
# ==========================================
print("Generating Fig 4: Specific Differences...", flush=True)

df_mg_ms = df[masks['MS'] | masks['MG']].copy()
df_mg_ms['Group'] = np.where(df_mg_ms['Type'].isin(['general', 'eye-type']), 'MG', 'MS')

fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=False)

amino_acids_fig4 = ['CIT_conc', 'GABA_conc', 'AAA_conc']

for i, aa in enumerate(amino_acids_fig4):
    ax = axes[i]
    sns.boxplot(data=df_mg_ms, x='Group', y=aa, palette={'MS':'salmon', 'MG':'teal'}, ax=ax)
    ax.set_title(aa.replace('_conc', ''), fontsize=12, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('Concentration [nmol/ml]' if i == 0 else '')
    
    # Add p-value
    ms_vals = df_mg_ms[df_mg_ms['Group']=='MS'][aa].dropna()
    mg_vals = df_mg_ms[df_mg_ms['Group']=='MG'][aa].dropna()
    _, pval = mannwhitneyu(ms_vals, mg_vals, alternative='two-sided')
    
    y_max = df_mg_ms[aa].max()
    y_min = df_mg_ms[aa].min()
    y_range = y_max - y_min
    ax.annotate(f'p = {pval:.2e}', xy=(0.5, y_max + 0.1*y_range), 
                fontsize=10, ha='center', fontweight='bold')

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig4_Specific_Diffs.png'))
plt.close()
print(f"Saved: Fig4_Specific_Diffs.png", flush=True)

# ==========================================
# FIG 5: Females Only
# ==========================================
print("Generating Fig 5: Female Specific...", flush=True)

df_fem = df_mg_ms[df_mg_ms['Sex'] == 'Female']
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for i, aa in enumerate(['CIT_conc', 'GABA_conc', 'AAA_conc']):
    ax = axes[i]
    sns.boxplot(data=df_fem, x='Group', y=aa, palette={'MS':'salmon', 'MG':'teal'}, ax=ax)
    ax.set_title(aa.replace('_conc', ''), fontsize=12, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('Concentration [nmol/ml]' if i == 0 else '')

legend_elements = [Patch(facecolor='salmon', edgecolor='black', label='MS'),
                   Patch(facecolor='teal', edgecolor='black', label='MG')]
fig.legend(handles=legend_elements, loc='upper right')

fig.suptitle('Specific Amino Acid Differences (Females Only)', fontsize=14, fontweight='bold')
plt.tight_layout(rect=[0, 0, 0.9, 0.95])
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig5_Female_Specific.png'))
plt.close()
print(f"Saved: Fig5_Female_Specific.png", flush=True)


# ==========================================
# GENERATE JSON DATA
# ==========================================
print("Generating JSON data...", flush=True)

json_data = {
    'metadata': {'title': 'Differential Diagnosis Biomarkers'},
    'fig4': {
        'traces': [
            {'name': 'MS', 'aa': 'ARG', 'y': df_mg_ms[df_mg_ms['Group']=='MS']['ARG_conc'].dropna().tolist(), 'color': '#fa8072'},
            {'name': 'MG', 'aa': 'ARG', 'y': df_mg_ms[df_mg_ms['Group']=='MG']['ARG_conc'].dropna().tolist(), 'color': '#008080'},
            {'name': 'MS', 'aa': 'PRO', 'y': df_mg_ms[df_mg_ms['Group']=='MS']['PRO_conc'].dropna().tolist(), 'color': '#fa8072'},
            {'name': 'MG', 'aa': 'PRO', 'y': df_mg_ms[df_mg_ms['Group']=='MG']['PRO_conc'].dropna().tolist(), 'color': '#008080'},
            {'name': 'MS', 'aa': 'CIT', 'y': df_mg_ms[df_mg_ms['Group']=='MS']['CIT_conc'].dropna().tolist(), 'color': '#fa8072'},
            {'name': 'MG', 'aa': 'CIT', 'y': df_mg_ms[df_mg_ms['Group']=='MG']['CIT_conc'].dropna().tolist(), 'color': '#008080'}
        ]
    },
    'fig5': {
        'traces': [
            {'name': 'MS', 'aa': 'ARG', 'y': df_fem[df_fem['Group']=='MS']['ARG_conc'].dropna().tolist(), 'color': '#fa8072'},
            {'name': 'MG', 'aa': 'ARG', 'y': df_fem[df_fem['Group']=='MG']['ARG_conc'].dropna().tolist(), 'color': '#008080'},
            {'name': 'MS', 'aa': 'PRO', 'y': df_fem[df_fem['Group']=='MS']['PRO_conc'].dropna().tolist(), 'color': '#fa8072'},
            {'name': 'MG', 'aa': 'PRO', 'y': df_fem[df_fem['Group']=='MG']['PRO_conc'].dropna().tolist(), 'color': '#008080'},
            {'name': 'MS', 'aa': 'CIT', 'y': df_fem[df_fem['Group']=='MS']['CIT_conc'].dropna().tolist(), 'color': '#fa8072'},
            {'name': 'MG', 'aa': 'CIT', 'y': df_fem[df_fem['Group']=='MG']['CIT_conc'].dropna().tolist(), 'color': '#008080'}
        ]
    }
}

json_path = os.path.join(OUTPUT_DIR, 'biomarkers_data.json')
with open(json_path, 'w', encoding='utf-8') as f:
    json.dump(json_data, f, indent=2)
print(f"Saved: {json_path}", flush=True)

html_content = generate_biomarkers_html(json_filename='biomarkers_data.json')

html_path = os.path.join(OUTPUT_DIR, 'biomarkers.html')
with open(html_path, 'w', encoding='utf-8') as f:
    f.write(html_content)
print(f"Saved: {html_path}", flush=True)

print("Node 6 completed successfully.", flush=True)

