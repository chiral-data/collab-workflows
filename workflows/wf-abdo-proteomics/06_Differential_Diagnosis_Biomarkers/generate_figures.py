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
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_FILE = "data_standardized.pkl"
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")

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
    sns.boxplot(data=df_mg_ms, x='Group', y=aa, hue='Group', palette={'MS':'salmon', 'MG':'teal'}, ax=ax, legend=False)
    ax.set_title(aa.replace('_conc', ''), fontsize=12, fontweight='bold')
    ax.set_xlabel('')
    ax.set_ylabel('Concentration [nmol/ml]' if i == 0 else '')
    
    # Add p-value
    ms_vals = df_mg_ms[df_mg_ms['Group']=='MS'][aa].dropna()
    mg_vals = df_mg_ms[df_mg_ms['Group']=='MG'][aa].dropna()
    if len(ms_vals) > 0 and len(mg_vals) > 0:
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

# Note: Using 'Female' as per existing data convention
df_fem = df_mg_ms[df_mg_ms['Sex'] == 'Female'].copy()
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for i, aa in enumerate(['CIT_conc', 'GABA_conc', 'AAA_conc']):
    ax = axes[i]
    sns.boxplot(data=df_fem, x='Group', y=aa, hue='Group', palette={'MS':'salmon', 'MG':'teal'}, ax=ax, legend=False)
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
# GENERATE ENRICHED JSON DATA
# ==========================================
print("Generating Enriched JSON data...", flush=True)

def create_enriched_data(df_source, amino_acids):
    """Generates data with statistical tests included"""
    subplots = []
    for aa in amino_acids:
        clean_name = aa.replace('_conc', '')
        
        # Extract Data
        ms_data = df_source[df_source['Group']=='MS'][aa].dropna()
        mg_data = df_source[df_source['Group']=='MG'][aa].dropna()
        
        # Calculate P-Value (Mann-Whitney U)
        p_val = 1.0
        if len(ms_data) > 0 and len(mg_data) > 0:
            stat, p_val = mannwhitneyu(ms_data, mg_data)
            
        subplots.append({
            'title': clean_name,
            'p_value': p_val,
            'traces': [
                {
                    'name': 'MS', 
                    'y': ms_data.tolist(), 
                    'color': '#fa8072' # Salmon
                },
                {
                    'name': 'MG', 
                    'y': mg_data.tolist(), 
                    'color': '#008080' # Teal
                }
            ]
        })
    return subplots

json_data = {
    'metadata': {'title': 'Differential Diagnosis Biomarkers'},
    'fig4': {
        'title': 'General Cohort (MS vs MG)',
        'subplots': create_enriched_data(df_mg_ms, amino_acids_fig4)
    },
    'fig5': {
        'title': 'Female Subgroup (MS vs MG)',
        'subplots': create_enriched_data(df_fem, amino_acids_fig4)
    }
}

json_path = os.path.join(OUTPUT_DIR, 'biomarkers_data.json')
with open(json_path, 'w', encoding='utf-8') as f:
    json.dump(json_data, f, indent=2)
print(f"Success! Enriched JSON saved to {json_path}", flush=True)

# ==========================================
# GENERATE HTML DASHBOARD
# ==========================================
html_content = generate_biomarkers_html(json_filename='biomarkers_data.json')

html_path = os.path.join(OUTPUT_DIR, 'biomarkers.html')
with open(html_path, 'w', encoding='utf-8') as f:
    f.write(html_content)
print(f"Saved: {html_path}", flush=True)

print("Node 6 completed successfully.", flush=True)
