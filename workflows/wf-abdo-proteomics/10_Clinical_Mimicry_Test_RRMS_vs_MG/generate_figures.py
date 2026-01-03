"""Node 10: Clinical Mimicry Test - Generate Fig 10 & 11"""
import os
import sys
import json
import io
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

from html_generator import generate_mimicry_html

# Configuration
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_FILE = "data_standardized.pkl"
AA_COLS_FILE = "aa_cols.txt"
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'sans-serif'

print(">>> NODE 10: CLINICAL MIMICRY TEST (RRMS vs MG)...")

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
# FIG 10: Total AA RRMS vs MG (Static)
# ==========================================
print("Generating Fig 10: Total AA RRMS vs MG...")

df_rm = df[masks['RRMS'] | masks['MG']].copy()
df_rm['Group'] = np.where(df_rm['Type']=='RRMS', 'RRMS', 'MG')

fig, ax = plt.subplots(figsize=(6, 6))
sns.boxplot(data=df_rm, x='Group', y='Total_AA', palette={'RRMS':'grey', 'MG':'skyblue'}, ax=ax)
sns.boxplot(data=df_rm, x='Group', y='Total_AA', hue='Group', palette={'RRMS':'grey', 'MG':'skyblue'}, ax=ax, legend=False)
ax.set_ylabel('Concentration [nmol/ml]', fontsize=12)
ax.set_xlabel('', fontsize=12)
ax.set_title('Total Amino Acid Concentration: MS-RRMS vs MG', fontsize=12, fontweight='bold')

# Calculate p-value
rrms_data = df_rm[df_rm['Group']=='RRMS']['Total_AA'].dropna()
mg_data = df_rm[df_rm['Group']=='MG']['Total_AA'].dropna()
_, pval = mannwhitneyu(rrms_data, mg_data, alternative='two-sided')

y_ticks = [-5, -2.5, 0, 2.5, 5, 7.5, 10, 12.5]
ax.set_yticks(y_ticks)

y_max = df_rm['Total_AA'].max()
annotation_y_pos = max(y_max, y_ticks[-1]) * 1.05
ax.annotate('p = {:.1e}'.format(pval), xy=(0.5, annotation_y_pos), fontsize=11, ha='center', fontweight='bold')
ax.set_ylim(y_ticks[0], annotation_y_pos * 1.1)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig10_TotalAA_RRMS_MG.png'), dpi=150, bbox_inches='tight')
plt.close()

# ==========================================
# FIG 11: Grid (Static) - Keep for reference
# ==========================================
print("Generating Fig 11: Static Grid...")

df_melt_11 = df_rm.melt(id_vars='Group', value_vars=aa_cols, var_name='Amino Acid', value_name='Concentration')
df_melt_11['Amino Acid'] = df_melt_11['Amino Acid'].str.replace('_conc', '', regex=False)

g = sns.FacetGrid(df_melt_11, col='Amino Acid', col_wrap=6, sharey=False, height=2.5, aspect=1,
                  gridspec_kws={'hspace': 0.5, 'wspace': 0.3})

g.map_dataframe(sns.boxplot, x='Group', y='Concentration', palette={'RRMS':'grey', 'MG':'skyblue'},
                hue='Group', legend=False,
                showfliers=True, 
                flierprops={'marker': 'o', 'markerfacecolor': 'black', 'markersize': 3})

g.set_titles("{col_name}")
g.set_xlabels("")

y_ticks = [-4, -2, 0, 2, 4, 6, 8, 10, 12]
for i, ax in enumerate(g.axes.flat):
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)
    ax.set_yticks(y_ticks)
    if i % 6 != 0:
        ax.set_yticklabels([])
    ax.set_ylabel('')

g.figure.supylabel('Concentration [nmol/ml]')
g.figure.suptitle('Fig 11: RRMS vs MG Grid', y=1.02, fontsize=16)

plt.tight_layout(rect=[0.03, 0, 1, 0.98])
g.savefig(os.path.join(OUTPUT_DIR, 'Fig11_RRMS_vs_MG_Grid.png'))
plt.close()
print("Saved: Fig11_RRMS_vs_MG_Grid.png")

# ==========================================
# GENERATE JSON DATA (Node 4 Style)
# ==========================================
print("Generating JSON data...")

json_data = {
    'metadata': {
        'title': 'Clinical Mimicry Test: RRMS vs MG',
        'amino_acids': [aa.replace('_conc', '') for aa in aa_cols]
    },
    'fig10': { # Hero Plot
        'title': 'Total Amino Acid Load',
        'variable': 'Total AA',
        'traces': [] 
    },
    'fig11': { # Grid Plot
        'type': 'box',
        'variable': 'Amino Acid',
        'traces': []
    }
}

# Populate Fig 10 (Total)
json_data['fig10']['traces'] = [{
    'name': 'Total AA',
    'RRMS': {'y': rrms_data.tolist()},
    'MG': {'y': mg_data.tolist()},
    'stats': {'p_value': pval}
}]

# Populate Fig 11 (Grid)
for aa in aa_cols:
    aa_clean = aa.replace('_conc', '')
    
    rrms_vals = df_rm[df_rm['Group']=='RRMS'][aa].dropna()
    mg_vals = df_rm[df_rm['Group']=='MG'][aa].dropna()
    
    # Calc stats
    p_val_aa = 1.0
    if len(rrms_vals) > 0 and len(mg_vals) > 0:
        _, p_val_aa = mannwhitneyu(rrms_vals, mg_vals, alternative='two-sided')
        
    json_data['fig11']['traces'].append({
        'aa': aa_clean,
        'RRMS': {'y': rrms_vals.tolist()},
        'MG': {'y': mg_vals.tolist()},
        'stats': {'p_value': p_val_aa}
    })

json_path = os.path.join(OUTPUT_DIR, 'mimicry_data.json')
with io.open(json_path, 'w', encoding='utf-8') as f:
    f.write(json.dumps(json_data, ensure_ascii=False))

# Generate HTML
html_content = generate_mimicry_html(json_filename='mimicry_data.json')

html_path = os.path.join(OUTPUT_DIR, 'mimicry.html')
with io.open(html_path, 'w', encoding='utf-8') as f:
    f.write(html_content)

print("Node 10 completed successfully.")
