"""Node 10: Clinical Mimicry Test - Generate Fig 10 & 11"""
import os
import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

from html_generator import generate_mimicry_html

# Configuration
INPUT_FILE = "data_standardized.pkl"
AA_COLS_FILE = "aa_cols.txt"
OUTPUT_DIR = "outputs"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'sans-serif'

print(">>> NODE 10: CLINICAL MIMICRY TEST (RRMS vs MG)...", flush=True)

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
# FIG 10: Total AA RRMS vs MG
# ==========================================
print("Generating Fig 10: Total AA RRMS vs MG...", flush=True)

df_rm = df[masks['RRMS'] | masks['MG']].copy()
df_rm['Group'] = np.where(df_rm['Type']=='RRMS', 'RRMS', 'MG')

fig, ax = plt.subplots(figsize=(6, 6))
sns.boxplot(data=df_rm, x='Group', y='Total_AA', palette={'RRMS':'grey', 'MG':'skyblue'}, ax=ax)
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
ax.annotate(f'p = {pval:.1e}', xy=(0.5, annotation_y_pos), fontsize=11, ha='center', fontweight='bold')
ax.set_ylim(y_ticks[0], annotation_y_pos * 1.1)

plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig10_TotalAA_RRMS_MG.png'), dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: Fig10_TotalAA_RRMS_MG.png", flush=True)

# ==========================================
# FIG 11: RRMS vs MG Grid
# ==========================================
print("Generating Fig 11: RRMS vs MG Grid...", flush=True)

df_melt_11 = df_rm.melt(id_vars='Group', value_vars=aa_cols, var_name='Amino Acid', value_name='Concentration')
df_melt_11['Amino Acid'] = df_melt_11['Amino Acid'].str.replace('_conc', '', regex=False)

g = sns.FacetGrid(df_melt_11, col='Amino Acid', col_wrap=6, sharey=False, height=2.5, aspect=1,
                  gridspec_kws={'hspace': 0.5, 'wspace': 0.3})

g.map_dataframe(sns.boxplot, x='Group', y='Concentration', palette={'RRMS':'lightgrey', 'MG':'lightblue'},
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
g.figure.suptitle('P2 Fig 9: RRMS vs MG', y=1.02, fontsize=16)

plt.tight_layout(rect=[0.03, 0, 1, 0.98])
g.savefig(os.path.join(OUTPUT_DIR, 'Fig11_RRMS_vs_MG_Grid.png'))
plt.close()
print(f"Saved: Fig11_RRMS_vs_MG_Grid.png", flush=True)


# ==========================================
# GENERATE JSON/HTML DATA
# ==========================================
print("Generating JSON and HTML...", flush=True)

# Helper for Fig 10
traces_fig10 = []
for group in ['RRMS', 'MG']:
    y_vals = df_rm[df_rm['Group'] == group]['Total_AA'].dropna().tolist()
    traces_fig10.append({
        'name': group,
        'y': y_vals,
        'color': 'grey' if group == 'RRMS' else 'skyblue'
    })

# Helper for Fig 11
subplots_fig11 = []
# aa_cols is already defined above
aa_names = [c.replace('_conc', '') for c in aa_cols]

for aa_col, aa_name in zip(aa_cols, aa_names):
    traces = []
    for group in ['RRMS', 'MG']:
        y_vals = df_rm[df_rm['Group'] == group][aa_col].dropna().tolist()
        traces.append({
            'name': group,
            'y': y_vals,
            'color': 'lightgrey' if group == 'RRMS' else 'lightblue'
        })
    subplots_fig11.append({
        'title': aa_name,
        'traces': traces
    })

json_data = {
    'metadata': {'title': 'RRMS vs MG Mimicry'},
    'fig10': {
        'title': 'Total Amino Acid Concentration: MS-RRMS vs MG',
        'yaxis': 'Concentration [nmol/ml]',
        'traces': traces_fig10
    },
    'fig11': {
        'subplots': subplots_fig11
    }
}

json_path = os.path.join(OUTPUT_DIR, 'mimicry_data.json')
with open(json_path, 'w', encoding='utf-8') as f:
    json.dump(json_data, f, indent=2)
print(f"Saved: {json_path}", flush=True)

html_content = generate_mimicry_html(json_filename='mimicry_data.json')

html_path = os.path.join(OUTPUT_DIR, 'mimicry.html')
with open(html_path, 'w', encoding='utf-8') as f:
    f.write(html_content)
print(f"Saved: {html_path}", flush=True)

print("Node 10 completed successfully.", flush=True)
