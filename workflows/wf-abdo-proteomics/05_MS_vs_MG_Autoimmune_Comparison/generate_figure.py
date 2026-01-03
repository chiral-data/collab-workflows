"""Node 5: MS vs MG Autoimmune Comparison - Generate Enriched JSON"""
import os
import sys
import json
import pandas as pd
import numpy as np
from scipy import stats  # Required for statistical tests

# Configuration
INPUT_FILE = "data_standardized.pkl"
AA_COLS_FILE = "aa_cols.txt"
OUTPUT_DIR = "outputs"
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

print(">>> NODE 5: MS VS MG AUTOIMMUNE COMPARISON (Enriched)...", flush=True)
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

# Create a combined dataframe for easy processing
df_p2 = df.copy()
# Create the specific groups for this analysis
df_p2['Group'] = 'Other'
df_p2.loc[masks['MS'] | masks['MG'], 'Group'] = 'MS+MG'
df_p2.loc[masks['Control'], 'Group'] = 'Controls'

# Initialize JSON structure
json_data = {
    'metadata': {
        'title': 'MS+MG vs Controls',
        'amino_acids': [aa.split(' ')[0] for aa in aa_cols]
    },
    'fig3': {
        'subplots': []
    }
}

# Generate Data
for aa in aa_cols:
    aa_clean = aa.split(' ')[0] # Clean name
    
    # Extract data for the two groups
    ms_mg_data = df_p2[df_p2['Group']=='MS+MG'][aa].dropna()
    ctrl_data = df_p2[df_p2['Group']=='Controls'][aa].dropna()
    
    # Calculate Statistics (T-test)
    # We use equal_var=False (Welch's t-test) which is safer for biological data
    p_val = 1.0
    if len(ms_mg_data) > 1 and len(ctrl_data) > 1:
        t_stat, p_val = stats.ttest_ind(ms_mg_data, ctrl_data, equal_var=False)
        
    ms_mg_vals = ms_mg_data.tolist()
    ctrl_vals = ctrl_data.tolist()
    
    json_data['fig3']['subplots'].append({
    'title': aa_clean,
    'p_value': p_val,  # Save p-value for the web to use
    'traces': [
            {'name': 'MS+MG', 'y': ms_mg_vals, 'color': '#e74c3c'},   # Red
            {'name': 'Controls', 'y': ctrl_vals, 'color': '#2ecc71'}  # Green
        ]
    })

# Save JSON
json_path = os.path.join(OUTPUT_DIR, 'autoimmune_data.json')
with open(json_path, 'w', encoding='utf-8') as f:
    json.dump(json_data, f, indent=2)

print(f"Success! Enriched JSON saved to {json_path}")
"""Node 5: MS vs MG Autoimmune Comparison - Generate Enriched JSON"""
import os
import sys
import json
import pandas as pd
import numpy as np
from scipy import stats  # Required for statistical tests

# Configuration
INPUT_FILE = "data_standardized.pkl"
AA_COLS_FILE = "aa_cols.txt"
OUTPUT_DIR = "outputs"
if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

print(">>> NODE 5: MS VS MG AUTOIMMUNE COMPARISON (Enriched)...", flush=True)
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

# Create a combined dataframe for easy processing
df_p2 = df.copy()
# Create the specific groups for this analysis
df_p2['Group'] = 'Other'
df_p2.loc[masks['MS'] | masks['MG'], 'Group'] = 'MS+MG'
df_p2.loc[masks['Control'], 'Group'] = 'Controls'

# Initialize JSON structure
json_data = {
    'metadata': {
        'title': 'MS+MG vs Controls',
        'amino_acids': [aa.split(' ')[0] for aa in aa_cols]
    },
    'fig3': {
        'subplots': []
    }
}

# Generate Data
for aa in aa_cols:
    aa_clean = aa.split(' ')[0] # Clean name
    
    # Extract data for the two groups
    ms_mg_data = df_p2[df_p2['Group']=='MS+MG'][aa].dropna()
    ctrl_data = df_p2[df_p2['Group']=='Controls'][aa].dropna()
    
    # Calculate Statistics (T-test)
    # We use equal_var=False (Welch's t-test) which is safer for biological data
    p_val = 1.0
    if len(ms_mg_data) > 1 and len(ctrl_data) > 1:
        t_stat, p_val = stats.ttest_ind(ms_mg_data, ctrl_data, equal_var=False)
        
    ms_mg_vals = ms_mg_data.tolist()
    ctrl_vals = ctrl_data.tolist()
    
    json_data['fig3']['subplots'].append({
        'title': aa_clean,
        'p_value': p_val,  # Save p-value for the web to use
        'traces': [
            {'name': 'MS+MG', 'y': ms_mg_vals, 'color': '#e74c3c'},   # Red
            {'name': 'Controls', 'y': ctrl_vals, 'color': '#2ecc71'}  # Green
        ]
    })

# Save JSON
json_path = os.path.join(OUTPUT_DIR, 'autoimmune_data.json')
with open(json_path, 'w', encoding='utf-8') as f:
    json.dump(json_data, f, indent=2)

print(f"Success! Enriched JSON saved to {json_path}")
"""Node 5: MS vs MG Autoimmune Comparison - Generate Fig 3"""
import os
import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from html_generator import generate_autoimmune_html

# Configuration
INPUT_FILE = "data_standardized.pkl"
AA_COLS_FILE = "aa_cols.txt"
OUTPUT_DIR = "outputs"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'sans-serif'

print(">>> NODE 5: MS VS MG AUTOIMMUNE COMPARISON...", flush=True)

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

# ==========================================
# FIG 3: MS+MG vs Control
# ==========================================
print("Generating Fig 3: MS+MG vs Controls...", flush=True)

df_p2 = df[masks['MS'] | masks['MG'] | masks['Control']].copy()
df_p2['Group'] = np.where(df_p2['Status']=='control', 'Controls', 'MS+MG')

df_melt = df_p2.melt(id_vars='Group', value_vars=aa_cols, var_name='Amino Acid', value_name='LogC')
df_melt['Amino Acid'] = df_melt['Amino Acid'].str.replace('_conc', '', regex=False)

g = sns.FacetGrid(df_melt, col='Amino Acid', col_wrap=6, sharey=False, height=2.5, aspect=1,
                  gridspec_kws={'hspace': 0.5, 'wspace': 0.3})

g.map_dataframe(sns.boxplot, x='Group', y='LogC', palette={'MS+MG':'skyblue', 'Controls':'grey'}, 
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
g.figure.suptitle('P2 Fig 1: MS+MG vs Controls', y=1.02, fontsize=16)

plt.tight_layout(rect=[0.03, 0, 1, 0.98])
g.savefig(os.path.join(OUTPUT_DIR, 'Fig3_MS_MG_vs_Control.png'))
plt.close()
print(f"Saved: Fig3_MS_MG_vs_Control.png", flush=True)


# ==========================================
# GENERATE JSON DATA
# ==========================================
print("Generating JSON data...", flush=True)

json_data = {
    'metadata': {
        'title': 'MS+MG vs Controls',
        'amino_acids': [aa.replace('_conc', '') for aa in aa_cols]
    },
    'fig3': {
        'subplots': []
    }
}

# Fig 3 data
for aa in aa_cols:
    aa_clean = aa.replace('_conc', '')
    ms_mg_vals = df_p2[df_p2['Group']=='MS+MG'][aa].dropna().tolist()
    ctrl_vals = df_p2[df_p2['Group']=='Controls'][aa].dropna().tolist()
    
    json_data['fig3']['subplots'].append({
        'title': aa_clean,
        'traces': [
            {'name': 'MS+MG', 'y': ms_mg_vals, 'color': '#87ceeb'},
            {'name': 'Controls', 'y': ctrl_vals, 'color': '#808080'}
        ]
    })

json_path = os.path.join(OUTPUT_DIR, 'autoimmune_data.json')
with open(json_path, 'w', encoding='utf-8') as f:
    json.dump(json_data, f, indent=2, ensure_ascii=False)
print(f"Saved: {json_path}", flush=True)

# ==========================================
# GENERATE HTML
# ==========================================
html_content = generate_autoimmune_html(json_filename='autoimmune_data.json')

html_path = os.path.join(OUTPUT_DIR, 'autoimmune.html')
with open(html_path, 'w', encoding='utf-8') as f:
    f.write(html_content)
print(f"Saved: {html_path}", flush=True)

print("Node 5 completed successfully.", flush=True)
