"""Node 9: Subtype Trajectories - Generate Fig 9"""
import os
import sys
import json
import io
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

from html_generator import generate_trajectories_html

# Configuration
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_FILE = "data_standardized.pkl"
AA_COLS_FILE = "aa_cols.txt"
OUTPUT_DIR = "outputs"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

plt.rcParams['font.family'] = 'sans-serif'

# Load data
df = pd.read_pickle(INPUT_FILE)
with open(AA_COLS_FILE, 'r') as f:
    aa_cols = [line.strip() for line in f.readlines()]

# Define masks
ms_types = ['RRMS', 'SPMS', 'PPMS']
mg_types = ['general', 'eye-type']
masks = {
    'MS': df['Type'].isin(ms_types),
    'MG': df['Type'].isin(mg_types)
}

df_mg_ms = df[masks['MS'] | masks['MG']].copy()

# ==========================================
# FIG 9: Duration Grid by Subtype (Static Image)
# ==========================================
# (Keep existing static generation as user liked it)

df_fig9 = df_mg_ms.copy()

def get_subtype(row):
    if row['Type'] == 'PPMS':
        return 'PPMS'
    elif row['Type'] == 'SPMS':
        return 'SPMS'
    elif row['Type'] == 'RRMS':
        return 'RRMS'
    elif row['Type'] == 'general':
        return 'GMG'
    elif row['Type'] == 'eye-type':
        return 'OMG'
    return 'Unknown'

df_fig9['Subtype'] = df_fig9.apply(get_subtype, axis=1)

subtype_colors = {'PPMS': '#E74C3C', 'SPMS': '#3498DB', 'RRMS': '#2ECC71', 'GMG': '#9B59B6', 'OMG': '#F39C12'}
subtype_order = ['PPMS', 'SPMS', 'RRMS', 'GMG', 'OMG']

aa_clean_names = [aa.replace('_conc', '') for aa in aa_cols]

n_cols = 10
n_rows = 3
fig, axes = plt.subplots(n_rows, n_cols, figsize=(24, 9))
axes = axes.flatten()

duration_min = df_fig9['Duration'].min()
duration_max = df_fig9['Duration'].max()

for idx, (aa_col, aa_name) in enumerate(zip(aa_cols, aa_clean_names)):
    ax = axes[idx]
    
    for subtype in subtype_order:
        subtype_data = df_fig9[df_fig9['Subtype'] == subtype]
        if len(subtype_data) > 0:
            ax.scatter(subtype_data['Duration'], subtype_data[aa_col], 
                      color=subtype_colors[subtype], label=subtype, alpha=0.6, s=30)
            
            if len(subtype_data) > 1:
                valid_data = subtype_data[['Duration', aa_col]].dropna()
                if len(valid_data) > 1:
                    slope, intercept, r_val, p_val, std_err = linregress(valid_data['Duration'], valid_data[aa_col])
                    x_line = np.linspace(duration_min, duration_max, 100)
                    y_line = slope * x_line + intercept
                    ax.plot(x_line, y_line, color=subtype_colors[subtype], linewidth=1.5, alpha=0.7)
    
    ax.set_title(aa_name, fontsize=10, fontweight='bold')
    ax.set_xticks([0, 20, 40])
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.grid(True, alpha=0.3)

# Use last subplot for legend
axes[-1].axis('off')
legend_handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=subtype_colors[st], 
                              markersize=10, label=st) for st in subtype_order]
axes[-1].legend(handles=legend_handles, loc='center', fontsize=12, title='Disease Subtype', 
                title_fontsize=12, frameon=True, fancybox=True, shadow=True)

fig.suptitle('Amino Acid Concentrations vs Disease Duration by Subtype', fontsize=14, fontweight='bold', y=0.98)
fig.text(0.02, 0.5, 'Concentration [nmol/ml]', va='center', rotation='vertical', fontsize=12, fontweight='bold')
fig.text(0.5, 0.02, 'Duration [years]', ha='center', fontsize=12, fontweight='bold')

plt.tight_layout(rect=[0.04, 0.04, 1, 0.96])
plt.savefig(os.path.join(OUTPUT_DIR, 'Fig9_Duration_Grid.png'), dpi=150, bbox_inches='tight')
plt.close()


# ==========================================
# GENERATE JSON DATA (Node 4 Style)
# ==========================================

def get_regression_stats(x, y):
    """Calculates linear regression stats for web visualization"""
    if len(x) < 2:
        return None
    slope, intercept, r_value, p_value, std_err = linregress(x, y)
    # Consistent scale for lines
    line_x = [float(min(x)), float(max(x))]
    line_y = [slope * xi + intercept for xi in line_x]
    
    return {
        'slope': float(slope),
        'intercept': float(intercept),
        'r_squared': float(r_value**2),
        'p_value': float(p_value),
        'line_x': line_x,
        'line_y': line_y
    }

json_data = {
    'metadata': {
        'title': 'Subtype Trajectories',
        'amino_acids': aa_clean_names
    },
    'fig9': {
        'type': 'scatter', 
        'variable': 'Duration',
        'traces': []
    }
}

for aa_col, aa_name in zip(aa_cols, aa_clean_names):
    trace_item = {'aa': aa_name}
    
    for subtype in subtype_order:
        subtype_data = df_fig9[df_fig9['Subtype'] == subtype]
        if len(subtype_data) > 0:
            valid_data = subtype_data[['Duration', aa_col]].dropna()
            
            # Scatter Data
            x_vals = valid_data['Duration'].tolist()
            y_vals = valid_data[aa_col].tolist()
            
            # Stats
            stats_obj = get_regression_stats(valid_data['Duration'].values, valid_data[aa_col].values)
            
            trace_item[subtype] = {
                'x': x_vals,
                'y': y_vals,
                'stats': stats_obj
            }
            
    json_data['fig9']['traces'].append(trace_item)


json_path = os.path.join(OUTPUT_DIR, 'trajectories_data.json')
with io.open(json_path, 'w', encoding='utf-8') as f:
    f.write(json.dumps(json_data, ensure_ascii=False))

# Generate HTML
html_content = generate_trajectories_html(json_filename='trajectories_data.json')

html_path = os.path.join(OUTPUT_DIR, 'trajectories.html')
with io.open(html_path, 'w', encoding='utf-8') as f:
    f.write(html_content)

print("Node 9 completed.")
