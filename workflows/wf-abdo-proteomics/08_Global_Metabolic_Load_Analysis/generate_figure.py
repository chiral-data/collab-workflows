"""Node 8: Global Metabolic Load Analysis - Generate Fig 8"""
import os
import sys
import json
import io
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

from html_generator import generate_metabolic_load_html

# Configuration
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
INPUT_FILE = "data_standardized.pkl" 
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

sns.set_style("whitegrid")
plt.rcParams['font.family'] = 'sans-serif'

print(">>> NODE 8: GLOBAL METABOLIC LOAD ANALYSIS...")

# Load data
try:
    df = pd.read_pickle(INPUT_FILE)
except Exception as e:
    print("Warning: Input file not found. Ensure data_standardized.pkl exists.")
    sys.exit(1)

# Define masks
ms_types = ['RRMS', 'SPMS', 'PPMS']
mg_types = ['general', 'eye-type']
masks = {
    'MS': df['Type'].isin(ms_types),
    'MG': df['Type'].isin(mg_types)
}

df_mg_ms = df[masks['MS'] | masks['MG']].copy()

# ==========================================
# FIG 8: Total AA by Type
# ==========================================
print("Generating Fig 8: Total AA by Type...")

df_mg_ms['Plot_Type'] = np.where(df_mg_ms['Type'].isin(['general', 'eye-type']), 'GMG', df_mg_ms['Type'])

# Static Plot (Legacy Support)
plt.figure(figsize=(8,6))
type_order = ['PPMS', 'SPMS', 'RRMS', 'GMG']
type_colors = {'PPMS': '#1f77b4', 'SPMS': '#ff7f0e', 'RRMS': '#2ca02c', 'GMG': '#d62728'}

sns.boxplot(data=df_mg_ms, x='Plot_Type', y='Total_AA', order=type_order,
            palette=type_colors, showfliers=True, width=0.5, linewidth=1.5,
            flierprops={'marker': 'o', 'markerfacecolor': 'white', 'markeredgecolor': 'black', 'markersize': 5})

plt.ylabel('Concentration [nmol/ml]')
plt.xlabel('')
plt.yticks([-4, -2, 0, 2, 4, 6, 8, 10, 12])
plt.title('Total Amino Acid Concentration by Disease Type', fontsize=12, fontweight='bold')

plt.savefig(os.path.join(OUTPUT_DIR, 'Fig8_TotalAA_Type.png'))
plt.close()
print("Saved: Fig8_TotalAA_Type.png")


# ==========================================
# STATISTICAL ANALYSIS (Kruskal-Wallis)
# ==========================================
print("Calculating Statistics...")
groups = []
group_names = []
for t in type_order:
    vals = df_mg_ms[df_mg_ms['Plot_Type'] == t]['Total_AA'].dropna().values
    if len(vals) > 0:
        groups.append(vals)
        group_names.append(t)

p_value = 1.0
test_name = "N/A"
if len(groups) > 1:
    try:
        stat, p_value = stats.kruskal(*groups)
        test_name = "Kruskal-Wallis"
        print("Kruskal-Wallis p-value: " + str(p_value))
    except Exception as e:
        print("Stats error: " + str(e))

# ==========================================
# GENERATE JSON/HTML DATA
# ==========================================
print("Generating JSON and HTML...")

traces = []
for t in type_order:
    y_vals = df_mg_ms[df_mg_ms['Plot_Type'] == t]['Total_AA'].dropna().tolist()
    # Adding jittered x-values for strip plot could be done in JS, 
    # but sending raw y-data allows JS to handle Box/Violin logic perfectly.
    traces.append({
        'name': t,
        'y': y_vals,
        'color': type_colors.get(t, '#333')
    })

json_data = {
    'metadata': {'title': 'Global Metabolic Load'},
    'stats': {
        'test': test_name,
        'p_value': p_value,
        'p_value_fmt': "< 0.001" if p_value < 0.001 else "{:.4f}".format(p_value)
    },
    'fig8': {
        'title': 'Total Amino Acid Concentration by Disease Type',
        'yaxis': 'Concentration [nmol/ml]',
        'traces': traces
    }
}

json_path = os.path.join(OUTPUT_DIR, 'metabolic_load_data.json')
with io.open(json_path, 'w', encoding='utf-8') as f:
    f.write(json.dumps(json_data, ensure_ascii=False))
print("Saved: " + json_path)

html_content = generate_metabolic_load_html('metabolic_load_data.json')

html_path = os.path.join(OUTPUT_DIR, 'metabolic_load.html')
with io.open(html_path, 'w', encoding='utf-8') as f:
    f.write(html_content)
print("Saved: " + html_path)

print("Node 8 completed successfully.")
