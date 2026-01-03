"""Node 4: Clinical Confounders & Validation - Generate Enriched JSON"""
import os
import sys
import json
import pandas as pd
import numpy as np
from scipy import stats  # Added for regression calculations

# Configuration
INPUT_FILE = "data_standardized.pkl"
AA_COLS_FILE = "aa_cols.txt"
OUTPUT_DIR = "outputs"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

print(">>> NODE 4: CLINICAL CONFOUNDERS & VALIDATION (Enriched)...", flush=True)

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

# Initialize JSON structure
json_data = {
    'metadata': {
        'title': 'Clinical Confounders & Validation',
        'amino_acids': aa_cols
    },
    'fig2a': {'variable': 'Age', 'traces': []},
    'fig2b': {'variable': 'Disease Duration', 'traces': []},
    'fig2c': {'variable': 'EDSS', 'traces': []}
}

def get_regression_stats(x, y):
    """Calculates linear regression stats for web visualization"""
    if len(x) < 2:
        return None
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
    # Create line coordinates for plotting
    line_x = [min(x), max(x)]
    line_y = [slope * xi + intercept for xi in line_x]
    
    return {
        'slope': slope,
        'intercept': intercept,
        'r_squared': r_value**2,
        'p_value': p_value,
        'line_x': line_x,
        'line_y': line_y
    }

for aa in aa_cols:
    # Use the clean name (first part before space if exists)
    aa_clean = aa.split(' ')[0]
    
    # --- Fig 2A: Age (MS vs Control) ---
    ms_data = df_p1a[df_p1a['Group']=='MS'][[aa, 'Age']].dropna()
    ctrl_data = df_p1a[df_p1a['Group']=='Control'][[aa, 'Age']].dropna()
    
    # Calculate Stats
    ms_stats = get_regression_stats(ms_data['Age'], ms_data[aa])
    ctrl_stats = get_regression_stats(ctrl_data['Age'], ctrl_data[aa])
    
    json_data['fig2a']['traces'].append({
        'aa': aa_clean,
        'MS': {
            'x': ms_data['Age'].tolist(), 
            'y': ms_data[aa].tolist(),
            'stats': ms_stats
        },
        'Control': {
            'x': ctrl_data['Age'].tolist(), 
            'y': ctrl_data[aa].tolist(),
            'stats': ctrl_stats
        }
    })
    
    # --- Fig 2B: Duration (MS only) ---
    dur_data = df_p1a[df_p1a['Group']=='MS'][[aa, 'Duration']].dropna()
    dur_stats = get_regression_stats(dur_data['Duration'], dur_data[aa])
    
    json_data['fig2b']['traces'].append({
        'aa': aa_clean,
        'MS': {
            'x': dur_data['Duration'].tolist(), 
            'y': dur_data[aa].tolist(),
            'stats': dur_stats
        }
    })
    
    # --- Fig 2C: EDSS (MS only) ---
    edss_data = df_p1a[df_p1a['Group']=='MS'][[aa, 'EDSS']].dropna()
    edss_stats = get_regression_stats(edss_data['EDSS'], edss_data[aa])
    
    json_data['fig2c']['traces'].append({
        'aa': aa_clean,
        'MS': {
            'x': edss_data['EDSS'].tolist(), 
            'y': edss_data[aa].tolist(),
            'stats': edss_stats
        }
    })

# Save JSON
json_path = os.path.join(OUTPUT_DIR, 'confounders_data.json')
with open(json_path, 'w', encoding='utf-8') as f:
    json.dump(json_data, f, indent=2)

print(f"Success! Enriched JSON saved to {json_path}")
print("Run html_generator.py next.")
