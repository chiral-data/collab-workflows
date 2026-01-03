"""Node 2: Cohort Demographics - Generate Tables 1 & 2"""
import os
import sys
import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Import local HTML generator
from html_generator import generate_demographics_html

# Configuration
INPUT_FILE = "data_standardized.pkl"
OUTPUT_DIR = "outputs"

if not os.path.exists(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

print(">>> NODE 2: COHORT DEMOGRAPHICS...", flush=True)

# Load data
df = pd.read_pickle(INPUT_FILE)
print(f"Loaded {len(df)} records", flush=True)

# Define masks
ms_types = ['RRMS', 'SPMS', 'PPMS']
mg_types = ['general', 'eye-type']

masks = {
    'MS': df['Type'].isin(ms_types),
    'MG': df['Type'].isin(mg_types),
    'Control': df['Status'] == 'control',
    'RRMS': df['Type'] == 'RRMS',
    'SPMS': df['Type'] == 'SPMS',
    'PPMS': df['Type'] == 'PPMS'
}

# ==========================================
# TABLE 1: MS vs Controls
# ==========================================
print("Generating Table 1...", flush=True)

# Calculate MS statistics
ms_data = df[masks['MS']]
n_ms = len(ms_data)
ms_males = len(ms_data[ms_data['Sex'] == 'Male'])
ms_females = len(ms_data[ms_data['Sex'] == 'Female'])
ms_age_mean = ms_data['Age'].mean()
ms_age_sd = ms_data['Age'].std()
ms_duration_mean = ms_data['Duration'].mean()
ms_duration_sd = ms_data['Duration'].std()
ms_duration_min = ms_data['Duration'].min()
ms_duration_max = ms_data['Duration'].max()
ms_edss_median = ms_data['EDSS'].median()
ms_edss_q1 = ms_data['EDSS'].quantile(0.25)
ms_edss_q3 = ms_data['EDSS'].quantile(0.75)

# Calculate Control statistics
ctrl_data = df[masks['Control']]
n_ctrl = len(ctrl_data)
ctrl_males = len(ctrl_data[ctrl_data['Sex'] == 'Male'])
ctrl_females = len(ctrl_data[ctrl_data['Sex'] == 'Female'])
ctrl_age_mean = ctrl_data['Age'].mean()
ctrl_age_sd = ctrl_data['Age'].std()

# Build table data
table_data = [
    ['Number of subjects (n)', str(n_ms), str(n_ctrl)],
    ['Sex (male/female)', f'{ms_males}/{ms_females}', f'{ctrl_males}/{ctrl_females}'],
    ['Age, years', f'{ms_age_mean:.1f} ± {ms_age_sd:.1f}', f'{ctrl_age_mean:.1f} ± {ctrl_age_sd:.1f}'],
    ['Disease duration (years) [range]', f'{ms_duration_mean:.1f} ± {ms_duration_sd:.1f} [{ms_duration_min:.0f}–{ms_duration_max:.0f}]', 'n.a.'],
    ['Median EDSS score (IQR)', f'{ms_edss_median:.1f} ({ms_edss_q1:.1f}–{ms_edss_q3:.1f})', 'n.a.'],
    ['MS type n (%)', '', ''],
    ['RRMS', f'{len(df[masks["RRMS"]])} ({len(df[masks["RRMS"]])/n_ms*100:.0f}%)', ''],
    ['SPMS', f'{len(df[masks["SPMS"]])} ({len(df[masks["SPMS"]])/n_ms*100:.0f}%)', ''],
    ['PPMS', f'{len(df[masks["PPMS"]])} ({len(df[masks["PPMS"]])/n_ms*100:.0f}%)', ''],
    ['Type of DMT n (%)', '', ''],
]

# Add DMT rows
dmt_counts = df[masks['MS']]['DMT_Clean'].value_counts()
for drug, count in dmt_counts.items():
    if pd.notna(drug):
        table_data.append([drug, f'{count} ({count/n_ms*100:.0f}%)', ''])

# Create DataFrame
df_table = pd.DataFrame(table_data, columns=['', 'MS', 'HCs'])

# Render as Image
fig, ax = plt.subplots(figsize=(10, len(table_data) * 0.5 + 1.5))
ax.axis('tight')
ax.axis('off')

the_table = ax.table(cellText=df_table.values,
                     colLabels=df_table.columns,
                     cellLoc='left',
                     loc='center',
                     colWidths=[0.5, 0.25, 0.25])

the_table.auto_set_font_size(False)
the_table.set_fontsize(10)
the_table.scale(1, 2)

# Format cells
for (i, j), cell in the_table.get_celld().items():
    if i == 0:
        cell.set_facecolor('#40466e')
        cell.set_text_props(weight='bold', color='white')
    else:
        cell.set_facecolor('#f1f1f2' if i % 2 == 0 else 'white')
    cell.set_edgecolor('black')
    cell.set_linewidth(0.5)

plt.title("Table 1. Demographic and clinical data of study participants", 
          y=0.98, fontsize=12, weight='bold', loc='left', pad=20)

save_path = os.path.join(OUTPUT_DIR, 'Table1_Demographics.png')
plt.savefig(save_path, bbox_inches='tight', dpi=300, facecolor='white')
plt.close()
print(f"Saved: {save_path}", flush=True)

# ==========================================
# TABLE 2: MS vs MG
# ==========================================
print("Generating Table 2...", flush=True)

# Calculate MG statistics
mg_data = df[masks['MG']]
n_mg = len(mg_data)
mg_males = len(mg_data[mg_data['Sex'] == 'Male'])
mg_females = len(mg_data[mg_data['Sex'] == 'Female'])
mg_age_mean = mg_data['Age'].mean()
mg_age_sd = mg_data['Age'].std()
mg_age_min = mg_data['Age'].min()
mg_age_max = mg_data['Age'].max()
mg_duration_mean = mg_data['Duration'].mean()
mg_duration_sd = mg_data['Duration'].std()
mg_duration_min = mg_data['Duration'].min()
mg_duration_max = mg_data['Duration'].max()

# Count MG subtypes
n_gmg = len(df[df['Type'] == 'general'])
n_omg = len(df[df['Type'] == 'eye-type'])

# Build table data
table2_data = [
    ['Number of subjects (n)', str(n_ms), str(n_mg)],
    ['Sex (male/female)', f'{ms_males}/{ms_females} ({ms_males/n_ms*100:.1f}%/{ms_females/n_ms*100:.1f}%)', 
     f'{mg_males}/{mg_females} ({mg_males/n_mg*100:.1f}%/{mg_females/n_mg*100:.1f}%)'],
    ['Age, years (mean ± SD)\nMin–Max', f'{ms_age_mean:.1f} ± {ms_age_sd:.2f}\n{ms_data["Age"].min():.0f}–{ms_data["Age"].max():.0f}', 
     f'{mg_age_mean:.2f} ± {mg_age_sd:.2f}\n{mg_age_min:.0f}–{mg_age_max:.0f}'],
    ['Disease duration (years)\n[range]', f'{ms_duration_mean:.1f} ± {ms_duration_sd:.1f}\n[{ms_duration_min:.0f}–{ms_duration_max:.0f}]', 
     f'{mg_duration_mean:.2f} ± {mg_duration_sd:.2f}\n[{mg_duration_min:.0f}–{mg_duration_max:.0f}]'],
    ['Median EDSS score (IQR)', f'{ms_edss_median:.1f} ({ms_edss_q1:.1f}–{ms_edss_q3:.1f})', 'n.a.'],
    ['MS type n (%)', '', ''],
    ['RRMS', f'{len(df[masks["RRMS"]])} ({len(df[masks["RRMS"]])/n_ms*100:.0f}%)', ''],
    ['SPMS', f'{len(df[masks["SPMS"]])} ({len(df[masks["SPMS"]])/n_ms*100:.0f}%)', ''],
    ['PPMS', f'{len(df[masks["PPMS"]])} ({len(df[masks["PPMS"]])/n_ms*100:.0f}%)', ''],
    ['MG type n (%)', '', ''],
    ['GMG', '', f'{n_gmg} ({n_gmg/n_mg*100:.1f}%)'],
    ['OMG', '', f'{n_omg} ({n_omg/n_mg*100:.1f}%)'],
]

df_table2 = pd.DataFrame(table2_data, columns=['', 'MS', 'MG'])

# Render Table 2
fig, ax = plt.subplots(figsize=(10, len(table2_data) * 0.5 + 2))
ax.axis('tight')
ax.axis('off')

the_table2 = ax.table(cellText=df_table2.values,
                      colLabels=df_table2.columns,
                      cellLoc='left',
                      loc='center',
                      colWidths=[0.5, 0.25, 0.25])

the_table2.auto_set_font_size(False)
the_table2.set_fontsize(9)
the_table2.scale(1, 2.2)

# Format cells
for (i, j), cell in the_table2.get_celld().items():
    if i == 0:
        cell.set_facecolor('#40466e')
        cell.set_text_props(weight='bold', color='white')
    else:
        cell.set_facecolor('#f1f1f2' if i % 2 == 0 else 'white')
    cell.set_edgecolor('black')
    cell.set_linewidth(0.5)

plt.title("Table 2. Demographic and clinical data of patients", 
          y=0.98, fontsize=12, weight='bold', loc='left', pad=20)

save_path2 = os.path.join(OUTPUT_DIR, 'Table2_MS_MG_Demographics.png')
plt.savefig(save_path2, bbox_inches='tight', dpi=300, facecolor='white')
plt.close()
print(f"Saved: {save_path2}", flush=True)

# ==========================================
# GENERATE JSON DATA
# ==========================================
print("Generating JSON data...", flush=True)

json_data = {
    'metadata': {
        'title': 'Cohort Demographics',
        'total_samples': len(df),
        'ms_count': n_ms,
        'mg_count': n_mg,
        'control_count': n_ctrl
    },
    'table1': {
        'columns': df_table.columns.tolist(),
        'rows': df_table.values.tolist()
    },
    'table2': {
        'columns': df_table2.columns.tolist(),
        'rows': df_table2.values.tolist()
    }
}

json_path = os.path.join(OUTPUT_DIR, 'demographics_data.json')
with open(json_path, 'w', encoding='utf-8') as f:
    json.dump(json_data, f, indent=2, ensure_ascii=False)
print(f"Saved: {json_path}", flush=True)

# ==========================================
# GENERATE HTML VISUALIZATION
# ==========================================
print("Generating HTML visualization...", flush=True)

html_content = generate_demographics_html(json_filename="demographics_data.json")

html_path = os.path.join(OUTPUT_DIR, 'demographics.html')
with open(html_path, 'w', encoding='utf-8') as f:
    f.write(html_content)
print(f"Saved: {html_path}", flush=True)

print("Node 2 completed successfully.", flush=True)

