"""Node 6: Differential Diagnosis Biomarkers - Generate Fig 4 & 5 (Python 3)"""
import os
import sys
import json
import traceback
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu
from matplotlib.patches import Patch

# Import local HTML generator
try:
    from html_generator import generate_biomarkers_html
except ImportError as e:
    print(f"Error importing html_generator: {e}")
    sys.exit(1)

def main():
    try:
        # Configuration
        SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
        INPUT_FILE = "inputs/data_standardized.pkl"
        OUTPUT_DIR = "outputs"

        if not os.path.exists(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR, exist_ok=True)
            print(f"Created output directory: {OUTPUT_DIR}")

        sns.set_style("whitegrid")
        plt.rcParams['font.family'] = 'sans-serif'

        print(">>> NODE 6: DIFFERENTIAL DIAGNOSIS BIOMARKERS...")

        # Load data
        if not os.path.exists(INPUT_FILE):
             print(f"ERROR: Input file {INPUT_FILE} not found in {os.getcwd()}")
             sys.exit(1)

        df = pd.read_pickle(INPUT_FILE)
        print(f"Loaded {len(df)} records")

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
        print("Generating Fig 4: Specific Differences...")

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
            
            # Add p-value (Static Image)
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
        out_path_4 = os.path.join(OUTPUT_DIR, 'Fig4_Specific_Diffs.png')
        plt.savefig(out_path_4)
        plt.close()
        print(f"Saved: {out_path_4}")

        # ==========================================
        # FIG 5: Females Only
        # ==========================================
        print("Generating Fig 5: Female Specific...")

        # Note: Using 'Female' as per existing data convention
        df_fem = df_mg_ms[df_mg_ms['Sex'] == 'Female'].copy()
        fig, axes = plt.subplots(1, 3, figsize=(15, 5))

        for i, aa in enumerate(amino_acids_fig4):
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
        out_path_5 = os.path.join(OUTPUT_DIR, 'Fig5_Female_Specific.png')
        plt.savefig(out_path_5)
        plt.close()
        print(f"Saved: {out_path_5}")


        # ==========================================
        # GENERATE JSON DATA
        # ==========================================
        print("Generating JSON data...")

        json_data = {
            'metadata': {
                'title': 'Differential Diagnosis Biomarkers'
            },
            'fig4': {
                'title': 'General Cohort (MS vs MG)',
                'variable': 'Group',
                'traces': []
            },
            'fig5': {
                'title': 'Female Subgroup (MS vs MG)',
                'variable': 'Group',
                'traces': []
            }
        }

        # Helper to process cohorts
        def process_cohort(cohort_df, key):
            for aa in amino_acids_fig4:
                aa_clean = aa.replace('_conc', '')
                ms_vals = cohort_df[cohort_df['Group']=='MS'][aa].dropna()
                mg_vals = cohort_df[cohort_df['Group']=='MG'][aa].dropna()
                
                # Calc Stats
                p_val = 1.0
                if len(ms_vals) > 0 and len(mg_vals) > 0:
                    try:
                        _, p_val = mannwhitneyu(ms_vals, mg_vals, alternative='two-sided')
                    except Exception:
                        p_val = 1.0

                json_data[key]['traces'].append({
                    'aa': aa_clean,
                    'MS': {'y': ms_vals.tolist()},
                    'MG': {'y': mg_vals.tolist()},
                    'stats': {'p_value': float(p_val)}
                })

        # Process Fig 4 (General)
        process_cohort(df_mg_ms, 'fig4')
        
        # Process Fig 5 (Female)
        process_cohort(df_fem, 'fig5')

        json_path = os.path.join(OUTPUT_DIR, 'biomarkers_data.json')
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(json_data, f, indent=2, ensure_ascii=False)
        print(f"Saved JSON: {json_path}")

        # Generate HTML
        html_content = generate_biomarkers_html(json_filename="biomarkers_data.json")

        html_path = os.path.join(OUTPUT_DIR, 'biomarkers.html')
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"Saved HTML: {html_path}")

        print("Node 6 completed successfully.")

    except Exception:
        print("CRITICAL ERROR IN GENERATE_FIGURES.PY:")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
