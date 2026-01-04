"""Node 10: Clinical Mimicry Test - Generate Fig 10 & 11 (Python 3)"""
import os
import sys
import json
import traceback
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import mannwhitneyu

# Import local HTML generator
try:
    from html_generator import generate_mimicry_html
except ImportError as e:
    print(f"Error importing html_generator: {e}")
    sys.exit(1)

def main():
    try:
        # Configuration
        SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
        INPUT_FILE = "data_standardized.pkl"
        AA_COLS_FILE = "aa_cols.txt"
        OUTPUT_DIR = "outputs"

        if not os.path.exists(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR, exist_ok=True)
            print(f"Created output directory: {OUTPUT_DIR}")

        sns.set_style("whitegrid")
        plt.rcParams['font.family'] = 'sans-serif'

        print(">>> NODE 10: CLINICAL MIMICRY TEST (RRMS vs MG)...")

        # Load data
        if not os.path.exists(INPUT_FILE):
             print(f"ERROR: Input file {INPUT_FILE} not found in {os.getcwd()}")
             sys.exit(1)

        df = pd.read_pickle(INPUT_FILE)
        
        if not os.path.exists(AA_COLS_FILE):
             print(f"ERROR: Column file {AA_COLS_FILE} not found in {os.getcwd()}")
             sys.exit(1)
             
        with open(AA_COLS_FILE, 'r', encoding='utf-8') as f:
            aa_cols = [line.strip() for line in f.readlines()]
        
        print(f"Loaded {len(df)} records")

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

        # Calculate statistics for Total AA
        rrms_total = df_rm[df_rm['Group']=='RRMS']['Total_AA'].dropna()
        mg_total = df_rm[df_rm['Group']=='MG']['Total_AA'].dropna()
        
        p_val_total = 1.0
        if len(rrms_total) > 0 and len(mg_total) > 0:
            try:
                _, p_val_total = mannwhitneyu(rrms_total, mg_total, alternative='two-sided')
            except Exception:
                p_val_total = 1.0

        fig, ax = plt.subplots(figsize=(6, 6))
        # Removed redundant call without hue
        sns.boxplot(data=df_rm, x='Group', y='Total_AA', hue='Group', palette={'RRMS':'grey', 'MG':'skyblue'}, ax=ax, legend=False)
        ax.set_ylabel('Concentration [nmol/ml]', fontsize=12)
        ax.set_xlabel('', fontsize=12)
        ax.set_title('Total Amino Acid Concentration: MS-RRMS vs MG', fontsize=12, fontweight='bold')
        
        # Annotation
        y_ticks = [-5, -2.5, 0, 2.5, 5, 7.5, 10, 12.5]
        ax.set_yticks(y_ticks)
        y_max = df_rm['Total_AA'].max()
        annotation_y_pos = max(y_max, y_ticks[-1]) * 1.05
        ax.annotate(f'p = {p_val_total:.1e}', xy=(0.5, annotation_y_pos), fontsize=11, ha='center', fontweight='bold')
        ax.set_ylim(y_ticks[0], annotation_y_pos * 1.1)

        plt.tight_layout()
        out_path_10 = os.path.join(OUTPUT_DIR, 'Fig10_TotalAA_RRMS_MG.png')
        plt.savefig(out_path_10, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved: {out_path_10}")

        # ==========================================
        # FIG 11: Grid (Static)
        # ==========================================
        print("Generating Fig 11: Static Grid...")

        df_melt_11 = df_rm.melt(id_vars='Group', value_vars=aa_cols, var_name='Amino Acid', value_name='Concentration')
        df_melt_11['Amino Acid'] = df_melt_11['Amino Acid'].str.replace('_conc', '', regex=False)

        g = sns.FacetGrid(df_melt_11, col='Amino Acid', col_wrap=6, sharey=False, height=2.5, aspect=1)

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
            
        out_path_11 = os.path.join(OUTPUT_DIR, 'Fig11_RRMS_vs_MG_Grid.png')
        g.savefig(out_path_11)
        plt.close()
        print(f"Saved: {out_path_11}")


        # ==========================================
        # GENERATE JSON DATA
        # ==========================================
        print("Generating JSON data...")

        json_data = {
            'metadata': {
                'title': 'Clinical Mimicry Test: RRMS vs MG',
                'amino_acids': [aa.replace('_conc', '') for aa in aa_cols]
            },
            'fig10': {
                'title': 'Total Amino Acid Overview',
                'variable': 'Group',
                'traces': [
                    {
                        'name': 'Total AA',
                        'RRMS': {'y': rrms_total.tolist()},
                        'MG': {'y': mg_total.tolist()},
                        'stats': {'p_value': float(p_val_total)}
                    }
                ]
            },
            'fig11': {
                'title': 'Detailed Amino Acid Analysis',
                'variable': 'Amino Acid', 
                'traces': []
            }
        }

        # Individual Amino Acids (Fig 11)
        for aa in aa_cols:
            aa_clean = aa.replace('_conc', '')
            rrms_vals = df_rm[df_rm['Group']=='RRMS'][aa].dropna()
            mg_vals = df_rm[df_rm['Group']=='MG'][aa].dropna()
            
            # Calc Stats (Mann-Whitney U)
            p_val = 1.0
            if len(rrms_vals) > 0 and len(mg_vals) > 0:
                try:
                    _, p_val = mannwhitneyu(rrms_vals, mg_vals, alternative='two-sided')
                except Exception:
                    p_val = 1.0

            json_data['fig11']['traces'].append({
                'aa': aa_clean,
                'RRMS': {'y': rrms_vals.tolist()},
                'MG': {'y': mg_vals.tolist()},
                'stats': {'p_value': float(p_val)}
            })

        json_path = os.path.join(OUTPUT_DIR, 'mimicry_data.json')
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(json_data, f, indent=2, ensure_ascii=False)
        print(f"Saved JSON: {json_path}")

        # Generate HTML
        html_content = generate_mimicry_html(json_filename="mimicry_data.json")

        html_path = os.path.join(OUTPUT_DIR, 'clinical_mimicry.html')
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"Saved HTML: {html_path}")

        print("Node 10 completed successfully.")

    except Exception:
        print("CRITICAL ERROR IN GENERATE_FIGURES.PY:")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
