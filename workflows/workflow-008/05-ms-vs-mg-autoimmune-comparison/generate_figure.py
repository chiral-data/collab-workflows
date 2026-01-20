"""Node 5: MS vs MG Autoimmune Comparison - Generate Enriched JSON (Python 3)"""
import os
import sys
import json
import traceback
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Import local HTML generator
try:
    from html_generator import generate_autoimmune_html
except ImportError as e:
    print(f"Error importing html_generator: {e}")
    sys.exit(1)

def main():
    try:
        # Configuration
        SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
        INPUT_FILE = "inputs/data_standardized.pkl"
        AA_COLS_FILE = "inputs/aa_cols.txt"
        OUTPUT_DIR = "outputs"

        if not os.path.exists(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR, exist_ok=True)
            print(f"Created output directory: {OUTPUT_DIR}")

        # Set style
        sns.set_style("whitegrid")
        plt.rcParams['font.family'] = 'sans-serif'

        print(">>> NODE 5: MS VS MG AUTOIMMUNE COMPARISON...")

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
        ms_types = ['RRMS', 'SPMS', 'PPMS']
        mg_types = ['general', 'eye-type']
        masks = {
            'MS': df['Type'].isin(ms_types),
            'MG': df['Type'].isin(mg_types),
            'Control': df['Status'] == 'control'
        }

        # Create combined dataframe
        df_p2 = df.copy()
        df_p2['Group'] = 'Other'
        df_p2.loc[masks['MS'] | masks['MG'], 'Group'] = 'MS+MG'
        df_p2.loc[masks['Control'], 'Group'] = 'Controls'

        # Filter for relevant groups
        df_plot = df_p2[df_p2['Group'].isin(['MS+MG', 'Controls'])].copy()

        # ==========================================
        # HERO PLOT: Total AA Comparison
        # ==========================================
        print("Generating Hero Plot Data (Total AA)...")
        
        # Calc Stats for Total AA
        ms_mg_total = df_plot[df_plot['Group']=='MS+MG']['Total_AA'].dropna()
        ctrl_total = df_plot[df_plot['Group']=='Controls']['Total_AA'].dropna()
        
        p_val_total = 1.0
        if len(ms_mg_total) > 1 and len(ctrl_total) > 1:
            try:
                # Welch's t-test
                _, p_val_total = stats.ttest_ind(ms_mg_total, ctrl_total, equal_var=False)
            except Exception:
                p_val_total = 1.0

        # Static Plot for Fig 3 (Legacy support)
        print("Generating Static Fig 3...")
        df_melt = df_plot.melt(id_vars='Group', value_vars=aa_cols, var_name='Amino Acid', value_name='LogC')
        df_melt['Amino Acid'] = df_melt['Amino Acid'].str.replace('_conc', '', regex=False)

        g = sns.FacetGrid(df_melt, col='Amino Acid', col_wrap=6, sharey=False, height=2.5, aspect=1)
        g.map_dataframe(sns.boxplot, x='Group', y='LogC', 
                        hue='Group', legend=False,
                        palette={'MS+MG':'skyblue', 'Controls':'grey'}, 
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
        
        out_path_3 = os.path.join(OUTPUT_DIR, 'Fig3_MS_MG_vs_Control.png')
        g.savefig(out_path_3)
        plt.close()
        print(f"Saved: {out_path_3}")

        # ==========================================
        # GENERATE JSON DATA
        # ==========================================
        print("Generating JSON data...")

        json_data = {
            'metadata': {
                'title': 'MS+MG vs Controls',
                'description': 'Comparison of autoimmune cohort (MS and MG combined) versus healthy controls.'
            },
            'hero': {
                'title': 'Total Amino Acid Overview',
                'variable': 'Group',
                'traces': [
                    {
                        'name': 'Total AA',
                        'MS+MG': {'y': ms_mg_total.tolist()},
                        'Controls': {'y': ctrl_total.tolist()},
                        'stats': {'p_value': float(p_val_total)}
                    }
                ]
            },
            'grid': {
                'title': 'Individual Amino Acids',
                'variable': 'Amino Acid',
                'traces': []
            }
        }

        # Individual Amino Acids (Grid)
        for aa in aa_cols:
            aa_clean = aa.replace('_conc', '')
            
            ms_mg_data = df_plot[df_plot['Group']=='MS+MG'][aa].dropna()
            ctrl_data = df_plot[df_plot['Group']=='Controls'][aa].dropna()
            
            # Calculate Statistics (Welch's T-test)
            p_val = 1.0
            if len(ms_mg_data) > 1 and len(ctrl_data) > 1:
                try:
                    _, p_val = stats.ttest_ind(ms_mg_data, ctrl_data, equal_var=False)
                except Exception:
                    p_val = 1.0
            
            json_data['grid']['traces'].append({
                'aa': aa_clean,
                'MS+MG': {'y': ms_mg_data.tolist()},
                'Controls': {'y': ctrl_data.tolist()},
                'stats': {'p_value': float(p_val)}
            })

        json_path = os.path.join(OUTPUT_DIR, 'autoimmune_data.json')
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(json_data, f, indent=2, ensure_ascii=False)
        print(f"Saved JSON: {json_path}")

        # Generate HTML
        html_content = generate_autoimmune_html(json_filename='autoimmune_data.json')
        
        html_path = os.path.join(OUTPUT_DIR, 'autoimmune.html')
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"Saved HTML: {html_path}")

        print("Node 5 completed successfully.")

    except Exception:
        print("CRITICAL ERROR IN GENERATE_FIGURE.PY:")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
