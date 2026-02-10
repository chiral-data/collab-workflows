"""Node 3: MS Pathology Overview - Generate Fig 1A & 1B (Python 3)"""
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
    from html_generator import generate_pathology_html
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

        sns.set_style("whitegrid")
        plt.rcParams['font.family'] = 'sans-serif'

        print(">>> NODE 3: MS PATHOLOGY OVERVIEW...")

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
        masks = {
            'MS': df['Type'].isin(ms_types),
            'Control': df['Status'] == 'control'
        }

        # ==========================================
        # FIG 1A: MS vs Control (Static + Data)
        # ==========================================
        print("Generating Fig 1A: MS vs Control...")

        df_p1a = df[masks['MS'] | masks['Control']].copy()
        df_p1a['Group'] = np.where(df_p1a['Status']=='control', 'Control', 'MS')

        # Static Plot
        print("Generating Static Fig 1A...")
        df_melt = df_p1a.melt(id_vars='Group', value_vars=aa_cols, var_name='Amino Acid', value_name='LogC')
        df_melt['Amino Acid'] = df_melt['Amino Acid'].str.replace('_conc', '', regex=False)

        g = sns.FacetGrid(df_melt, col='Amino Acid', col_wrap=6, sharey=False, height=2.5, aspect=1)
        g.map_dataframe(sns.boxplot, x='Group', y='LogC', 
                        hue='Group', palette={'MS':'steelblue', 'Control':'salmon'}, legend=False,
                        showfliers=True, 
                        flierprops={'marker': 'o', 'markerfacecolor': 'black', 'markersize': 3})

        g.set_titles("{col_name}")
        g.set_xlabels("")

        y_ticks = [0, 10]
        for i, ax in enumerate(g.axes.flat):
            for spine in ax.spines.values():
                spine.set_edgecolor('black')
                spine.set_linewidth(1.5)
            ax.set_yticks(y_ticks)
            if i % 6 != 0:
                ax.set_yticklabels([])
            ax.set_ylabel('')

        out_path_1a = os.path.join(OUTPUT_DIR, 'Fig1A_MS_vs_Control.png')
        g.savefig(out_path_1a)
        plt.close()
        print(f"Saved: {out_path_1a}")

        # ==========================================
        # FIG 1B: MS Subtypes (Static + Data)
        # ==========================================
        print("Generating Fig 1B: MS Subtypes...")

        df_p1b = df[masks['MS']].copy()
        
        # Static Plot
        print("Generating Static Fig 1B...")
        df_melt_b = df_p1b.melt(id_vars='Type', value_vars=aa_cols, var_name='Amino Acid', value_name='LogC')
        df_melt_b['Amino Acid'] = df_melt_b['Amino Acid'].str.replace('_conc', '', regex=False)
        
        g = sns.FacetGrid(df_melt_b, col='Amino Acid', col_wrap=5, sharey=False, height=2.5, aspect=1)
        g.map_dataframe(sns.boxplot, x='Type', y='LogC', 
                        order=['RRMS', 'SPMS', 'PPMS'],
                        hue='Type', legend=False,
                        palette={'RRMS':'lightgreen', 'SPMS':'orange', 'PPMS':'purple'},
                        showfliers=True,
                        flierprops={'marker': 'o', 'markerfacecolor': 'black', 'markersize': 3})
        
        g.set_titles("{col_name}")
        g.set_xlabels("")
        
        y_ticks = [-5, 0, 5]
        for i, ax in enumerate(g.axes.flat):
            for spine in ax.spines.values():
                spine.set_edgecolor('black')
                spine.set_linewidth(1.5)
            ax.set_yticks(y_ticks)
            if i % 5 != 0: # 5 columns
                ax.set_yticklabels([])
            ax.set_ylabel('')
            
        out_path_1b = os.path.join(OUTPUT_DIR, 'Fig1B_MS_Subtypes.png')
        g.savefig(out_path_1b)
        plt.close()
        print(f"Saved: {out_path_1b}")


        # ==========================================
        # GENERATE JSON DATA
        # ==========================================
        print("Generating JSON data...")

        json_data = {
            'metadata': {
                'title': 'MS Pathology Overview',
                'description': 'Comparison of MS vs Controls (Fig 1A) and MS Subtypes (Fig 1B).'
            },
            'fig1a': {
                'title': 'MS vs Control',
                'traces': []
            },
            'fig1b': {
                'title': 'MS Subtypes',
                'traces': []
            }
        }

        for aa in aa_cols:
            aa_clean = aa.replace('_conc', '')
            
            # --- Fig 1A: MS vs Control ---
            ms_vals = df_p1a[df_p1a['Group']=='MS'][aa].dropna()
            ctrl_vals = df_p1a[df_p1a['Group']=='Control'][aa].dropna()
            
            p_val_a = 1.0
            if len(ms_vals) > 0 and len(ctrl_vals) > 0:
                try:
                    # Mann-Whitney U for MS vs Control
                    _, p_val_a = stats.mannwhitneyu(ms_vals, ctrl_vals, alternative='two-sided')
                except Exception:
                    p_val_a = 1.0
            
            json_data['fig1a']['traces'].append({
                'aa': aa_clean,
                'MS': {'y': ms_vals.tolist()},
                'Control': {'y': ctrl_vals.tolist()},
                'stats': {'p_value': float(p_val_a)}
            })

            # --- Fig 1B: MS Subtypes ---
            rrms_vals = df_p1b[df_p1b['Type']=='RRMS'][aa].dropna()
            spms_vals = df_p1b[df_p1b['Type']=='SPMS'][aa].dropna()
            ppms_vals = df_p1b[df_p1b['Type']=='PPMS'][aa].dropna()
            
            p_val_b = 1.0
            if len(rrms_vals) > 0 and len(spms_vals) > 0 and len(ppms_vals) > 0:
                try:
                    # Kruskal-Wallis for 3 groups
                    _, p_val_b = stats.kruskal(rrms_vals, spms_vals, ppms_vals)
                except Exception:
                    p_val_b = 1.0
            
            json_data['fig1b']['traces'].append({
                'aa': aa_clean,
                'RRMS': {'y': rrms_vals.tolist()},
                'SPMS': {'y': spms_vals.tolist()},
                'PPMS': {'y': ppms_vals.tolist()},
                'stats': {'p_value': float(p_val_b)}
            })

        json_path = os.path.join(OUTPUT_DIR, 'pathology_data.json')
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(json_data, f, indent=2, ensure_ascii=False)
        print(f"Saved JSON: {json_path}")

        # Generate HTML
        html_content = generate_pathology_html(json_filename='pathology_data.json')
        
        html_path = os.path.join(OUTPUT_DIR, 'pathology_overview.html')
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"Saved HTML: {html_path}")

        print("Node 3 completed successfully.")

    except Exception:
        print("CRITICAL ERROR IN GENERATE_FIGURES.PY:")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
