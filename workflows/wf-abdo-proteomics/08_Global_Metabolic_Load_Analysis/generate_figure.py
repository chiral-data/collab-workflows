"""Node 8: Global Metabolic Load Analysis - Generate Fig 8 (Python 3)"""
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
    from html_generator import generate_metabolic_load_html
except ImportError as e:
    print(f"Error importing html_generator: {e}")
    sys.exit(1)

def main():
    try:
        # Configuration
        SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
        INPUT_FILE = "data_standardized.pkl" 
        OUTPUT_DIR = "outputs"

        if not os.path.exists(OUTPUT_DIR):
            os.makedirs(OUTPUT_DIR, exist_ok=True)
            print(f"Created output directory: {OUTPUT_DIR}")

        sns.set_style("whitegrid")
        plt.rcParams['font.family'] = 'sans-serif'

        print(">>> NODE 8: GLOBAL METABOLIC LOAD ANALYSIS...")

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

        df_mg_ms = df[masks['MS'] | masks['MG']].copy()
        df_mg_ms['Plot_Type'] = np.where(df_mg_ms['Type'].isin(['general', 'eye-type']), 'GMG', df_mg_ms['Type'])

        # ==========================================
        # FIG 8: Total AA by Type (Static)
        # ==========================================
        print("Generating Fig 8: Total AA by Type...")

        plt.figure(figsize=(8,6))
        type_order = ['PPMS', 'SPMS', 'RRMS', 'GMG']
        type_colors = {'PPMS': '#1f77b4', 'SPMS': '#ff7f0e', 'RRMS': '#2ca02c', 'GMG': '#d62728'}

        sns.boxplot(data=df_mg_ms, x='Plot_Type', y='Total_AA', order=type_order,
                    palette=type_colors, showfliers=True, width=0.5, linewidth=1.5,
                    flierprops={'marker': 'o', 'markerfacecolor': 'white', 'markeredgecolor': 'black', 'markersize': 5})

        plt.ylabel('Concentration [nmol/ml]')
        plt.xlabel('')
        
        # Calculate max y to place simple annotation
        y_max = df_mg_ms['Total_AA'].max()
        plt.ylim(top=y_max * 1.1)

        plt.title('Total Amino Acid Concentration by Disease Type', fontsize=12, fontweight='bold')

        out_path = os.path.join(OUTPUT_DIR, 'Fig8_TotalAA_Type.png')
        plt.savefig(out_path)
        plt.close()
        print(f"Saved: {out_path}")

        # ==========================================
        # STATISTICAL ANALYSIS (Kruskal-Wallis)
        # ==========================================
        print("Calculating Statistics...")
        groups = []
        for t in type_order:
            vals = df_mg_ms[df_mg_ms['Plot_Type'] == t]['Total_AA'].dropna().values
            if len(vals) > 0:
                groups.append(vals)

        p_value = 1.0
        test_name = "N/A"
        if len(groups) > 1:
            try:
                stat, p_value = stats.kruskal(*groups)
                test_name = "Kruskal-Wallis"
                print(f"Kruskal-Wallis p-value: {p_value:.4e}")
            except Exception as e:
                print(f"Stats error: {e}")

        # ==========================================
        # GENERATE JSON DATA
        # ==========================================
        print("Generating JSON...")

        json_data = {
            'metadata': {
                'title': 'Global Metabolic Load Analysis'
            },
            'fig8': {
                'title': 'Total Amino Acid by Disease Type',
                'variable': 'Total_AA',
                'traces': []
            }
        }
        
        # Create a single trace object that contains all group data
        # This structure matches the "Hero Plot" expectation
        trace_data = {
             'name': 'Total AA',
             'stats': {
                 'p_value': float(p_value),
                 'test': test_name
             }
        }
        
        # Add data for each group
        for t in type_order:
            y_vals = df_mg_ms[df_mg_ms['Plot_Type'] == t]['Total_AA'].dropna().tolist()
            trace_data[t] = {'y': y_vals}
            
        json_data['fig8']['traces'].append(trace_data)

        json_path = os.path.join(OUTPUT_DIR, 'metabolic_load_data.json')
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(json_data, f, indent=2, ensure_ascii=False)
        print(f"Saved JSON: {json_path}")

        # ==========================================
        # GENERATE HTML
        # ==========================================
        html_content = generate_metabolic_load_html(json_filename='metabolic_load_data.json')
        
        html_path = os.path.join(OUTPUT_DIR, 'metabolic_load.html')
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"Saved HTML: {html_path}")
        
        print("Node 8 completed successfully.")

    except Exception:
        print("CRITICAL ERROR IN GENERATE_FIGURE.PY:")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
