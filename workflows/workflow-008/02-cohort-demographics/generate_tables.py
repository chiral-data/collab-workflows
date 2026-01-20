"""Node 2: Cohort Demographics - Generate Tables & Data (Python 3)"""
import os
import sys
import json
import traceback
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Import local HTML generator
try:
    from html_generator import generate_demographics_html
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

        print(">>> NODE 2: COHORT DEMOGRAPHICS...")

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
            'MG': df['Type'].isin(mg_types),
            'Control': df['Status'] == 'control',
            'RRMS': df['Type'] == 'RRMS',
            'SPMS': df['Type'] == 'SPMS',
            'PPMS': df['Type'] == 'PPMS',
            'GMG': df['Type'] == 'general',
            'OMG': df['Type'] == 'eye-type'
        }

        # ==========================================
        # TABLE 1: MS vs Controls (Static Image)
        # ==========================================
        print("Generating Table 1...")

        ms_data = df[masks['MS']]
        ctrl_data = df[masks['Control']]
        
        # Calculations
        n_ms = len(ms_data)
        ms_males = len(ms_data[ms_data['Sex'] == 'Male'])
        ms_females = len(ms_data[ms_data['Sex'] == 'Female'])
        
        n_ctrl = len(ctrl_data)
        ctrl_males = len(ctrl_data[ctrl_data['Sex'] == 'Male'])
        ctrl_females = len(ctrl_data[ctrl_data['Sex'] == 'Female'])

        # Create Table Data (Simplified for brevity in static export, rich data in JSON)
        table_data = [
            ['Number of subjects (n)', str(n_ms), str(n_ctrl)],
            ['Sex (male/female)', f'{ms_males}/{ms_females}', f'{ctrl_males}/{ctrl_females}'],
            ['Age, years (Mean ± SD)', f"{ms_data['Age'].mean():.1f} ± {ms_data['Age'].std():.1f}", f"{ctrl_data['Age'].mean():.1f} ± {ctrl_data['Age'].std():.1f}"],
            ['Disease duration (Mean ± SD)', f"{ms_data['Duration'].mean():.1f} ± {ms_data['Duration'].std():.1f}", 'n.a.'],
            ['Median EDSS (IQR)', f"{ms_data['EDSS'].median():.1f} ({ms_data['EDSS'].quantile(0.25):.1f}–{ms_data['EDSS'].quantile(0.75):.1f})", 'n.a.'],
            ['DMT: None', f"{len(ms_data[ms_data['DMT_Clean']=='None'])}", 'n.a.']
        ]
        
        df_table = pd.DataFrame(table_data, columns=['Metric', 'MS', 'Controls'])
        
        # Save Static Image (Legacy Support)
        fig, ax = plt.subplots(figsize=(8, len(table_data) * 0.5 + 1))
        ax.axis('tight')
        ax.axis('off')
        the_table = ax.table(cellText=df_table.values, colLabels=df_table.columns, cellLoc='left', loc='center', colWidths=[0.5, 0.25, 0.25])
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(10)
        the_table.scale(1, 1.8)
        
        # Basic Styling
        for (i, j), cell in the_table.get_celld().items():
            if i == 0:
                cell.set_facecolor('#40466e')
                cell.set_text_props(weight='bold', color='white')
            else:
                cell.set_facecolor('#f1f1f2' if i % 2 == 0 else 'white')
        
        plt.title("Table 1. Demographic Overview", y=0.98, weight='bold')
        out_path1 = os.path.join(OUTPUT_DIR, 'Table1_Demographics.png')
        plt.savefig(out_path1, bbox_inches='tight', dpi=150)
        plt.close()


        # ==========================================
        # TABLE 2: MS vs MG (Static Image)
        # ==========================================
        print("Generating Table 2...")
        mg_data = df[masks['MG']]
        n_mg = len(mg_data)
        mg_males = len(mg_data[mg_data['Sex'] == 'Male'])
        mg_females = len(mg_data[mg_data['Sex'] == 'Female'])
        
        table2_data = [
            ['Number of subjects (n)', str(n_ms), str(n_mg)],
            ['Sex (male/female)', f'{ms_males}/{ms_females}', f'{mg_males}/{mg_females}'],
            ['Age (Mean ± SD)', f"{ms_data['Age'].mean():.1f} ± {ms_data['Age'].std():.1f}", f"{mg_data['Age'].mean():.1f} ± {mg_data['Age'].std():.1f}"],
            ['Duration (Mean ± SD)', f"{ms_data['Duration'].mean():.1f} ± {ms_data['Duration'].std():.1f}", f"{mg_data['Duration'].mean():.1f} ± {mg_data['Duration'].std():.1f}"],
        ]
        
        df_table2 = pd.DataFrame(table2_data, columns=['Metric', 'MS', 'MG'])
        
        # Save Static Image
        fig, ax = plt.subplots(figsize=(8, len(table2_data) * 0.5 + 1))
        ax.axis('tight')
        ax.axis('off')
        the_table2 = ax.table(cellText=df_table2.values, colLabels=df_table2.columns, cellLoc='left', loc='center')
        the_table2.auto_set_font_size(False)
        the_table2.set_fontsize(10)
        the_table2.scale(1, 1.8)
        
        for (i, j), cell in the_table2.get_celld().items():
            if i == 0:
                cell.set_facecolor('#40466e')
                cell.set_text_props(weight='bold', color='white')
            else:
                cell.set_facecolor('#f1f1f2' if i % 2 == 0 else 'white')

        plt.title("Table 2. MS vs MG Comparison", y=0.98, weight='bold')
        out_path2 = os.path.join(OUTPUT_DIR, 'Table2_MS_MG_Demographics.png')
        plt.savefig(out_path2, bbox_inches='tight', dpi=150)
        plt.close()


        # ==========================================
        # GENERATE JSON DATA (Rich Data for Dashboard)
        # ==========================================
        print("Generating JSON data...")
        
        # Helper to get distribution data
        def get_dist(data_subset, col):
            return data_subset[col].dropna().tolist()

        json_data = {
            'metadata': {
                'title': 'Cohort Demographics',
                'description': 'Interactive visualization of demographic distributions.'
            },
            'table1': {
                'columns': df_table.columns.tolist(),
                'rows': df_table.values.tolist()
            },
            'table2': {
                'columns': df_table2.columns.tolist(),
                'rows': df_table2.values.tolist()
            },
            # Raw Data for Plots
            'plots': {
                'age': {
                    'MS': get_dist(ms_data, 'Age'),
                    'Control': get_dist(ctrl_data, 'Age'),
                    'MG': get_dist(mg_data, 'Age')
                },
                'duration': {
                    'MS': get_dist(ms_data, 'Duration'),
                    'MG': get_dist(mg_data, 'Duration')
                    # Control has no duration
                },
                'edss': {
                    'MS': get_dist(ms_data, 'EDSS')
                    # MG/Control have no EDSS in this context usually, or n.a.
                },
                'sex': {
                    'MS': {'Male': ms_males, 'Female': ms_females},
                    'Control': {'Male': ctrl_males, 'Female': ctrl_females},
                    'MG': {'Male': mg_males, 'Female': mg_females}
                },
                'subtypes': {
                    'MS': {
                        'RRMS': len(df[masks['RRMS']]),
                        'SPMS': len(df[masks['SPMS']]),
                        'PPMS': len(df[masks['PPMS']])
                    },
                    'MG': {
                        'GMG': len(df[masks['GMG']]),
                        'OMG': len(df[masks['OMG']])
                    }
                }
            }
        }

        json_path = os.path.join(OUTPUT_DIR, 'demographics_data.json')
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(json_data, f, indent=2, ensure_ascii=False)
        print(f"Saved: {json_path}")

        # Generate HTML
        html_content = generate_demographics_html(json_filename="demographics_data.json")
        html_path = os.path.join(OUTPUT_DIR, 'demographics.html')
        with open(html_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        print(f"Saved: {html_path}")

        print("Node 2 completed successfully.")

    except Exception:
        print("CRITICAL ERROR IN GENERATE_TABLES.PY:")
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
