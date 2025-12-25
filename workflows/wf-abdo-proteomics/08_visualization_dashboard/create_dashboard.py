"""Node 8: Visualization Dashboard - Generate Data and Plots"""
import os
import sys
import json
import argparse
import subprocess
import http.server
import socketserver
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy.stats import gaussian_kde
import logging
from functools import lru_cache

# Set non-interactive backend
matplotlib.use('Agg')

# Setup paths
current_dir = Path(__file__).parent
root_dir = current_dir.parent
input_file = root_dir / "02_feature_engineering" / "data_cleaned.pkl"
std_file = root_dir / "03_standardization" / "data_standardized.pkl"
stats_file = root_dir / "05_statistics" / "comprehensive_stats.pkl"
adv_file = root_dir / "06_advanced_analysis" / "advanced_results.pkl"
melted_file = root_dir / "07_transform" / "data_melted.pkl"
metadata_file = root_dir / "global_params.json"
output_json = current_dir / "results.json"

# Load global parameters
try:
    with open(metadata_file, 'r') as f:
        config = json.load(f)
        if 'amino_acids_Conc' in config:
            AMINO_ACIDS = config['amino_acids_Conc']
        else:
            AMINO_ACIDS = config.get('amino_acids', [])
        DISEASE_TYPES = config.get('disease_types', [])
        SEXES = config.get('sexes', [])
        KDE_POINTS = config.get('kde_points', 100)
        PLOT_TYPE = config.get('plot_type', 'overview')
except FileNotFoundError:
    print(f"Warning: {metadata_file} not found. Using defaults.")
    AMINO_ACIDS = []
    DISEASE_TYPES = []
    SEXES = []
    KDE_POINTS = 100
    PLOT_TYPE = 'overview'

# Allow environment override for plot type
PLOT_TYPE = os.environ.get("PARAM_PLOT_TYPE", PLOT_TYPE).lower()


# Helper function to clean amino acid names for display
def clean_aa_name(aa):
    """Remove _conc suffix from amino acid names for cleaner display."""
    return aa.replace('_conc', '')


# --- JSON Generation Logic ---

@lru_cache(maxsize=256)
def compute_kde(data_tuple, x_min, x_max, n_points=100):
    data = np.array(data_tuple)
    if len(data) < 2:
        x_range = np.linspace(x_min, x_max, n_points)
        return x_range, np.zeros(n_points)
    try:
        kde = gaussian_kde(data)
        x_range = np.linspace(x_min, x_max, n_points)
        y_range = kde(x_range)
        return x_range, y_range
    except:
        x_range = np.linspace(x_min, x_max, n_points)
        return x_range, np.zeros(n_points)

class FigureGenerator:
    def __init__(self, df, df_std, conc_cols):
        self.df = df
        self.df_std = df_std
        self.conc_cols = conc_cols
        
    def generate_figure1(self):
        print("Generating Figure 1...", flush=True)
        case_data = self.df_std[self.df_std['status'] == 'case']
        control_data = self.df_std[self.df_std['status'] == 'control']
        fig1_data = {}
        for col in self.conc_cols:
            case_vals = case_data[col].dropna().values
            control_vals = control_data[col].dropna().values
            if len(case_vals) > 1 and len(control_vals) > 1:
                x_min = min(case_vals.min(), control_vals.min())
                x_max = max(case_vals.max(), control_vals.max())
                case_x, case_y = compute_kde(tuple(case_vals), x_min, x_max, KDE_POINTS)
                control_x, control_y = compute_kde(tuple(control_vals), x_min, x_max, KDE_POINTS)
                fig1_data[clean_aa_name(col)] = {
                    'case_x': case_x.tolist(), 'case_y': case_y.tolist(),
                    'control_x': control_x.tolist(), 'control_y': control_y.tolist()
                }
        return fig1_data

    def generate_figure2(self, mw_stats):
        print("Generating Figure 2...", flush=True)
        case_data = self.df_std[self.df_std['status'] == 'case']
        control_data = self.df_std[self.df_std['status'] == 'control']
        fig2_data = {}
        for col in self.conc_cols:
            fig2_data[clean_aa_name(col)] = {
                'case': case_data[col].dropna().tolist(),
                'control': control_data[col].dropna().tolist(),
                'p_value': mw_stats.get(col, {}).get('p_value', 1.0)
            }
        return fig2_data

    def generate_figure3(self):
        print("Generating Figure 3...", flush=True)
        self.df_std['mean_conc'] = self.df_std[self.conc_cols].mean(axis=1)
        fig3_data = {}
        for dtype in DISEASE_TYPES:
            subset = self.df_std[self.df_std['type'] == dtype]['mean_conc'].dropna()
            if len(subset) > 0:
                fig3_data[dtype] = subset.tolist()
        return fig3_data

    def generate_figure4(self):
        print("Generating Figure 4...", flush=True)
        valid_data = self.df_std[['age'] + self.conc_cols].dropna(subset=['age'])
        fig4_data = {'age': valid_data['age'].tolist()}
        for col in self.conc_cols:
            fig4_data[clean_aa_name(col)] = valid_data[col].tolist()
        return fig4_data

    def generate_figure6(self):
        print("Generating Figure 6...", flush=True)
        cases = self.df_std[self.df_std['status'] == 'case']
        controls = self.df_std[self.df_std['status'] == 'control']
        fig6_data = {}
        for col in self.conc_cols:
            case_vals = cases[col].dropna().values
            control_vals = controls[col].dropna().values
            if len(case_vals) > 1 and len(control_vals) > 1:
                x_min = min(case_vals.min(), control_vals.min())
                x_max = max(case_vals.max(), control_vals.max())
                x_range, case_y = compute_kde(tuple(case_vals), x_min, x_max, KDE_POINTS)
                _, control_y = compute_kde(tuple(control_vals), x_min, x_max, KDE_POINTS)
                fig6_data[clean_aa_name(col)] = {
                    'x': x_range.tolist(), 'y_case': case_y.tolist(), 'y_control': control_y.tolist()
                }
        return fig6_data

    def generate_figure7(self):
        print("Generating Figure 7...", flush=True)
        fig7_data = {}
        for col in self.conc_cols:
            clean_col = clean_aa_name(col)
            fig7_data[clean_col] = {}
            for sex in SEXES:
                subset = self.df_std[self.df_std['sex'] == sex][col].dropna()
                if len(subset) > 0:
                    fig7_data[clean_col][sex] = subset.tolist()
        return fig7_data

    def generate_figure8(self):
        print("Generating Figure 8...", flush=True)
        fig8_data = {}
        for col in self.conc_cols:
            clean_col = clean_aa_name(col)
            fig8_data[clean_col] = {}
            for dtype in DISEASE_TYPES:
                subset = self.df_std[self.df_std['type'] == dtype][col].dropna()
                if len(subset) > 0:
                    fig8_data[clean_col][dtype] = subset.tolist()
        return fig8_data


def clean_nan_values(obj):
    """Recursively replace NaN and Inf values with None in nested structures."""
    if isinstance(obj, dict):
        return {k: clean_nan_values(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [clean_nan_values(item) for item in obj]
    elif isinstance(obj, float):
        if np.isnan(obj) or np.isinf(obj):
            return None
        return obj
    elif isinstance(obj, (np.integer, np.floating)):
        if isinstance(obj, np.floating) and (np.isnan(obj) or np.isinf(obj)):
            return None
        return obj.item()  # Convert numpy types to Python types
    elif isinstance(obj, np.ndarray):
        return clean_nan_values(obj.tolist())
    elif isinstance(obj, np.bool_):
        return bool(obj)
    return obj

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            if np.isnan(obj) or np.isinf(obj):
                return None
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, np.bool_):
            return bool(obj)
        return super(NpEncoder, self).default(obj)

def run_generate_json():
    print("Loading data for JSON generation...", flush=True)
    if not input_file.exists():
        print(f"Error: {input_file} missing. Run previous steps.", flush=True)
        return

    df = pd.read_pickle(input_file)
    df_std = pd.read_pickle(std_file)
    stats_results = pd.read_pickle(stats_file)
    adv_results = pd.read_pickle(adv_file)

    generator = FigureGenerator(df, df_std, AMINO_ACIDS)
    mw_stats = stats_results.get('mann_whitney', {})

    # Clean amino acid names in statistical tests
    clean_mw_stats = {clean_aa_name(k): v for k, v in mw_stats.items()}
    clean_amino_acids = [clean_aa_name(aa) for aa in AMINO_ACIDS]

    results = {
        'metadata': {
            'total_samples': len(df),
            'total_amino_acids': len(AMINO_ACIDS),
            'amino_acids': clean_amino_acids,
            'case_count': int((df['status'] == 'case').sum()),
            'control_count': int((df['status'] == 'control').sum()),
            'significant_tests': sum(1 for r in mw_stats.values() if r.get('significant', False))
        },
        'fig1': generator.generate_figure1(),
        'fig2': generator.generate_figure2(mw_stats),
        'fig3': generator.generate_figure3(),
        'fig4': generator.generate_figure4(),
        'fig6': generator.generate_figure6(),
        'fig7': generator.generate_figure7(),
        'fig8': generator.generate_figure8(),
        'statistical_tests': clean_mw_stats,
        'additional_stats': {
            'sex_distribution': stats_results.get('sex_distribution'),
            'age_analysis': stats_results.get('age_analysis'),
            'disease_types': stats_results.get('disease_types')
        },
        'advanced_analysis': adv_results
    }

    # Clean all NaN values from the results
    results = clean_nan_values(results)

    with open(output_json, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=2, ensure_ascii=False, cls=NpEncoder)
    print(f"Saved: {output_json.name}", flush=True)


# --- Plotting Logic ---

def generate_overview_plot(melted):
    """Generate simple overview boxplot using matplotlib."""
    print("Generating overview boxplot...", flush=True)
    fig, ax = plt.subplots(figsize=(14, 6))

    aa_list = melted['AA'].unique()
    case_data = [melted[(melted['AA'] == aa) & (melted['status'] == 'case')]['conc'].values for aa in aa_list]
    ctrl_data = [melted[(melted['AA'] == aa) & (melted['status'] == 'control')]['conc'].values for aa in aa_list]

    positions_case = [i - 0.2 for i in range(len(aa_list))]
    positions_ctrl = [i + 0.2 for i in range(len(aa_list))]

    bp1 = ax.boxplot(case_data, positions=positions_case, widths=0.35, patch_artist=True)
    bp2 = ax.boxplot(ctrl_data, positions=positions_ctrl, widths=0.35, patch_artist=True)

    for patch in bp1['boxes']:
        patch.set_facecolor('red')
        patch.set_alpha(0.7)
    for patch in bp2['boxes']:
        patch.set_facecolor('blue')
        patch.set_alpha(0.7)

    ax.set_xticks(range(len(aa_list)))
    ax.set_xticklabels([aa.replace('_conc', '') for aa in aa_list], rotation=45, ha='right')
    ax.set_xlabel('Amino Acid')
    ax.set_ylabel('Standardized Concentration')
    ax.set_title('Amino Acid Concentrations: Case (red) vs Control (blue)')
    ax.legend([bp1['boxes'][0], bp2['boxes'][0]], ['Case', 'Control'], loc='upper right')

    plt.tight_layout()
    output_png = current_dir / 'boxplot_overview.png'
    plt.savefig(output_png, dpi=100)
    plt.close('all')
    print(f"Saved: {output_png.name}", flush=True)

def run_faceted_worker():
    """Worker function to generate faceted plots (isolated process)."""
    print("Loading data for faceted plot...", flush=True)
    melted = pd.read_pickle(melted_file)

    print("Generating faceted boxplots...", flush=True)
    sns.set(font_scale=1.0)
    g = sns.catplot(
        data=melted, x='status', y='conc', col='AA',
        kind='box', col_wrap=6, height=2.0,
        hue='status', palette={'case': 'red', 'control': 'blue'},
        legend=False
    )
    g.set_titles("{col_name}")
    g.tight_layout()
    output_png = current_dir / 'boxplot_faceted.png'
    g.savefig(output_png, dpi=150)
    plt.close('all')
    print(f"Saved: {output_png.name}", flush=True)

def generate_faceted_plot_subprocess():
    """Run faceted plot generation in subprocess (self-invoking)."""
    print("\n" + "-"*60, flush=True)
    print("NOTE: Running faceted plot in subprocess to handle potential memory/segfault issues.", flush=True)
    print("-"*60 + "\n", flush=True)

    try:
        # Calls this same script with --worker flag
        result = subprocess.run(
            [sys.executable, str(Path(__file__).resolve()), "--worker"],
            capture_output=True,
            text=True,
            timeout=120
        )
        if result.stdout:
            print(result.stdout, flush=True)
        if result.returncode != 0:
            print(f"WARNING: Faceted plot generation failed. Exit code: {result.returncode}", flush=True)
            if result.stderr:
                print(f"Error: {result.stderr[:500]}", flush=True)
            return False
        return True
    except subprocess.TimeoutExpired:
        print("WARNING: Faceted plot generation timed out.", flush=True)
        return False
    except Exception as e:
        print(f"WARNING: Faceted plot generation failed: {e}", flush=True)
        return False


def run_plotting():
    print("Loading melted data...", flush=True)
    if not melted_file.exists():
         print(f"Error: {melted_file} missing.", flush=True)
         return

    melted = pd.read_pickle(melted_file)

    if PLOT_TYPE in ['overview', 'both']:
        generate_overview_plot(melted)
    
    if PLOT_TYPE in ['faceted', 'both']:
        generate_faceted_plot_subprocess()

# --- Server Logic ---

def run_server(port=8080):
    handler = http.server.SimpleHTTPRequestHandler
    with socketserver.TCPServer(("", port), handler) as httpd:
        print(f"Serving dashboard at http://localhost:{port}/", flush=True)
        print("Press Ctrl+C to stop.", flush=True)
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nStopping server.", flush=True)


# --- Main ---

def setup_logging():
    log_file = root_dir / "workflow_execution.log"
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(str(log_file)),
            logging.StreamHandler(sys.stdout)
        ]
    )
    # Redirect stdout and stderr to logger
    class StreamToLogger(object):
        def __init__(self, logger, log_level=logging.INFO):
            self.logger = logger
            self.log_level = log_level
            self.linebuf = ''
        def write(self, buf):
            for line in buf.rstrip().splitlines():
                self.logger.log(self.log_level, line.rstrip())
        def flush(self):
            pass
    sys.stdout = StreamToLogger(logging.getLogger('STDOUT'), logging.INFO)
    sys.stderr = StreamToLogger(logging.getLogger('STDERR'), logging.ERROR)

def main():
    setup_logging()
    print("Starting dashboard generation...", flush=True)
    parser = argparse.ArgumentParser(description="Generate Dashboard Data & Plots")
    parser.add_argument("--worker", action="store_true", help="Run as worker process for faceted plots")
    parser.add_argument("--serve", action="store_true", help="Serve the dashboard after generation")
    args = parser.parse_args()

    if args.worker:
        run_faceted_worker()
        return

    # Normal execution
    run_generate_json()
    run_plotting()

    if args.serve:
        run_server()

if __name__ == "__main__":
    main()
