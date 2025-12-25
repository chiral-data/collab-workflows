"""Node 7: Visualization - Generate boxplots comparing Case vs Control"""
import os
import sys
import subprocess
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

# Get plot type from environment variable (default: overview)
plot_type = os.environ.get("PARAM_PLOT_TYPE", "overview").lower()

print(f"Plot type: {plot_type}", flush=True)
print("Loading melted data and statistics...", flush=True)
melted = pd.read_pickle("../06_transform/data_melted.pkl")
stats = pd.read_pickle("../05_statistics/statistics_results.pkl")


def generate_overview_plot():
    """Generate simple overview boxplot using matplotlib (ARM-compatible)."""
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
    plt.savefig('boxplot_overview.png', dpi=100)
    plt.close('all')
    print("Saved: boxplot_overview.png", flush=True)


def generate_faceted_plot():
    """Generate faceted boxplots in subprocess (catches segfault gracefully)."""
    print("\n" + "-"*60, flush=True)
    print("NOTE FOR MAC USERS (Apple Silicon M1/M2/M3):", flush=True)
    print("The faceted plot is memory-intensive and may cause crashes", flush=True)
    print("under x86 emulation. If this fails, only 'overview' plot will be used.", flush=True)
    print("-"*60 + "\n", flush=True)

    # Run faceted plot generation in subprocess to catch segfaults
    script_dir = os.path.dirname(os.path.abspath(__file__))
    faceted_script = os.path.join(script_dir, "generate_faceted.py")

    try:
        result = subprocess.run(
            [sys.executable, faceted_script],
            capture_output=True,
            text=True,
            timeout=120  # 2 minute timeout
        )

        # Print subprocess output
        if result.stdout:
            print(result.stdout, flush=True)

        if result.returncode != 0:
            print("\n" + "="*70, flush=True)
            print("WARNING: Faceted plot generation failed!", flush=True)
            print("="*70, flush=True)

            if result.returncode == -11:  # SIGSEGV
                print("\nSegmentation fault detected (common on Mac ARM with x86 emulation)", flush=True)
            else:
                print(f"\nExit code: {result.returncode}", flush=True)

            if result.stderr:
                print(f"Error: {result.stderr[:500]}", flush=True)

            print("\nSOLUTION: Set plot_type='overview' in global_params.json", flush=True)
            print("The overview plot was still generated successfully.", flush=True)
            print("="*70 + "\n", flush=True)
            return False

        return True

    except subprocess.TimeoutExpired:
        print("\nWARNING: Faceted plot generation timed out.", flush=True)
        print("The overview plot was still generated successfully.", flush=True)
        return False
    except Exception as e:
        print(f"\nWARNING: Faceted plot generation failed: {e}", flush=True)
        print("The overview plot was still generated successfully.", flush=True)
        return False


# Generate plots based on plot_type parameter
if plot_type in ['overview', 'both']:
    generate_overview_plot()

if plot_type in ['faceted', 'both']:
    generate_faceted_plot()

# Print summary
print("\n=== Analysis Complete ===", flush=True)
sig = stats[stats['significant']]
print(f"Significant amino acids ({len(sig)}):", flush=True)
for _, row in sig.iterrows():
    print(f"  {row['amino_acid']}: p={row['p_value']:.4f}", flush=True)
