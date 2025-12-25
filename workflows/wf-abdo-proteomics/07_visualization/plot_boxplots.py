"""Node 7: Visualization - Generate boxplots comparing Case vs Control"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

print("Loading melted data and statistics...", flush=True)
melted = pd.read_pickle("../06_transform/data_melted.pkl")
stats = pd.read_pickle("../05_statistics/statistics_results.pkl")

# Simple overview boxplot only (avoid memory-heavy faceted plots for ARM emulation)
print("Generating overview boxplot...", flush=True)
fig, ax = plt.subplots(figsize=(14, 6))

# Group data for boxplot
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

# Print summary
print("\n=== Analysis Complete ===", flush=True)
sig = stats[stats['significant']]
print(f"Significant amino acids ({len(sig)}):", flush=True)
for _, row in sig.iterrows():
    print(f"  {row['amino_acid']}: p={row['p_value']:.4f}", flush=True)
