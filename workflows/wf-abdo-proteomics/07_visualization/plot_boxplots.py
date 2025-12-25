"""Node 7: Visualization - Generate boxplots comparing Case vs Control"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

print("Loading melted data and statistics...", flush=True)
melted = pd.read_pickle("../06_transform/data_melted.pkl")
stats = pd.read_pickle("../05_statistics/statistics_results.pkl")

# Overview boxplot
print("Generating overview boxplot...", flush=True)
plt.figure(figsize=(14, 6))
sns.boxplot(x='AA', y='conc', hue='status', data=melted)
plt.xticks(rotation=45, ha='right')
plt.title('Amino Acid Concentrations: Case vs Control')
plt.xlabel('Amino Acid')
plt.ylabel('Standardized Concentration')
plt.tight_layout()
plt.savefig('boxplot_overview.png', dpi=300)
plt.close()
print("Saved: boxplot_overview.png", flush=True)

# Faceted boxplots
print("Generating faceted boxplots...", flush=True)
sns.set(font_scale=1.0)
g = sns.catplot(
    data=melted, x='status', y='conc', col='AA',
    kind='box', col_wrap=6, height=2.5,
    palette={'case': 'red', 'control': 'blue'}
)
g.set_titles("{col_name}")
plt.tight_layout()
plt.savefig('boxplot_faceted.png', dpi=300)
plt.close()
print("Saved: boxplot_faceted.png", flush=True)

# Print summary
print("\n=== Analysis Complete ===", flush=True)
sig = stats[stats['significant']]
print(f"Significant amino acids ({len(sig)}):", flush=True)
for _, row in sig.iterrows():
    print(f"  {row['amino_acid']}: p={row['p_value']:.4f}", flush=True)
