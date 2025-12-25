"""Subprocess script to generate faceted boxplot (isolated for segfault handling)"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

print("Loading data for faceted plot...", flush=True)
melted = pd.read_pickle("../06_transform/data_melted.pkl")

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
g.savefig('boxplot_faceted.png', dpi=150)
plt.close('all')
print("Saved: boxplot_faceted.png", flush=True)
