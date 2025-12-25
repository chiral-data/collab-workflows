"""Node 5: Differential Expression Analysis - Mann-Whitney U tests"""
import pandas as pd
from scipy.stats import mannwhitneyu

AMINO_ACIDS = [
    'SER_conc', 'GLN_conc', 'ARG_conc', 'CIT_conc', 'ASN_conc', '1MHIS_conc', '3MHIS_conc', 'HYP_conc', 'GLY_conc',
    'THR_conc', 'ALA_conc', 'GABA_conc', 'SAR_conc', 'BAIB_conc', 'ABA_conc', 'ORN_conc', 'MET_conc', 'PRO_conc',
    'LYS_conc', 'ASP_conc', 'HIS_conc', 'VAL_conc', 'TRP_conc', 'AAA_conc', 'LEU_conc', 'PHE_conc', 'ILE_conc',
    'C-C_conc', 'TYR_conc'
]

print("Loading segmented data...", flush=True)
case = pd.read_pickle("../04_segmentation/case.pkl")
control = pd.read_pickle("../04_segmentation/control.pkl")

print("Performing Mann-Whitney U tests for each amino acid...", flush=True)
results = []
for var in AMINO_ACIDS:
    case_vals = case[var].dropna()
    control_vals = control[var].dropna()
    stat, pval = mannwhitneyu(case_vals, control_vals, alternative='two-sided')
    results.append({'amino_acid': var, 'U_statistic': stat, 'p_value': pval})

results_df = pd.DataFrame(results)
results_df['significant'] = results_df['p_value'] < 0.05
results_df = results_df.sort_values('p_value')

sig_count = results_df['significant'].sum()
print(f"Significant amino acids (p<0.05): {sig_count}", flush=True)
print("\nTop 10 most significant:", flush=True)
print(results_df.head(10).to_string(index=False), flush=True)

results_df.to_csv("statistics_results.csv", index=False)
results_df.to_pickle("statistics_results.pkl")
print("\nSaved: statistics_results.csv, statistics_results.pkl", flush=True)
