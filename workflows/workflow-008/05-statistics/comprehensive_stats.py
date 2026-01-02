"""Node 5: Comprehensive Statistical Analysis"""
import os
import pandas as pd
import scikit_posthocs as sp
from scipy.stats import mannwhitneyu, chi2_contingency, shapiro, ttest_ind, kruskal

# Setup paths
case_file = "case.pkl"
control_file = "control.pkl"
data_file = "data_cleaned.pkl"
output_csv = "statistics_results.csv"
output_pkl = "statistics_results.pkl"
comprehensive_pkl = "comprehensive_stats.pkl"

# Load parameters from environment variables
DEFAULT_AMINO_ACIDS = "SER,GLN,ARG,CIT,ASN,1MHIS,3MHIS,HYP,GLY,THR,ALA,GABA,SAR,BAIB,ABA,ORN,MET,PRO,LYS,ASP,HIS,VAL,TRP,AAA,LEU,PHE,ILE,C-C,TYR"
AMINO_ACIDS = os.environ.get("PARAM_AMINO_ACIDS", DEFAULT_AMINO_ACIDS).split(",")
print(f"Using {len(AMINO_ACIDS)} amino acids from environment", flush=True)

print("Loading data...", flush=True)
case = pd.read_pickle(case_file)
control = pd.read_pickle(control_file)
df = pd.read_pickle(data_file)

# 1. Mann-Whitney U Tests
print("Performing Mann-Whitney U tests...", flush=True)
mw_results = {}
results_list = []
for var in AMINO_ACIDS:
    case_vals = case[var].dropna()
    control_vals = control[var].dropna()
    stat, pval = mannwhitneyu(case_vals, control_vals, alternative='two-sided')
    mw_results[var] = {'U_statistic': stat, 'p_value': pval, 'significant': pval < 0.05}
    results_list.append({'amino_acid': var, 'U_statistic': stat, 'p_value': pval})

# Save basic results for backward compatibility
results_df = pd.DataFrame(results_list)
results_df['significant'] = results_df['p_value'] < 0.05
results_df = results_df.sort_values('p_value')
results_df.to_csv(output_csv, index=False)
results_df.to_pickle(output_pkl)

# 2. Sex Distribution (Chi-Square)
print("Performing Chi-Square test for sex distribution...", flush=True)
contingency_table = pd.crosstab(df['status'], df['sex'])
chi2, p, dof, expected = chi2_contingency(contingency_table)
sex_stats = {
    'chi2': float(chi2),
    'p_value': float(p),
    'table': contingency_table.to_dict()
}

# 3. Age Analysis (T-Test & Shapiro)
print("Performing Age analysis...", flush=True)
case_age = df[df['status'] == 'case']['age'].dropna()
control_age = df[df['status'] == 'control']['age'].dropna()
shapiro_stat, shapiro_p = shapiro(case_age)
t_stat, t_p = ttest_ind(case_age, control_age, equal_var=False)
age_stats = {
    'shapiro': {'statistic': float(shapiro_stat), 'p_value': float(shapiro_p)},
    't_test': {'statistic': float(t_stat), 'p_value': float(t_p)}
}

# 4. Disease Type Analysis (Kruskal-Wallis & Post-hoc)
print("Performing Disease Type analysis...", flush=True)
melted = df.melt(id_vars=['type'], value_vars=AMINO_ACIDS, value_name='conc')
groups = [group['conc'].dropna() for name, group in melted.groupby('type')]
kw_stat, kw_p = kruskal(*groups)

posthoc_dict = {}
if 'EDSS' in df.columns:
    edss_df = df[df['type'].isin(['PPMS', 'SPMS', 'RRMS'])].dropna(subset=['EDSS'])
    if not edss_df.empty:
        try:
            posthoc = sp.posthoc_conover(edss_df, val_col='EDSS', group_col='type', p_adjust='holm')
            posthoc_dict = posthoc.to_dict()
        except Exception as e:
            posthoc_dict = {"error": str(e)}

disease_stats = {
    'kruskal_wallis_conc': {'statistic': float(kw_stat), 'p_value': float(kw_p)},
    'posthoc_edss': posthoc_dict
}

# Save comprehensive results
comprehensive_results = {
    'mann_whitney': mw_results,
    'sex_distribution': sex_stats,
    'age_analysis': age_stats,
    'disease_types': disease_stats
}

pd.to_pickle(comprehensive_results, comprehensive_pkl)
print(f"Saved: {comprehensive_pkl}", flush=True)
