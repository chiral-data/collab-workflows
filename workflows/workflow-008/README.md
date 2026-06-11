---
doc_id: workflow-008
domain: metabolomics
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  Serum amino acid profiling for autoimmune neurological disease differential
  diagnosis. Analyzes LC-MS/MS amino acid concentrations across Multiple
  Sclerosis, Myasthenia Gravis, and healthy controls through an 11-node
  statistical and visualization pipeline.
tags: [metabolomics, amino-acids, multiple-sclerosis, myasthenia-gravis, biomarkers, lc-msms, clinical]
---

# Workflow 008: Serum Amino Acid Profiles in MS and MG

Serum amino acid profiling for differential diagnosis of autoimmune neurological diseases. This workflow processes LC-MS/MS amino acid concentration data from patients with Multiple Sclerosis (MS), Myasthenia Gravis (MG), and healthy controls through an 11-node analytical pipeline, producing publication-quality figures and interactive HTML dashboards for biomarker discovery and clinical decision support.

## Overview

The workflow addresses the diagnostic gray zone between MS and MG — two T-cell-mediated autoimmune diseases that share early symptoms (fatigue, diplopia) but attack different targets (CNS myelin in MS vs. neuromuscular junction in MG). It identifies peripheral metabolic biomarkers (e.g., 3-methylhistidine, citrulline) that can differentiate these conditions when clinical presentation is ambiguous.

The input dataset contains serum concentrations of 29 amino acids and derivatives measured via LC-MS/MS, plus clinical metadata (disease subtype, EDSS score, disease duration, age, sex) for 208 samples across three groups: MS (RRMS, SPMS, PPMS subtypes, n~121), MG (generalized and ocular subtypes, n~27), and healthy controls (n~52). Data were collected under controlled conditions (morning fasting blood draws, steroid washout, relapse-free period) using the EZ:faast LC-MS kit on a Shimadzu LCMS-8045 triple quadrupole (Rzepinski et al., 2023; Koslinski et al., 2024).

Node 01 preprocesses and Z-score standardizes the raw concentrations. Nodes 02–11 run in parallel, each performing a distinct statistical analysis with dual output: static PNG figures (matplotlib/seaborn) and interactive HTML dashboards (Plotly.js).

## When to use this workflow

Use this workflow when you have serum amino acid concentration data and want to perform comparative metabolomics analysis across disease groups, identify potential biomarkers, or validate clinical confounders. The pipeline is designed for the specific MS/MG/control dataset but the analytical framework (Z-score normalization, Mann-Whitney U tests, correlation clustering) applies to similar multi-group metabolomics studies.

Do not use this workflow for protein structure prediction — use workflow-012 (Boltz-2). For protein function annotation, use workflow-015 (mDeepFRI). For ADMET property prediction of small molecules, use workflow-014 (ADMET-AI). For QSAR modeling, use workflow-005.

## Architecture and data flow

```text
                              ┌─> [02: Cohort Demographics]
                              ├─> [03: MS Pathology Overview]
                              ├─> [04: Clinical Confounders]
                              ├─> [05: MS vs MG Autoimmune]
[01: Data Ingestion] ─────────├─> [06: Differential Diagnosis]
                              ├─> [07: Metabolic Network Clustering]
                              ├─> [08: Global Metabolic Load]
                              ├─> [09: Subtype Trajectories]
                              ├─> [10: Clinical Mimicry Test]
                              └─> [11: Pathway Coherence]
```

Node 01 runs first. Nodes 02–11 all depend only on Node 01 and run in parallel.

## Input requirements

- **Dataset:** `database-multiple-sclerosis-myasthenia.csv` — tab-separated file with 208 rows and 38 columns containing patient demographics, clinical variables, and 29 amino acid concentrations in nmol/mL.
- **Column headers are in Polish** — Node 01 translates them to English equivalents.
- **Key columns:**
  - Demographics: `ID Pacjenta` (patient ID), `wiek` (age), `Plec` (gender), `miejsce` (recruitment site)
  - Clinical: `status` (case/control), `postać` (disease form: RRMS/SPMS/PPMS/general/eye-type), `Czas trwania` (disease duration), `EDSS` (disability score), `Lek` (medication)
  - Biomarkers: 29 amino acid `*_conc` columns (ARG, CIT, ORN, TRP, 3MHIS, 1MHIS, BAIB, GABA, AAA, ASP, GLN, LEU, ILE, VAL, PHE, TYR, MET, C-C, LYS, HIS, THR, PRO, ALA, SER, GLY, ASN, HYP, SAR, ABA)

### Data quality notes

- Some patient IDs appear in multiple rows (technical replicates) — Node 01 handles deduplication
- `postać` is NaN for healthy controls and ~8 cases
- `Lek` (medication) has significant missing values
- Raw concentrations span a 100-fold range (glutamine ~600 vs GABA ~5 nmol/mL), requiring Z-score standardization

### Biomarker panels

| Panel | Amino acids | Biological relevance |
|-------|-------------|---------------------|
| Immune regulators | ARG, CIT, ORN, TRP | NO production, urea cycle, kynurenine pathway |
| Muscle catabolism | 3MHIS, 1MHIS, BAIB, GABA | Myofibrillar breakdown, exercise-induced myokine signaling |
| Excitatory/oxidative | AAA, ASP, GLN | Glutamate excitotoxicity, oxidative stress |
| BCAAs | LEU, ILE, VAL | Energy metabolism, immune cell fuel |

## Workflow nodes

### Node 01: Data Ingestion and Preprocessing

**Goal:** Load, clean, translate, and Z-score standardize raw LC-MS/MS amino acid concentrations.

**Process:** Ingests the CSV with Polish column headers, maps them to English equivalents, calculates `Total_AA` (sum of all 29 amino acids per patient), and applies Z-score standardization (mean=0, SD=1) to all amino acid concentrations. Exports the standardized DataFrame as a pickle file and the list of amino acid column names.

**Scientific notes:** Z-score standardization is essential because raw concentrations span two orders of magnitude. Without it, high-abundance amino acids (glutamine, alanine) dominate statistical models, masking dysregulation in trace metabolites that serve as critical biomarkers (3-MHIS for muscle breakdown, BAIB as a contraction-induced myokine). This matches the methodology of the source publications.

**Outputs:**
- `data_standardized.pkl` — pandas DataFrame with standardized concentrations
- `aa_cols.txt` — list of 29 amino acid column names
- `preprocessing_data.json` — raw vs standardized distribution data
- `preprocessing_report.html` — interactive dashboard with side-by-side boxplots

### Node 02: Cohort Demographics

**Goal:** Generate baseline population statistics for study cohorts.

**Process:** Extracts age, gender, disease duration statistics. Computes mean +/- SD for continuous variables and sex ratios. Generates two demographic tables: MS vs Controls (Paper 1 cohort) and MS vs MG (Paper 2 cohort).

**Scientific notes:** Establishes cohort comparability and validates age/sex matching between groups, critical for ensuring metabolic signatures are disease-specific rather than demographic artifacts.

**Outputs:**
- `Table1_Demographics.png`, `Table2_MS_MG_Demographics.png` — publication-quality demographic tables
- `demographics_data.json`, `demographics.html` — interactive demographic tables

### Node 03: MS Pathology Overview

**Goal:** Identify global metabolic dysregulation in MS vs controls and across MS subtypes.

**Process:** Generates 29-panel boxplot grids: (A) all MS vs healthy controls, (B) stratified by MS subtype (RRMS, SPMS, PPMS). Statistical annotation via Mann-Whitney U tests.

**Scientific notes:** Detects key immune dysregulations — arginine and tryptophan depression indicates kynurenine pathway activation. IFN-gamma-induced IDO (indoleamine 2,3-dioxygenase) catabolizes tryptophan, producing kynurenines that suppress T-cell responses. BAIB elevation in PPMS suggests skeletal muscle wasting in progressive phenotypes.

**Outputs:**
- `Fig1A_MS_vs_Control.png`, `Fig1B_MS_Subtypes.png` — 29-amino-acid boxplot grids
- `pathology_data.json`, `pathology.html` — interactive grid with zoom/pan and statistical overlays

### Node 04: Clinical Confounders and Validation

**Goal:** Validate that metabolic changes are disease-driven rather than artifacts of aging, disease chronicity, or disability.

**Process:** Linear regression analysis of each amino acid against age, disease duration, and EDSS score. Generates scatter plots with regression lines, 95% CI shading, and correlation coefficient/p-value annotations.

**Scientific notes:** Weak age correlations confirm metabolic shifts are not normal aging. Strong EDSS correlations for specific amino acids (e.g., 3-MHIS) validate biomarker relevance to disability progression rather than disease duration alone.

**Outputs:**
- `Fig2A_Age_Grid.png`, `Fig2B_Duration_Grid.png`, `Fig2C_EDSS_Grid.png` — correlation grids
- `confounders_data.json`, `confounders.html` — interactive scatter plots with regression lines

### Node 05: MS vs MG Autoimmune Comparison

**Goal:** Identify shared autoimmune metabolic dysregulation common to both MS and MG.

**Process:** Pools MS and MG patients into a combined autoimmune group. Generates 29-panel boxplot grid comparing (MS+MG) vs controls with Mann-Whitney U tests and FDR correction.

**Scientific notes:** Reveals universal T-cell-mediated autoimmune metabolic patterns (arginine/tryptophan dysregulation) that transcend organ-specific pathology. Arginine depletion reflects iNOS overstimulation consuming arginine for proinflammatory NO production, a hallmark of systemic autoimmune metabolic reprogramming.

**Outputs:**
- `Fig3_MS_MG_vs_Control.png` — shared autoimmune metabolic signature
- `autoimmune_data.json`, `autoimmune.html` — interactive comparison

### Node 06: Differential Diagnosis Biomarkers

**Goal:** Identify amino acids that specifically distinguish MS from MG, with sex-stratified validation.

**Process:** Targeted boxplots for discriminatory amino acids (ARG, PRO, CIT) in MS vs MG (full cohort and female-only subset). Mann-Whitney U tests with Bonferroni correction.

**Scientific notes:** Citrulline, a byproduct of NOS-mediated NO synthesis and a key urea cycle intermediate, shows differential levels between MS and MG. Female-only analysis controls for sex-hormone effects on amino acid metabolism.

**Outputs:**
- `Fig4_Specific_Diffs.png`, `Fig5_Female_Specific.png` — biomarker panels
- `biomarkers_data.json`, `biomarkers.html` — interactive biomarker explorer

### Node 07: Metabolic Network Clustering

**Goal:** Visualize metabolic pathway dysregulation via unsupervised hierarchical clustering.

**Process:** Computes Pearson correlation matrices (29x29 amino acids) separately for patients (MS+MG) and healthy controls. Generates clustered heatmaps with dendrograms.

**Scientific notes:** Demonstrates "loss of homeostasis" in autoimmunity — controls show tight metabolic coordination (strong positive correlations), while patients exhibit weakened, fragmented correlation structures indicating dysregulated metabolic networks.

**Outputs:**
- `Fig6_Corr_MS_MG.png`, `Fig7_Corr_Control.png` — clustered correlation heatmaps
- `clustering_data.json`, `clustering.html` — interactive heatmaps with dendrogram navigation

### Node 08: Global Metabolic Load Analysis

**Goal:** Assess total amino acid burden across disease subtypes.

**Process:** Boxplot comparison of total serum amino acid concentrations stratified by MS subtypes (RRMS, SPMS, PPMS), MG subtypes (generalized, ocular), and controls. Kruskal-Wallis H-test with post-hoc pairwise comparisons.

**Scientific notes:** Validates hypermetabolism in progressive phenotypes — SPMS and PPMS show elevated total amino acid levels, suggesting increased protein catabolism or impaired peripheral utilization in advanced neurodegeneration.

**Outputs:**
- `Fig8_TotalAA_Type.png` — total amino acid load by subtype
- `metabolic_load_data.json`, `metabolic_load.html` — interactive load comparison

### Node 09: Subtype Trajectories

**Goal:** Simulate longitudinal disease progression by plotting amino acid concentrations vs disease duration.

**Process:** Scatter plots of disease duration (x-axis) vs each of 29 amino acids (y-axis) with separate regression lines per MS/MG subtype.

**Scientific notes:** Amino acids with steep positive slopes in progressive subtypes (SPMS/PPMS) indicate accelerating metabolic dysregulation over time. Stable trajectories in RRMS suggest relapse-remission metabolic resilience. This cross-sectional analysis simulates longitudinal trends.

**Outputs:**
- `Fig9_Duration_Grid.png` — 29-panel duration trajectory grid
- `trajectories_data.json`, `trajectories.html` — interactive trajectory explorer with subtype filtering

### Node 10: Clinical Mimicry Test (RRMS vs MG)

**Goal:** Distinguish early-stage MS (RRMS) from MG in the diagnostic gray zone.

**Process:** Total amino acid load comparison and 29-panel amino acid grid comparing RRMS vs MG specifically. Mann-Whitney U tests for each amino acid.

**Scientific notes:** 3-Methylhistidine (3-MHIS) emerges as a critical early discriminator. 3-MHIS is formed by post-translational methylation of histidine in actin and myosin and cannot be reused for protein synthesis, making it a specific marker of myofibrillar protein breakdown. It is elevated in MG due to increased muscle proteolysis at the neuromuscular junction but remains normal in RRMS where pathology is confined to CNS myelin.

**Outputs:**
- `Fig10_TotalAA_RRMS_MG.png`, `Fig11_RRMS_vs_MG_Grid.png` — early-stage differentiation
- `mimicry_data.json`, `mimicry.html` — interactive diagnostic decision-support tool

### Node 11: Pathway Coherence Analysis

**Goal:** Contrast metabolic network organization between MS and MG.

**Process:** Side-by-side triangular correlation heatmaps (MS vs MG) with identical amino acid ordering for direct visual comparison. Pearson correlation with hierarchical clustering.

**Scientific notes:** MS shows a highly fragmented correlation structure (reflecting diffuse CNS demyelination effects on systemic metabolism) while MG retains more organized peripheral metabolic networks (reflecting localized neuromuscular pathology). Validates that MS and MG exhibit distinct systemic metabolic architectures despite shared autoimmune mechanisms.

**Outputs:**
- `Fig12_Split_Corr.png` — split-panel correlation comparison
- `coherence_data.json`, `coherence.html` — interactive side-by-side network explorer

## Parameters

### data_file

- **Type:** string
- **Default:** `"database-multiple-sclerosis-myasthenia.csv"`
- **Node:** 01
- **Description:** Input CSV/TSV file with raw amino acid concentrations and clinical metadata.

Nodes 02–11 have no user-configurable parameters. All analysis settings (statistical tests, correction methods, figure layouts) are hardcoded for reproducibility of the published results.

## Outputs and interpretation

### Output structure

Each analysis node (02–11) generates three output types:
- **PNG files** — publication-ready figures (300 DPI, matplotlib/seaborn)
- **JSON files** — structured plot data (coordinates, statistics, metadata) for web rendering
- **HTML files** — self-contained interactive dashboards (Plotly.js 2.27.0 with zoom, pan, hover tooltips)

Total: 38 files per complete run (4 from Node 01 + 18 PNG + 10 JSON + 10 HTML from Nodes 02–11).

### Output summary by node

| Node | PNG outputs | JSON/HTML outputs |
|------|-------------|-------------------|
| 01 | — | preprocessing_data.json, preprocessing_report.html |
| 02 | Table1_Demographics, Table2_MS_MG_Demographics | demographics_data.json, demographics.html |
| 03 | Fig1A_MS_vs_Control, Fig1B_MS_Subtypes | pathology_data.json, pathology.html |
| 04 | Fig2A_Age_Grid, Fig2B_Duration_Grid, Fig2C_EDSS_Grid | confounders_data.json, confounders.html |
| 05 | Fig3_MS_MG_vs_Control | autoimmune_data.json, autoimmune.html |
| 06 | Fig4_Specific_Diffs, Fig5_Female_Specific | biomarkers_data.json, biomarkers.html |
| 07 | Fig6_Corr_MS_MG, Fig7_Corr_Control | clustering_data.json, clustering.html |
| 08 | Fig8_TotalAA_Type | metabolic_load_data.json, metabolic_load.html |
| 09 | Fig9_Duration_Grid | trajectories_data.json, trajectories.html |
| 10 | Fig10_TotalAA_RRMS_MG, Fig11_RRMS_vs_MG_Grid | mimicry_data.json, mimicry.html |
| 11 | Fig12_Split_Corr | coherence_data.json, coherence.html |

### Statistical methods

- **Group comparisons:** Mann-Whitney U test (non-parametric, does not assume normality)
- **Multiple testing correction:** FDR (Nodes 03, 05) or Bonferroni (Node 06) depending on analysis type
- **Correlations:** Pearson correlation with 95% CI for confounder analysis (Node 04) and network clustering (Nodes 07, 11)
- **Multi-group comparisons:** Kruskal-Wallis H-test with post-hoc pairwise tests (Node 08)

## Quick start

### Running with Docker

All nodes use the same container image:

```bash
docker pull ghcr.io/chiral-data/proteomics:2025_12_31
```

### Running on Silva

1. Select "MS_MG Amino Acid Analysis" from the workflow list
2. Ensure the input CSV is present (or use the bundled dataset)
3. Click Run — Node 01 preprocesses data, then Nodes 02–11 run in parallel

### Running individual nodes

```bash
cd 03-ms-pathology-overview
docker run -v $(pwd):/workspace ghcr.io/chiral-data/proteomics:2025_12_31 bash run.sh
```

Execution time: ~2–3 minutes for the complete workflow (with cached Docker image).

## Dataset documentation

### Cohort composition

- **Multiple Sclerosis (MS):** ~121 unique patients. Subtypes: RRMS (relapsing-remitting), SPMS (secondary progressive), PPMS (primary progressive). Median EDSS ~6.0.
- **Myasthenia Gravis (MG):** ~27 unique patients. Subtypes: generalized and ocular.
- **Healthy controls:** ~52 participants with no evidence of CNS or peripheral nervous system disorders.

### Sample collection

- Blood drawn 7–9 AM following overnight fast
- Stored at -80C until batch analysis
- Participants avoided supplementation and protein-rich meals for 7 days prior
- Excluded: patients on steroids within 3 months, patients with recent relapse within 3 months

### Analytical chemistry

- **Kit:** EZ:faast LC-MS Free Amino Acid Analysis Kit (Phenomenex)
- **Platform:** Shimadzu Nexera XR HPLC + LCMS-8045 Triple Quadrupole
- **Ionization:** Positive ion mode ESI+
- **Sample preparation:** Solid-phase extraction, chemical derivatization, liquid-liquid extraction
- **Internal standards:** Homoarginine, methionine-d3, homophenylalanine

### Data provenance

The dataset is an educational subset curated by Prof. Emilia Daghir-Wojtkowiak from previously published clinical research. Formal permission has been granted for use in this workflow.

## References

- Rzepinski, L., Koslinski, P., Kowalewski, M., Koba, M. & Maciejek, Z. "Serum amino acid profiling in differentiating clinical outcomes of multiple sclerosis." *Neurologia i Neurochirurgia Polska* 57(5):414–422, 2023. DOI: https://doi.org/10.5603/PJNNS.a2023.0054
- Koslinski, P., Rzepinski, L., Koba, M., Maciejek, Z., Kowalewski, M. & Daghir-Wojtkowiak, E. "Comparative Analysis of Serum Amino Acid Profiles in Patients with Myasthenia Gravis and Multiple Sclerosis." *J. Clin. Med.* 13(14):4083, 2024. DOI: https://doi.org/10.3390/jcm13144083
