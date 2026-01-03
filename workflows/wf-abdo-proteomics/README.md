# Serum Amino Acid Profiles in MS and MG

This project analyzes clinical and metabolomics data to identify metabolic differences between patients with **Multiple Sclerosis (MS)**, **Myasthenia Gravis (MG)**, and healthy controls. The workflow processes raw LC-MS/MS amino acid profiles through an 11-node analytical pipeline to uncover significant biomarkers distinguishing autoimmune neurological diseases.

---

## Part 1: Dataset Documentation

## 1. Overview
This dataset contains serum amino acid profiles and clinical metadata for patients with Multiple Sclerosis (MS), Myasthenia Gravis (MG), and healthy controls. The data supports research into metabolic biomarkers for autoimmune neurological diseases and has been utilized in two published studies focusing on clinical outcome differentiation and comparative disease profiling.

### Clinical Rationale
This dataset specifically targets the **"Diagnostic Gray Zone"** between Multiple Sclerosis (MS) and Myasthenia Gravis (MG).

* **Clinical Mimicry**: Both diseases share early symptoms like fatigue and diplopia (double vision), leading to potential misdiagnosis.
* **Shared Pathology**: Both are T-cell mediated autoimmune disorders, yet they attack different targets (Myelin in MS vs. Neuromuscular Junction in MG).
* **Goal**: The dataset is designed to identify peripheral metabolic biomarkers (e.g., 3-MHIS vs. Citrulline) that can differentiate these conditions when clinical presentation is ambiguous.

### Key Features
* **Disease Groups**: Multiple Sclerosis (RRMS, SPMS, PPMS), Myasthenia Gravis (Generalized, Ocular), and Healthy Controls.
* **Biomarkers**: Concentrations of 29 amino acids and derivatives measured via LC-MS/MS.
* **Clinical Data**: Disease duration, EDSS scores, age, gender, and medication history.

## 2. File Information
* **Filename**: `database-multiple-sclerosis-myasthenia.csv`
* **Format**: Tab-separated values (TSV/CSV)
* **Rows**: 208
* **Columns**: 38
* **Encoding**: UTF-8 (implied)

## 3. Data Dictionary
The dataset uses Polish column headers. The table below provides the mapping to English descriptions and variable types.

### Identifiers & Demographics
| Column Name (Raw) | Description | Data Type | Notes |
| :--- | :--- | :--- | :--- |
| `ID Pacjenta` | Patient ID | String | Unique identifier (e.g., "1 Sanitas"). *Note: Some IDs appear in multiple rows (replicates).* |
| `status` | Group Status | String | Categories: `case` (Patients), `control` (Healthy Controls). |
| `miejsce` | Recruitment Site | String | e.g., "Sanitas", "Borne". |
| `wiek` | Age | Float | Age in years (Range: 19–81). |
| `Plec` | Gender | String | `Female`, `Male`. |

### Clinical Variables
| Column Name (Raw) | Description | Data Type | Notes |
| :--- | :--- | :--- | :--- |
| `postać` | Disease Form | String | **MS Types**: `RRMS` (Relapsing-Remitting), `SPMS` (Secondary Progressive), `PPMS` (Primary Progressive)<br>**MG Types**: `general` (Generalized), `eye-type` (Ocular) |
| `Czas trwania` | Disease Duration | Float | Duration of the disease in years. |
| `EDSS` | EDSS Score | Float | Expanded Disability Status Scale (MS only). |
| `Lek` | Medication | String | Disease-modifying therapies or other drugs (e.g., "mestinon", "Tysabri"). |

### Serum Biomarker Panels (nmol/mL)
The following columns represent serum concentrations determined by LC-MS/MS, grouped by biological function.

**Panel A: Immune Regulators & Urea Cycle**
*Relevant for NO production and immune suppression.*
* `ARG_conc` (Arginine) – Precursor for Nitric Oxide (NO)
* `CIT_conc` (Citrulline) – Marker of urea cycle and NO synthesis
* `ORN_conc` (Ornithine) – Product of Arginase activity
* `TRP_conc` (Tryptophan) – Modulates T-cell response via kynurenine pathway

**Panel B: Muscle Metabolism & Catabolism**
*Relevant for Myasthenia Gravis and progressive MS neuromuscular atrophy.*
* `3MHIS_conc` (3-Methyl-L-histidine) – Specific marker of myofibrillar protein breakdown
* `1MHIS_conc` (1-Methyl-L-histidine) – Anserine/Carnosine metabolite
* `BAIB_conc` (Beta-aminoisobutyric acid) – "Myokine" released by contracting skeletal muscle
* `GABA_conc` (Gamma-aminobutyric acid) – Neuromuscular tone regulator

**Panel C: Excitatory & Oxidative Stress Markers**
*Relevant for neurodegeneration and excitotoxicity.*
* `AAA_conc` (α-Aminoadipic acid) – Modulator of glutamate uptake; linked to oxidative stress
* `ASP_conc` (Aspartic acid) – Excitatory neurotransmitter
* `GLN_conc` (Glutamine) – Major fuel for immune cells

**Panel D: Standard Amino Acid Profile**
*Branched-Chain Amino Acids (BCAAs):*
* `LEU_conc` (Leucine)
* `ILE_conc` (Isoleucine)
* `VAL_conc` (Valine)

*Aromatic & Sulfur-Containing:*
* `PHE_conc` (Phenylalanine)
* `TYR_conc` (Tyrosine)
* `MET_conc` (Methionine)
* `C-C_conc` (Cysteine/Cystine)

*Other Essential & Non-Essential:*
* `LYS_conc` (Lysine)
* `HIS_conc` (Histidine)
* `THR_conc` (Threonine)
* `PRO_conc` (Proline)
* `ALA_conc` (Alanine)
* `SER_conc` (Serine)
* `GLY_conc` (Glycine)
* `ASN_conc` (Asparagine)

*Other Metabolites:*
* `HYP_conc` (4-Hydroxyproline)
* `SAR_conc` (Sarcosine)
* `ABA_conc` (α-Aminobutyric acid)

## 4. Cohort Composition

1.  **Multiple Sclerosis (MS)** (n=121 unique patients)
    * **Subtypes**: Secondary Progressive (SPMS), Relapsing-Remitting (RRMS), Primary Progressive (PPMS).
    * **Clinical**: Median EDSS ~6.0.
2.  **Myasthenia Gravis (MG)** (n=27-28 unique patients)
    * **Subtypes**: Generalized MG (`general`) and Ocular MG (`eye-type`).
3.  **Healthy Controls** (n=52-53 unique participants)
    * No evidence of central or peripheral nervous system disorders.

## 5. Methodological Notes & Study Design

### Sample Collection Standards
* **Timing**: Blood samples were drawn in the morning (7–9 a.m.) following an overnight fast.
* **Storage**: Samples were immediately centrifuged and stored at -80°C until batch analysis.
* **Dietary Restrictions (Crucial)**: Participants were explicitly instructed to maintain current dietary habits but avoid supplementation and protein-rich meals for 7 days prior to collection. This control minimizes diet-induced noise in the amino acid profiles.

### Exclusion Criteria
* **Steroid Washout**: Patients treated with steroids within 3 months preceding the blood draw were excluded to prevent drug-induced metabolic shifts.
* **Relapse-Free**: Patients experiencing a disease relapse within 3 months were excluded to capture baseline metabolic state.
* **Neurological Controls**: Healthy controls were screened to exclude any evidence of central or peripheral nervous system disorders.

### Analytical Chemistry
* **Sample Preparation**: A 3-step procedure involving:
    1. Solid-phase extraction (SPE) for protein removal
    2. Chemical derivatization to enhance ionization
    3. Liquid-liquid extraction to separate amino acids from interfering compounds
* **Instrumentation**:
    * **Kit**: EZ:faast™ LC-MS Free Amino Acid Analysis Kit (Phenomenex)
    * **Platform**: Shimadzu Nexera XR HPLC coupled with LCMS-8045 Triple Quadrupole
    * **Ionization**: Positive ion mode electrospray ionization (ESI+)
* **Quantification**: Internal standard method using:
    * Homoarginine
    * Methionone-d3 (deuterated)
    * Homophenylalanine

## 6. Data Quality & Usage Notes

### Technical Replicates
* **Observation**: The dataset contains 208 rows. Inspection reveals that some `ID Pacjenta` values are duplicated (e.g., ID "37 Sanitas" appears twice). These rows share identical demographic/clinical data but have varying amino acid concentrations, indicating technical replicates or repeated measures.
* **Recommendation**: Users should handle these duplicates (e.g., by averaging) before statistical analysis to avoid pseudoreplication bias.

### Missing Values
* `postać` (Disease Form) is `NaN` for healthy controls and a small subset of cases (n=8).
* `Lek` (Medication) contains significantly missing values, indicating either no treatment or incomplete medication history.
* **Impact**: Medication-stratified analyses may have reduced statistical power.

### Standardization Status
* **Current State**: The values in this dataset are **Raw Concentrations (nmol/mL)**.
* **Publication Methods**: The source publications explicitly state that data was "centered and standardized prior to analysis" to account for the high variance between high-abundance amino acids (like Glutamine ~600 nmol/mL) and trace metabolites (like GABA ~50 nmol/mL).
* **Action Required**: Users must apply **Z-score standardization** (as per Pipeline Node 1) to replicate the statistical significance and effect sizes reported in the referenced studies. Without standardization, low-abundance biomarkers will be statistically underweighted.

---

## Part 2: Technical Documentation (Pipeline Architecture)

The analysis is implemented as a **modular, parallelizable 11-node workflow** orchestrated by the `silva` command-line tool. It features dual output generation (static PNG visualizations + interactive JSON/HTML web presentations) and is designed for both local execution and remote deployment via TOML-based configuration.

### Architecture Overview

**Execution Model:**
* **Orchestration**: The `silva` tool manages the workflow execution based on dependencies defined in the node-level `job.toml` files.
* **Sequential Foundation**: Node 01 (Data Ingestion) is executed first as other nodes depend on its output.
* **Parallel Analysis**: Nodes 02-11 are executed in parallel after Node 01 completes.
* **Dual Output System**: Each analysis node generates:
  - **Static PNG images** for publication-quality figures.
  - **JSON data files** storing plot data in a structured format.
  - **Interactive HTML files** with Plotly.js visualizations for web-based exploration.

**Configuration System:**
* **Workflow Definition**: The workflow is defined by the collection of `job.toml` files within each node's directory. The `silva` tool discovers and executes these jobs.
* **Node Configuration**: Each node has a `job.toml` file specifying its inputs, outputs, and execution parameters.

**Technology Stack:**
* **Core**: Python 3.11 with a virtual environment (`.venv`).
* **Analysis**: pandas, numpy, scipy (statistical analysis).
* **Visualization**: matplotlib, seaborn (static plots), Plotly.js 2.27.0 (interactive).
* **Web Technologies**: Plain JavaScript ES6+ (async/await, classes, arrow functions).
* **Modular Architecture**: Each analysis node (02-11) contains its own `html_generator.py` file for self-contained HTML generation.
  - **Benefits**: Easier debugging, isolated development, and no cross-node dependencies.
  - **Size**: ~118 lines, ~12 KB per generator.
  - **Customization**: Node-specific rendering functions (e.g., `renderGridPlot()`, `renderBoxPlot()`, `renderTable()`).

---

### **Node 01: Data Ingestion and Preprocessing**

**Purpose:** Loads, cleans, and standardizes raw LC-MS/MS data to prepare for cross-group statistical comparison.

**Input:**
* `database-multiple-sclerosis-myasthenia.csv` (raw amino acid concentrations)

**Processing:**
1. **Data Loading**: Ingests CSV/TSV with Polish column headers
2. **Translation**: Maps Polish clinical variables to English equivalents
3. **Feature Engineering**: Calculates `Total_AA` (sum of all 29 amino acids per patient)
4. **Z-Score Standardization**: Transforms all amino acid concentrations to Mean=0, SD=1
   - Rationale: Accounts for 100-fold concentration differences (Glutamine ~600 vs GABA ~5 nmol/mL)
5. **Data Export**: Saves standardized DataFrame and amino acid column list

**Outputs:**
* `outputs/data_standardized.pkl` (65KB) - Pandas DataFrame with standardized concentrations
* `outputs/aa_cols.txt` (294B) - List of 29 amino acid column names for downstream processing

**Script:** `load_data.py` | **Execution:** `run.sh` via `.venv\Scripts\python.exe`

**Biological Context:** Without standardization, high-abundance amino acids dominate statistical models, masking dysregulation in trace metabolites that serve as critical biomarkers (e.g., 3-MHIS for muscle breakdown, GABA for neuromuscular tone).

---

### **Node 02: Cohort Demographics**

**Purpose:** Generates baseline population statistics for both study cohorts (MS vs Controls, MS vs MG).

**Input:**
* `../01_Data_Ingestion_and_Preprocessing/outputs/data_standardized.pkl`

**Processing:**
* Extracts Age, Gender, Disease Duration statistics
* Computes Mean ± SD for continuous variables
* Calculates Sex ratios (Female:Male)
* Generates styled dataframe tables

**Outputs:**
* `outputs/Table1_Demographics.png` - MS vs Controls demographics (Paper 1 cohort)
* `outputs/Table2_MS_MG_Demographics.png` - MS vs MG demographics (Paper 2 cohort)
* `outputs/demographics_data.json` (2.7KB) - Structured table data for web display
* `outputs/demographics.html` (7.7KB) - Interactive demographic tables with hover tooltips

**Script:** `generate_tables.py` | **Configuration:** `job.toml` | **HTML Generator:** `html_generator.py`

**HTML Implementation Details:**
* Self-contained module with `generate_demographics_html(json_filename='demographics_data.json')`
* Custom `renderTable()` function for demographic table rendering
* No dependencies on other nodes

**Biological Context:** Establishes cohort heterogeneity and validates age/sex matching between groups, critical for identifying disease-specific metabolic signatures independent of demographic confounders.

---

### **Node 03: MS Pathology Overview**

**Purpose:** Identifies global metabolic dysregulation in MS vs Controls and phenotypic differences across MS subtypes.

**Input:**
* `data_standardized.pkl`
* `aa_cols.txt` (29 amino acids)

**Processing:**
* **Figure 1A**: 29-panel boxplot grid comparing MS patients (all subtypes) vs Healthy Controls
* **Figure 1B**: 29-panel boxplot grid stratifying by MS subtype (RRMS, SPMS, PPMS)
* Statistical annotation via Mann-Whitney U tests
* JSON export includes all datapoints, medians, quartiles, and p-values

**Outputs:**
* `outputs/Fig1A_MS_vs_Control.png` (64KB) - Global MS metabolic shifts
* `outputs/Fig1B_MS_Subtypes.png` (67KB) - Subtype-specific patterns
* `outputs/pathology_data.json` (328KB) - Complete boxplot data for 58 comparisons (29 AA × 2 figures)
* `outputs/pathology.html` (9.4KB) - Interactive grid with zoom/pan, statistical overlays

**Script:** `generate_figures.py` | **HTML Generator:** `html_generator.py`

**HTML Implementation Details:**
* Self-contained module with `generate_pathology_html(json_filename='pathology_data.json')`
* Custom `renderGridPlot()` function with 5-column layout optimized for 29 amino acids
* CONFIG object includes color schemes for MS subtypes and statistical annotations

**Biological Context:** Detects key immune dysregulations (ARG↓, TRP↓ indicating kynurenine pathway activation) and identifies progressive subtype markers (BAIB elevation in PPMS suggesting muscle wasting).

---

### **Node 04: Clinical Confounders and Validation**

**Purpose:** Validates that metabolic changes are disease-driven rather than artifacts of aging, disease chronicity, or physical disability.

**Input:**
* `data_standardized.pkl`
* `aa_cols.txt`

**Processing:**
* Linear regression analysis: Age vs 29 amino acids
* Linear regression analysis: Disease Duration vs 29 amino acids
* Linear regression analysis: EDSS (disability score) vs 29 amino acids
* Scatter plots with regression lines and 95% CI shading
* Correlation coefficient (r) and p-value annotations

**Outputs:**
* `outputs/Fig2A_Age_Grid.png` - Age correlation grid (29 amino acids)
* `outputs/Fig2B_Duration_Grid.png` - Disease duration correlation grid
* `outputs/Fig2C_EDSS_Grid.png` - Disability correlation grid
* `outputs/confounders_data.json` (537KB) - Scatter plot coordinates, regression parameters, statistics
* `outputs/confounders.html` (9.5KB) - Interactive scatter plots with dynamic regression lines

**Script:** `generate_figures.py`

**Biological Context:** Distinguishes primary disease mechanisms from secondary effects. For example, weak Age correlations confirm metabolic shifts are not normal aging; strong EDSS correlations for specific amino acids (e.g., 3-MHIS) validate biomarker relevance to disability progression.

---

### **Node 05: MS vs MG Autoimmune Comparison**

**Purpose:** Identifies shared "autoimmune background" dysregulation common to both MS and MG when compared to healthy controls.

**Input:**
* `data_standardized.pkl`
* `aa_cols.txt`

**Processing:**
* Pools MS and MG patients into combined "Autoimmune" group
* Generates 29-panel boxplot grid: (MS+MG) vs Controls
* Mann-Whitney U tests with FDR correction

**Outputs:**
* `outputs/Fig3_MS_MG_vs_Control.png` - Shared autoimmune metabolic signature
* `outputs/autoimmune_data.json` - Boxplot data with statistical annotations
* `outputs/autoimmune.html` - Interactive comparison with toggle controls

**Script:** `generate_figure.py`

**Biological Context:** Reveals universal T-cell-mediated autoimmune metabolic patterns (e.g., ARG/TRP dysregulation) that transcend organ-specific pathology, supporting the concept of systemic metabolic reprogramming in autoimmunity.

---

### **Node 06: Differential Diagnosis Biomarkers**

**Purpose:** Identifies amino acids that specifically distinguish central demyelination (MS) from neuromuscular autoimmunity (MG), with sex-stratified validation.

**Input:**
* `data_standardized.pkl`

**Processing:**
* **Figure 4**: Targeted boxplots for discriminatory amino acids (ARG, PRO, CIT) in MS vs MG
* **Figure 5**: Same analysis restricted to female patients (controls for sex-hormone effects)
* Mann-Whitney U tests with Bonferroni correction

**Outputs:**
* `outputs/Fig4_Specific_Diffs.png` - Key biomarker panel (full cohort)
* `outputs/Fig5_Female_Specific.png` - Female-only validation
* `outputs/biomarkers_data.json` - Boxplot data for targeted amino acids
* `outputs/biomarkers.html` - Interactive biomarker explorer

**Script:** `generate_figures.py`

**Biological Context:** Highlights **Citrulline** (urea cycle, lower in MS) and **GABA** (neuromuscular tone, higher in MG) as clinical decision-support biomarkers for ambiguous early-stage presentations (e.g., fatigue + diplopia).

---

### **Node 07: Metabolic Network Clustering**

**Purpose:** Visualizes metabolic pathway dysregulation via unsupervised hierarchical clustering of amino acid correlations.

**Input:**
* `data_standardized.pkl`

**Processing:**
* Computes Pearson correlation matrices (29×29 amino acids)
* Generates clustered heatmaps separately for:
  - MS+MG patients
  - Healthy controls
* Hierarchical clustering reveals disrupted metabolic coordination

**Outputs:**
* `outputs/Fig6_Corr_MS_MG.png` - Patient correlation network
* `outputs/Fig7_Corr_Control.png` - Control correlation network
* `outputs/clustering_data.json` - Correlation matrices with cluster assignments
* `outputs/clustering.html` - Interactive heatmap with dendrogram navigation

**Script:** `generate_figures.py`

**Biological Context:** Demonstrates "loss of homeostasis" in autoimmunity – controls show tight metabolic coordination (strong positive correlations), while patients exhibit chaotic, weakened correlation structures indicating dysregulated metabolic networks.

---

### **Node 08: Global Metabolic Load Analysis**

**Purpose:** Assesses total amino acid burden across disease subtypes.

**Input:**
* `data_standardized.pkl` (`Total_AA` column)

**Processing:**
* Boxplot comparison of total serum amino acid concentrations
* Stratified by: MS subtypes (RRMS, SPMS, PPMS), MG subtypes (Generalized, Ocular), Controls
* Kruskal-Wallis H-test with post-hoc pairwise comparisons

**Outputs:**
* `outputs/Fig8_TotalAA_Type.png` - Total AA load by subtype
* `outputs/metabolic_load_data.json` - Boxplot data with statistics
* `outputs/metabolic_load.html` - Interactive load comparison

**Script:** `generate_figure.py`

**Biological Context:** Validates observations of hypermetabolism in progressive phenotypes (SPMS, PPMS show elevated total AA), suggesting increased protein catabolism or impaired peripheral utilization in advanced neurodegeneration.

---

### **Node 09: Subtype Trajectories**

**Purpose:** Simulates longitudinal disease progression by plotting amino acid concentrations vs disease duration, stratified by subtype.

**Input:**
* `data_standardized.pkl`
* `aa_cols.txt`

**Processing:**
* Scatter plots: Disease Duration (x-axis) vs 29 amino acids (y-axis)
* Separate regression lines for each MS/MG subtype
* Identifies amino acids with subtype-divergent trajectories

**Outputs:**
* `outputs/Fig9_Duration_Grid.png` - 29-panel duration trajectory grid
* `outputs/trajectories_data.json` - Scatter data with subtype-specific regression parameters
* `outputs/trajectories.html` - Interactive trajectory explorer with subtype filtering

**Script:** `generate_figure.py`

**Biological Context:** Reveals prognostic markers – amino acids with steep positive slopes in progressive subtypes (SPMS/PPMS) indicate accelerated metabolic dysregulation over time, while stable trajectories in RRMS suggest relapse-remission metabolic resilience.

---

### **Node 10: Clinical Mimicry Test (RRMS vs MG)**

**Purpose:** Distinguishes early-stage MS (RRMS) from MG in the diagnostic gray zone where symptoms overlap.

**Input:**
* `data_standardized.pkl` (filtered for RRMS and MG only)

**Processing:**
* **Figure 10**: Total AA load comparison (RRMS vs MG)
* **Figure 11**: 29-panel amino acid grid (RRMS vs MG)
* Mann-Whitney U tests for each amino acid

**Outputs:**
* `outputs/Fig10_TotalAA_RRMS_MG.png` - Global metabolic load comparison
* `outputs/Fig11_RRMS_vs_MG_Grid.png` - Full amino acid profile comparison
* `outputs/mimicry_data.json` - Complete boxplot data for early-stage differentiation
* `outputs/mimicry.html` - Interactive diagnostic decision-support tool

**Script:** `generate_figures.py`

**Biological Context:** Identifies **3-Methylhistidine (3-MHIS)** as a critical early discriminator – elevated in MG due to muscle protein breakdown at neuromuscular junction, normal in RRMS where pathology is confined to CNS myelin. Supports differential diagnosis when clinical presentation is ambiguous.

---

### **Node 11: Pathway Coherence Analysis**

**Purpose:** Contrasts metabolic network organization between MS and MG to reveal disease-specific metabolic architectures.

**Input:**
* `data_standardized.pkl`
* `aa_cols.txt`

**Processing:**
* Side-by-side triangular correlation heatmaps (MS vs MG)
* Identical amino acid ordering for direct visual comparison
* Pearson correlation with hierarchical clustering

**Outputs:**
* `outputs/Fig12_Split_Corr.png` - Split-panel correlation comparison (MS | MG)
* `outputs/coherence_data.json` - Dual correlation matrices with cluster metrics
* `outputs/coherence.html` - Interactive side-by-side network explorer

**Script:** `generate_figure.py`

**Biological Context:** Demonstrates fundamental metabolic architectural differences – MS shows highly fragmented correlation structure (chaotic central demyelination effects), while MG retains more organized peripheral metabolic networks (localized neuromuscular pathology). Validates that MS and MG, despite shared autoimmune mechanisms, exhibit distinct systemic metabolic dysregulation patterns.

---

### Workflow Execution

**Run Complete Workflow:**
```bash
silva .
```
The `silva` command discovers all `job.toml` files, resolves dependencies, and executes nodes in the correct order. Node 01 runs first, then Nodes 02-11 execute in parallel inside Docker containers.

**How Silva Works:**
1. Reads `.chiral/workflow.toml` in the root directory for workflow-level configuration and dependencies
2. Scans each node's `.chiral/job.toml` for inputs, outputs, and container settings
3. Builds a dependency graph based on the `[dependencies]` section in `workflow.toml`
4. Pulls or reuses the Docker image (`chiral.sakuracr.jp/proteomics:2025_12_31`)
5. Executes each node's `run.sh` script inside the container
6. Copies input files between nodes as specified in `job.toml`
7. Collects outputs to a timestamped folder (e.g., `C:\Windows\TEMP\silva-2026-01-03-...`)
8. Cleans up containers after workflow completion

**Node Configuration (`.chiral/job.toml`):**
Each node contains a `.chiral/job.toml` file that specifies:
* `name` - Human-readable node name
* `description` - What the node does
* `inputs` - Files required from other nodes
* `outputs` - Files produced by this node
* `container.image` - Docker image to use
* `scripts.run` - Shell script to execute (typically `run.sh`)

**Output Structure:**
Each node generates 3 output types in its `outputs/` directory:
1. **PNG files** - Publication-ready figures (300 DPI, matplotlib/seaborn)
2. **JSON files** - Structured data (plot coordinates, statistics, metadata)
3. **HTML files** - Interactive visualizations (Plotly.js, responsive design)

**Total Outputs:** 38 files per complete workflow execution
* Node 01: 2 files (data + metadata)
* Nodes 02-11: 36 files (18 PNG + 10 JSON + 10 HTML)

**Modular HTML Generation Infrastructure:**
Each analysis node (02-11) contains its own `html_generator.py` module with:
* **Self-Contained Design**: No shared dependencies between nodes
* **Node-Specific Functions**: 
  - Node 02: `generate_demographics_html()` with `renderTable()`
  - Node 03: `generate_pathology_html()` with 5-column grid layout for 29 amino acids
  - Nodes 04-11: Dual-mode generators with `renderGridPlot()` and `renderBoxPlot()`
* **Configuration Objects**: Each generator includes its own `CONFIG` object with color schemes, layout parameters, and plot styling
* **Plotly.js Integration**: CDN v2.27.0 for interactive plotting
* **Features**: Zoom, pan, hover tooltips, data export, responsive design
* **Debugging Benefits**: Smaller files (~118 lines vs 572-line shared template) for easier troubleshooting

**How to Modify Visualizations:**
1. Navigate to node directory: `cd 03_MS_Pathology_Overview/`
2. Edit `html_generator.py` to customize rendering functions or CONFIG settings
3. Regenerate HTML: `bash run.sh`
4. Changes affect only this node - no impact on other nodes

---

## Part 3: Quick Start Guide

### Prerequisites

**Required Software:**
* **Docker**: Must be installed and running (containers execute the analysis)
* **Silva**: Workflow orchestration tool (`silva` command must be available in PATH)

**Input Data:**
The database file should be present in Node 01 directory:
```
01_Data_Ingestion_and_Preprocessing/database-multiple-sclerosis-myasthenia.csv
```

**Note:** Python packages are pre-installed in the Docker image (`chiral.sakuracr.jp/proteomics:2025_12_31`). No local Python environment setup is required.

### Execution Commands

**Run Complete Workflow:**
```bash
silva .
```
This command executes all 11 nodes with proper dependency ordering. Outputs are collected to a timestamped folder.

**Run Individual Nodes (for development/debugging):**
```bash
# Navigate to node directory and run the script directly
cd 03_MS_Pathology_Overview
docker run -v $(pwd):/workspace chiral.sakuracr.jp/proteomics:2025_12_31 bash run.sh
```

### Output Summary

Each node generates files in its `outputs/` directory:

| Node | PNG Outputs | JSON/HTML Outputs |
|------|-------------|-------------------|
| 01 | data_standardized.pkl, aa_cols.txt | N/A |
| 02 | Table1_Demographics.png, Table2_MS_MG_Demographics.png | demographics_data.json, demographics.html |
| 03 | Fig1A_MS_vs_Control.png, Fig1B_MS_Subtypes.png | pathology_data.json, pathology.html |
| 04 | Fig2A_Age_Grid.png, Fig2B_Duration_Grid.png, Fig2C_EDSS_Grid.png | confounders_data.json, confounders.html |
| 05 | Fig3_MS_MG_vs_Control.png | autoimmune_data.json, autoimmune.html |
| 06 | Fig4_Specific_Diffs.png, Fig5_Female_Specific.png | biomarkers_data.json, biomarkers.html |
| 07 | Fig6_Corr_MS_MG.png, Fig7_Corr_Control.png | clustering_data.json, clustering.html |
| 08 | Fig8_TotalAA_Type.png | metabolic_load_data.json, metabolic_load.html |
| 09 | Fig9_Duration_Grid.png | trajectories_data.json, trajectories.html |
| 10 | Fig10_TotalAA_RRMS_MG.png, Fig11_RRMS_vs_MG_Grid.png | mimicry_data.json, mimicry.html |
| 11 | Fig12_Split_Corr.png | coherence_data.json, coherence.html |

**Total Outputs:** 38 files (18 PNG + 10 JSON + 10 HTML)

**Execution Time:** ~2-3 minutes for complete workflow (with Docker image cached)

### Interactive Visualizations

**View HTML Outputs:**
1. Navigate to node output folder: `cd XX_Node_Name/outputs/`
2. Open HTML file in browser: `[node_name].html`
3. Features: Zoom, pan, hover tooltips, data export

**Using Python HTTP Server:**
```bash
cd XX_Node_Name/outputs/
python -m http.server 8000
# Open browser to http://localhost:8000/[node_name].html
```

### Troubleshooting

**Docker not running:**
```bash
# Ensure Docker Desktop is running before executing silva
docker info
```

**Database file not found:**
```bash
# Check if file exists
ls 01_Data_Ingestion_and_Preprocessing/database-multiple-sclerosis-myasthenia.csv
```

**Silva command not found:**
```bash
# Ensure silva is installed and in PATH
# Contact your administrator for installation instructions
```

**Node fails - Check container logs:**
```bash
# Silva outputs detailed logs during execution
# Look for error messages in the terminal output
# Common issues: missing input files, Python errors in scripts
```

**Clean outputs (start fresh):**
```bash
# Remove all outputs
rm -rf */outputs/*

# Remove specific node
rm -rf 03_MS_Pathology_Overview/outputs/*
```

**Modify HTML visualizations:**
```bash
# Each node has its own html_generator.py - edit to customize
cd 03_MS_Pathology_Overview
nano html_generator.py  # or use your preferred editor

# Modify CONFIG object (colors, layout) or rendering functions
# Re-run the workflow to regenerate
cd .. && silva .
```

**HTML Generator Structure:**
* **Node 02**: Custom `renderTable()` for demographics
* **Node 03**: Custom `renderGridPlot()` with 5-column layout for 29 amino acids
* **Nodes 04-11**: Generic dual-mode generators with `renderGridPlot()` and `renderBoxPlot()`
* Each includes its own CONFIG object with color schemes and styling

### Workflow Structure

```
wf-abdo-proteomics/
├── .chiral/
│   └── workflow.toml              # Main workflow configuration (dependencies)
├── README.md                      # This documentation file
│
├── 01_Data_Ingestion_and_Preprocessing/
│   ├── .chiral/
│   │   └── job.toml               # Node configuration for silva
│   ├── database-multiple-sclerosis-myasthenia.csv  # Input dataset
│   ├── load_data.py               # Data loading and standardization script
│   └── run.sh                     # Execution script (called by silva)
│
├── 02_Cohort_Demographics/
│   ├── .chiral/
│   │   └── job.toml               # Node configuration (inputs, outputs, image)
│   ├── generate_tables.py         # Demographics analysis script
│   ├── html_generator.py          # Node-specific HTML generator
│   └── run.sh                     # Execution script
│
└── [03-11]_*/                     # Additional analysis nodes
    ├── .chiral/
    │   └── job.toml               # Node configuration for silva
    ├── generate_*.py              # Analysis script (generate_figure.py or generate_figures.py)
    ├── html_generator.py          # Node-specific HTML generator (~12 KB each)
    └── run.sh                     # Execution script
```

**Modular Architecture Notes:**
* Each node (02-11) has its own `html_generator.py` module - no shared dependencies
* Benefits: Easier debugging (~118 lines per file vs 572-line shared template)
* Isolated development: Changes to one node don't affect others
* Self-contained: Each generator includes its own CONFIG, rendering functions, and Plotly.js integration

### Data Flow

```
Node 01 (Preprocessing)
    ↓
    ├─→ Node 02 (Demographics)
    ├─→ Node 03 (MS Pathology)
    ├─→ Node 04 (Confounders)
    ├─→ Node 05 (MS vs MG)
    ├─→ Node 06 (Biomarkers)
    ├─→ Node 07 (Network)
    ├─→ Node 08 (Metabolic Load)
    ├─→ Node 09 (Trajectories)
    ├─→ Node 10 (RRMS vs MG)
    └─→ Node 11 (Coherence)
```

**Execution Model:** Silva orchestrates the workflow - Node 01 runs first, then Nodes 02-11 execute in parallel inside Docker containers. Run with `silva .` from the workflow directory.

---

## Part 4: Data Sources and Acknowledgments

### Origin
The datasets utilized in this study are educational subsets curated from previously published clinical and omics research regarding Multiple Sclerosis (MS) and Myasthenia Gravis (MG). These subsets were specifically filtered and prepared by **Prof. Emilia Daghir-Wojtkowiak** to serve as educational material.

> **Permission**: Formal permission has been granted by Prof. Daghir-Wojtkowiak to utilize this modified dataset for the scope of this volunteer work at Chiral.

### Authorship
*   **Data Curation**: Faculty (Prof. Emilia Daghir-Wojtkowiak)
*   **Analytical Implementation**: Abdelrahman Mohamed Taha MAHMOUD (Code development, pipeline design, statistical analysis)

### References
The original data sources correspond to the following publications:

1.  **For MS Outcomes**:
**Serum amino acid profiling in differentiating clinical outcomes of multiple sclerosis.**
*Neurologia i Neurochirurgia Polska*, 57(5), 414–422.
[https://doi.org/10.5603/PJNNS.a2023.0054](https://doi.org/10.5603/PJNNS.a2023.0054)

2.  **For MS vs. MG Comparison**: 
**Comparative Analysis of Serum Amino Acid Profiles in Patients with Myasthenia Gravis and Multiple Sclerosis**
*J. Clin. Med. 2024, 13(14), 4083;*
[https://doi.org/10.3390/jcm13144083](https://doi.org/10.3390/jcm13144083)
