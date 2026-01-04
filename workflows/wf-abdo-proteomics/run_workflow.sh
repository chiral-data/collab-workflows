#!/bin/bash
# Master Workflow Execution Script
# MS and MG Amino Acid Profiling Analysis

set -e

echo "=========================================="
echo "MS & MG Amino Acid Analysis Workflow"
echo "=========================================="
echo ""

# Check if database file exists
if [ ! -f "01_Data_Ingestion_and_Preprocessing/database-multiple-sclerosis-myasthenia.csv" ]; then
    echo "ERROR: Database file not found!"
    echo "Please place 'database-multiple-sclerosis-myasthenia.csv' in the 01_Data_Ingestion_and_Preprocessing folder"
    exit 1
fi

# Node 01: Data Ingestion and Preprocessing
echo "[1/11] Running Node 01: Data Ingestion and Preprocessing..."
cd 01_Data_Ingestion_and_Preprocessing
bash run.sh
cd ..
echo "✓ Node 01 completed"
echo ""

# Node 02: Cohort Demographics
echo "[2/11] Running Node 02: Cohort Demographics..."
cd 02_Cohort_Demographics
bash run.sh
cd ..
echo "✓ Node 02 completed"
echo ""

# Nodes 03-11 can run in parallel, but we'll run them sequentially for safety
# Node 03: MS Pathology Overview
echo "[3/11] Running Node 03: MS Pathology Overview..."
cd 03_MS_Pathology_Overview
bash run.sh
cd ..
echo "✓ Node 03 completed"
echo ""

# Node 04: Clinical Confounders and Validation
echo "[4/11] Running Node 04: Clinical Confounders and Validation..."
cd 04_Clinical_Confounders_and_Validation
bash run.sh
cd ..
echo "✓ Node 04 completed"
echo ""

# Node 05: MS vs MG Autoimmune Comparison
echo "[5/11] Running Node 05: MS vs MG Autoimmune Comparison..."
cd 05_MS_vs_MG_Autoimmune_Comparison
bash run.sh
cd ..
echo "✓ Node 05 completed"
echo ""

# Node 06: Differential Diagnosis Biomarkers
echo "[6/11] Running Node 06: Differential Diagnosis Biomarkers..."
cd 06_Differential_Diagnosis_Biomarkers
bash run.sh
cd ..
echo "✓ Node 06 completed"
echo ""

# Node 07: Metabolic Network Clustering
echo "[7/11] Running Node 07: Metabolic Network Clustering..."
cd 07_Metabolic_Network_Clustering
bash run.sh
cd ..
echo "✓ Node 07 completed"
echo ""

# Node 08: Global Metabolic Load Analysis
echo "[8/11] Running Node 08: Global Metabolic Load Analysis..."
cd 08_Global_Metabolic_Load_Analysis
bash run.sh
cd ..
echo "✓ Node 08 completed"
echo ""

# Node 09: Subtype Trajectories
echo "[9/11] Running Node 09: Subtype Trajectories..."
cd 09_Subtype_Trajectories
bash run.sh
cd ..
echo "✓ Node 09 completed"
echo ""

# Node 10: Clinical Mimicry Test
echo "[10/11] Running Node 10: Clinical Mimicry Test (RRMS vs MG)..."
cd 10_Clinical_Mimicry_Test_RRMS_vs_MG
bash run.sh
cd ..
echo "✓ Node 10 completed"
echo ""

# Node 11: Pathway Coherence Analysis
echo "[11/11] Running Node 11: Pathway Coherence Analysis..."
cd 11_Pathway_Coherence_Analysis
bash run.sh
cd ..
echo "✓ Node 11 completed"
echo ""

echo "=========================================="
echo "✓ ALL NODES COMPLETED SUCCESSFULLY!"
echo "=========================================="
echo ""
echo "Output files have been generated in each node's 'output' folder:"
echo "  - Node 01: Preprocessed data"
echo "  - Node 02: Tables 1 & 2"
echo "  - Node 03: Figures 1A & 1B"
echo "  - Node 04: Figures 2A, 2B, 2C"
echo "  - Node 05: Figure 3"
echo "  - Node 06: Figures 4 & 5"
echo "  - Node 07: Figures 6 & 7"
echo "  - Node 08: Figure 8"
echo "  - Node 09: Figure 9"
echo "  - Node 10: Figures 10 & 11"
echo "  - Node 11: Figure 12"
echo ""
