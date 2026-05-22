#!/bin/bash
set -e

# Node 5: Interface Analysis
# Analyzes the best-scoring docking pose to identify antibody-antigen binding interface

echo "=========================================="
echo "  Node 5: Interface Analysis"
echo "=========================================="

# Activate conda environment if available
if command -v conda &> /dev/null; then
    CONDA_PATH=$(conda info --base)
    source "$CONDA_PATH/etc/profile.d/conda.sh"
    conda activate diffdock_abag
fi

# Get the directory where this script is located
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

ROOT_DIR="$(dirname "$PWD")"

# Define paths
RECEPTOR_PDB="$(find inputs -type f -name processed_antibody.pdb | head -n 1)"
BEST_POSE_PDB="$(find inputs -type f -name rank1.pdb | head -n 1)"
CONFIDENCE_JSON="$(find inputs -type f -name confidence_scores.json | head -n 1)"

if [ -z "$RECEPTOR_PDB" ] && [ -f "$ROOT_DIR/2_diffdock_prep/outputs/processed_antibody.pdb" ]; then
    RECEPTOR_PDB="$ROOT_DIR/2_diffdock_prep/outputs/processed_antibody.pdb"
fi
if [ -z "$BEST_POSE_PDB" ] && [ -f "$ROOT_DIR/4_diffdock_inference/outputs/rank1.pdb" ]; then
    BEST_POSE_PDB="$ROOT_DIR/4_diffdock_inference/outputs/rank1.pdb"
fi
if [ -z "$CONFIDENCE_JSON" ] && [ -f "$ROOT_DIR/4_diffdock_inference/outputs/confidence_scores.json" ]; then
    CONFIDENCE_JSON="$ROOT_DIR/4_diffdock_inference/outputs/confidence_scores.json"
fi

if [ -z "$RECEPTOR_PDB" ] || [ -z "$BEST_POSE_PDB" ] || [ -z "$CONFIDENCE_JSON" ]; then
    echo "Error: missing required inputs for Node 5"
    exit 1
fi
OUTPUT_DIR="outputs"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Run science analysis
echo ""
echo "Running interface analysis..."
python 9_analysis_science.py \
    --receptor_pdb "$RECEPTOR_PDB" \
    --best_pose_pdb "$BEST_POSE_PDB" \
    --confidence_json "$CONFIDENCE_JSON" \
    --output_dir "$OUTPUT_DIR"

# Generate HTML report
echo ""
echo "Generating interactive HTML report..."
python 10_analysis_html.py \
    --data_json "$OUTPUT_DIR/data.json" \
    --contact_residues_json "$OUTPUT_DIR/contact_residues.json" \
    --interface_analysis_txt "$OUTPUT_DIR/interface_analysis.txt" \
    --output_html "$OUTPUT_DIR/report.html"

echo ""
echo "=========================================="
echo "  Node 5 Complete!"
echo "=========================================="
echo "Outputs saved to: $OUTPUT_DIR"
echo "  - interface_analysis.txt (detailed contacts)"
echo "  - contact_residues.json (machine-readable)"
echo "  - final_complex.pdb (combined structure)"
echo "  - data.json (summary statistics)"
echo "  - report.html (interactive dashboard)"
echo "=========================================="