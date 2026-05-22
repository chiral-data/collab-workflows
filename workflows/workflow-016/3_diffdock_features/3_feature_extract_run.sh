#!/bin/bash
set -e
echo "=== Node 3: Feature extraction ==="

ROOT_DIR="$(dirname "$PWD")"

PROC_AB="$(find inputs -type f -name processed_antibody.pdb | head -n 1)"
PROC_AG="$(find inputs -type f -name processed_antigen.pdb | head -n 1)"

if [ -z "$PROC_AB" ] && [ -f "$ROOT_DIR/2_diffdock_prep/outputs/processed_antibody.pdb" ]; then
    PROC_AB="$ROOT_DIR/2_diffdock_prep/outputs/processed_antibody.pdb"
fi
if [ -z "$PROC_AG" ] && [ -f "$ROOT_DIR/2_diffdock_prep/outputs/processed_antigen.pdb" ]; then
    PROC_AG="$ROOT_DIR/2_diffdock_prep/outputs/processed_antigen.pdb"
fi

if [ -z "$PROC_AB" ] || [ -z "$PROC_AG" ]; then
    echo "Error: missing processed inputs for Node 3"
    exit 1
fi

python 5_feature_extract_science.py \
        --antibody_pdb "$PROC_AB" \
        --antigen_pdb "$PROC_AG" \
    --output_dir "outputs"
python 6_feature_extract_html.py --data_json "outputs/data.json" --output_html "outputs/report.html"
echo "Node 3 finished"