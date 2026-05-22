#!/bin/bash
set -e
echo "=== Node 2: Structure preparation ==="

ROOT_DIR="$(dirname "$PWD")"

ANTIBODY_PDB="$(find inputs -type f -name antibody.pdb | head -n 1)"
ANTIGEN_PDB="$(find inputs -type f -name antigen.pdb | head -n 1)"

if [ -z "$ANTIBODY_PDB" ] && [ -f "$ROOT_DIR/1_complex_splitting/outputs/antibody.pdb" ]; then
    ANTIBODY_PDB="$ROOT_DIR/1_complex_splitting/outputs/antibody.pdb"
fi
if [ -z "$ANTIGEN_PDB" ] && [ -f "$ROOT_DIR/1_complex_splitting/outputs/antigen.pdb" ]; then
    ANTIGEN_PDB="$ROOT_DIR/1_complex_splitting/outputs/antigen.pdb"
fi

if [ -z "$ANTIBODY_PDB" ] || [ -z "$ANTIGEN_PDB" ]; then
    echo "Error: missing antibody/antigen inputs for Node 2"
    exit 1
fi

python 3_structure_prep_science.py \
        --antibody_pdb "$ANTIBODY_PDB" \
        --antigen_pdb "$ANTIGEN_PDB" \
    --output_dir "outputs"
python 4_structure_prep_html.py --data_json "outputs/data.json" --output_html "outputs/report.html"
echo "Node 2 finished"