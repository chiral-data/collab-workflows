#!/bin/bash
set -e

workflow-run python validate_ligand.py \
    --input "${PARAM_LIGAND_FILE}" \
    --output ligand.smi

# Generate visualization HTML
workflow-run python visualize.py ligand.smi ligand_viz.html
