#!/bin/bash
set -e

workflow-run python validate_ligand.py \
    --input "inputs/${PARAM_LIGAND_FILE}" \
    --output outputs/ligand.smi

# Generate visualization HTML
workflow-run python visualize.py outputs/ligand.smi outputs/ligand_viz.html
