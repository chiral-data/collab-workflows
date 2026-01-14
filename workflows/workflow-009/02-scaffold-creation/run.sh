#!/bin/bash
set -e

# Find the ligand SMILES file from previous step (contains atom map numbers)
LIGAND_FILE=$(ls *.smi 2>/dev/null | head -1)

if [ -z "$LIGAND_FILE" ]; then
    echo "Error: No SMILES file found from previous step"
    exit 1
fi

workflow-run python create_scaffold.py \
    --ligand "${LIGAND_FILE}" \
    --attachment-id "${PARAM_ATTACHMENT_ID}" \
    --output scaffold.pkl

# Generate visualization HTML
workflow-run python visualize.py scaffold.pkl scaffold_viz.html
