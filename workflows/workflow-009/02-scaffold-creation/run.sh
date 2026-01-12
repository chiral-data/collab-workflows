#!/bin/bash
set -e

# Find the ligand file from previous step
LIGAND_FILE=$(ls *.sdf 2>/dev/null | head -1)

if [ -z "$LIGAND_FILE" ]; then
    echo "Error: No SDF file found from previous step"
    exit 1
fi

workflow-run python create_scaffold.py \
    --ligand "${LIGAND_FILE}" \
    --attachment-id "${PARAM_ATTACHMENT_ID:-27}" \
    --output scaffold.pkl

# Generate visualization HTML
workflow-run python visualize.py scaffold.pkl scaffold_viz.html
