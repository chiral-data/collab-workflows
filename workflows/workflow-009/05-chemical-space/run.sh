#!/bin/bash
set -e

# Find inputs from previous steps
SCAFFOLD_FILE=$(ls inputs/scaffold.pkl 2>/dev/null | head -1)
PROTEIN_FILE=$(ls inputs/rec_final.pdb 2>/dev/null | head -1)

if [ -z "$SCAFFOLD_FILE" ]; then
    echo "Error: scaffold.pkl not found from previous step"
    exit 1
fi

if [ -z "$PROTEIN_FILE" ]; then
    echo "Error: rec_final.pdb not found from previous step"
    exit 1
fi

workflow-run python create_chemspace.py \
    --scaffold "${SCAFFOLD_FILE}" \
    --protein "${PROTEIN_FILE}" \
    --num-linkers "${PARAM_NUM_LINKERS:-5}" \
    --num-rgroups "${PARAM_NUM_RGROUPS:-5}" \
    --output outputs/chemspace.pkl

# Copy protein file to outputs for downstream nodes
cp "${PROTEIN_FILE}" outputs/rec_final.pdb

# Generate visualization HTML
workflow-run python visualize.py outputs/chemspace.pkl outputs/chemspace_viz.html
