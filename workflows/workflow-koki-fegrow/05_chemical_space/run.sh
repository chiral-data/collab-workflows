#!/bin/bash
set -e

# Find inputs from previous steps
SCAFFOLD_FILE=$(ls scaffold.pkl 2>/dev/null | head -1)
PROTEIN_FILE=$(ls rec_final.pdb 2>/dev/null | head -1)

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
    --num-linkers "${PARAM_NUM_LINKERS:-10}" \
    --num-rgroups "${PARAM_NUM_RGROUPS:-10}" \
    --output chemspace.pkl

# Generate visualization HTML
workflow-run python visualize.py chemspace.pkl chemspace_viz.html
