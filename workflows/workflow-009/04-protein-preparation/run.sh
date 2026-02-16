#!/bin/bash
set -e

# Find the protein file from previous step
PROTEIN_FILE=$(ls inputs/*.pdb 2>/dev/null | head -1)

if [ -z "$PROTEIN_FILE" ]; then
    echo "Error: No PDB file found from previous step"
    exit 1
fi

workflow-run python prepare_protein.py \
    --input "${PROTEIN_FILE}" \
    --output outputs/rec_final.pdb
