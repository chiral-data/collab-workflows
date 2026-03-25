#!/bin/bash
set -e
echo "Starting Node 01: Validate SMILES Inputs"
mkdir -p outputs

python validate.py \
    --input "${PARAM_INPUT_FILE}" \
    --smiles-column "${PARAM_SMILES_COLUMN:-smiles}"

echo "Node 01 completed"
