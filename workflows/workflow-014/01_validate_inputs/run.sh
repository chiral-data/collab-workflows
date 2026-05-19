#!/bin/bash
set -e
echo "Starting Node 01: Validate SMILES Inputs"
python validate.py \
    --input "${PARAM_INPUT_FILE:-./drugbank_approved.csv}" \
    --smiles-column "${PARAM_SMILES_COLUMN:-smiles}"

echo "Node 01 completed"
