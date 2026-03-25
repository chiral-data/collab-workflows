#!/bin/bash
set -e
echo "Starting Node 01: Validate SMILE Inputs"

python validate.py \
    --input "${PARAM_INPUT_FILE}" \
    --output standardized_molecules.csv

echo "Node 01 completed"
