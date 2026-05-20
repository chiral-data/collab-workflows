#!/bin/bash
set -e

echo "Starting Node 02: Preprocessing"

cp inputs/* . 2>/dev/null || true

# Find the validated YAML from node 01
YAML_FILE=$(ls boltz_input.yaml *.yaml 2>/dev/null | head -1)

if [ -z "$YAML_FILE" ]; then
    echo "Error: No YAML input file found"
    exit 1
fi

echo "Input file: $YAML_FILE"

python3 preprocess.py \
    --input "$YAML_FILE" \
    --boltz-output boltz_input.yaml \
    --chai-output chai_input.fasta

echo "Node 02 completed"
