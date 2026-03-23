#!/bin/bash
set -e

echo "Starting Node 01: Sequence Upload"

# Copy input files from silva's inputs/ directory to working directory
cp inputs/* . 2>/dev/null || true

# Find YAML input from input_files/
INPUT_FILE=$(ls *.yaml 2>/dev/null | head -1)

if [ -z "$INPUT_FILE" ]; then
    echo "Error: No input YAML file found in inputs/"
    exit 1
fi

python validate_input.py \
    --input "$INPUT_FILE" \
    --output validated_input.yaml

echo "Node 01 completed"
