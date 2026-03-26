#!/bin/bash
set -e

echo "Starting Node 01: Validate Inputs (Boltz + Chai)"

# Copy input files to working directory
cp ../input/prot.yaml . 2>/dev/null || true

# Detect input files
YAML_FILE=$(ls *.yaml *.yml 2>/dev/null | head -1)
FASTA_FILE=$(ls *.fasta *.fa 2>/dev/null | head -1)

if [ -z "$YAML_FILE" ] && [ -z "$FASTA_FILE" ]; then
    echo "Error: No YAML or FASTA input file found in inputs/"
    exit 1
fi

# Build args dynamically based on what's present
ARGS=""
if [ -n "$YAML_FILE" ]; then
    echo "Found YAML (Boltz): $YAML_FILE"
    ARGS="$ARGS --yaml $YAML_FILE"
fi
if [ -n "$FASTA_FILE" ]; then
    echo "Found FASTA (Chai): $FASTA_FILE"
    ARGS="$ARGS --fasta $FASTA_FILE"
fi

python3 validate_input.py $ARGS

echo "Node 01 completed"
