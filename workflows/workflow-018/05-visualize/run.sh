#!/bin/bash
set -e

echo "Starting Node 05: Visualization"

cp inputs/* . 2>/dev/null || true

# Verify inputs
if [ ! -f boltz_summary.json ]; then
    echo "Error: boltz_summary.json not found — did node 03 run?"
    exit 1
fi
if [ ! -f chai_summary.json ]; then
    echo "Error: chai_summary.json not found — did node 04 run?"
    exit 1
fi

# Detect ground truth structure from node 00-download (optional)
REFERENCE_PDB=$(ls *_reference.pdb 2>/dev/null | head -1)
REF_ARG=""
if [ -n "$REFERENCE_PDB" ]; then
    echo "Found reference structure: $REFERENCE_PDB"
    REF_ARG="--reference-pdb $REFERENCE_PDB"
fi

python3 generate_report.py \
    --boltz-summary boltz_summary.json \
    --chai-summary  chai_summary.json \
    --output        comparison_report.html \
    $REF_ARG

echo "Node 05 completed"
