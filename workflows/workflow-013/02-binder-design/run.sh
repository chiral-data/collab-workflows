#!/bin/bash
set -e

echo "Starting Node 02: BoltzGen Binder Design"

# Copy input files from silva's inputs/ directory to working directory
cp inputs/* . 2>/dev/null || true

# Find design spec from previous node
DESIGN_SPEC=$(ls design_spec.yaml 2>/dev/null | head -1)

if [ -z "$DESIGN_SPEC" ]; then
    echo "Error: No design_spec.yaml found from previous step"
    exit 1
fi

echo "Design spec: $DESIGN_SPEC"

python run_design.py \
    --design-spec "$DESIGN_SPEC" \
    --protocol "${PARAM_PROTOCOL:-protein-anything}" \
    --num-designs "${PARAM_NUM_DESIGNS:-50}" \
    --budget "${PARAM_BUDGET:-10}" \
    --output-dir "."

echo "Node 02 completed"
