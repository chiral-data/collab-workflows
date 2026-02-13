#!/bin/bash
set -e

echo "Starting Node 01: Sequence Upload"

python validate_input.py \
    --input "${PARAM_INPUT_FILE}" \
    --output validated_input.yaml

echo "Node 01 completed"
