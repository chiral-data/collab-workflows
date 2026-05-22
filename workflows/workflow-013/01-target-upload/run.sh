#!/bin/bash
set -e

echo "Starting Node 01: Target Upload"

python validate_target.py \
    --design-spec "${PARAM_DESIGN_SPEC}" \
    --target-structure "${PARAM_TARGET_STRUCTURE:-}" \
    --output-dir "."

echo "Node 01 completed"
