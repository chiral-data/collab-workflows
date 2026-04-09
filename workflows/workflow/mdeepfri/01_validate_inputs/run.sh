#!/bin/bash
set -e
echo "Starting Node 01: Validate Inputs"
python validate.py \
    --min-length "${PARAM_MIN_LENGTH:-10}" \
    --max-length "${PARAM_MAX_LENGTH:-5000}"
echo "Node 01 completed"
