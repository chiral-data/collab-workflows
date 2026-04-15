#!/bin/bash
set -e
echo "Starting Node 03: Generate Visualization"
python generate_report.py \
    --min-score "${PARAM_MIN_SCORE:-0.3}" \
    --top-n "${PARAM_TOP_N_TERMS:-10}"
echo "Node 03 completed"
