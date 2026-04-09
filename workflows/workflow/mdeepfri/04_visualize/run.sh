#!/bin/bash
set -e
echo "Starting Node 04: Generate Visualization"
python generate_report.py \
    --min-score "${PARAM_MIN_SCORE:-0.3}" \
    --top-n "${PARAM_TOP_N_TERMS:-10}"
echo "Node 04 completed"
