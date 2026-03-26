#!/bin/bash
set -e

echo "Starting Node 05: Visualization"

cp inputs/* . 2>/dev/null || true

# Pull summaries from prediction nodes
cp ../03-predict-boltz/boltz_summary.json . 2>/dev/null || true
cp ../04-predict-chai/chai_summary.json . 2>/dev/null || true

# Verify inputs
if [ ! -f boltz_summary.json ]; then
    echo "Error: boltz_summary.json not found — did node 03 run?"
    exit 1
fi
if [ ! -f chai_summary.json ]; then
    echo "Error: chai_summary.json not found — did node 04 run?"
    exit 1
fi

python3 generate_report.py \
    --boltz-summary boltz_summary.json \
    --chai-summary  chai_summary.json \
    --output        comparison_report.html

echo "Node 05 completed"
