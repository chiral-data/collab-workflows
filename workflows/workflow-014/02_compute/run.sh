#!/bin/bash
set -e
echo "Starting Node 02: Compute ADMET Predictions"
mkdir -p outputs

python compute_admet.py

echo "Node 02 completed"
