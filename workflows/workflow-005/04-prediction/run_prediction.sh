#!/bin/bash
set -e
echo "Starting Node 04: Prediction"
mkdir -p outputs

python run_prediction.py
python generate_prediction_report.py

echo "Node 04 completed"
