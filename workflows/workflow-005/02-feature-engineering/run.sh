#!/bin/bash
set -e
echo "Starting Node 02: Feature Engineering"
mkdir -p outputs

python run_feature_eng.py
python generate_feature_eng_report.py

echo "Node 02 completed"
