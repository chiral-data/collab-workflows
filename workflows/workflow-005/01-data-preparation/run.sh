#!/bin/bash
set -e
echo "Starting Node 01: Data Preparation"
mkdir -p outputs

python run_data_prep.py
python generate_data_prep_report.py

echo "Node 01 completed"
