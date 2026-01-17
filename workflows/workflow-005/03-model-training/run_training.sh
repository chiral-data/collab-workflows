#!/bin/bash
set -e
echo "Starting Node 03: Model Training"
mkdir -p outputs

python run_training.py
python generate_training_report.py

echo "Node 03 completed"
