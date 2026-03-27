#!/bin/bash
set -e
echo "Starting Node 04: Prediction"
mkdir -p outputs

python 7_predict_from_csv.py
python 8_predict_from_csv_HTML.py

echo "Node 04 completed"

