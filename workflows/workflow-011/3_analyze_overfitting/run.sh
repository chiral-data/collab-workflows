#!/bin/bash
set -e
echo "Starting Node 03: Analyze Overfitting"
mkdir -p outputs

python 5_analyze_overfitting.py
python 6_analyze_overfitting_HTML.py

echo "Node 03 completed"
