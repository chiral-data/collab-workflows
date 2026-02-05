#!/bin/bash
set -e
echo "Starting Node 01: Data Preparation"
mkdir -p outputs

python 1_data_prep.py
python 2_data_prep_HTML.py

echo "Node 01 completed"
