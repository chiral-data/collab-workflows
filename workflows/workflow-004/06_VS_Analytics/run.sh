#!/bin/bash
set -e
echo "Starting Node 06: VS Analytics"

mkdir -p inputs
cp *.pdbqt inputs/ 2>/dev/null || true
cp *.json inputs/ 2>/dev/null || true

python run_autodock_vina.py
python generate_scientific_dashboard.py

echo "Node 06 completed"