#!/bin/bash
set -e
echo "Starting Node 03: Receptor Preparation"

# Copy inputs from upstream (Silva copies them to root, script expects them)
mkdir -p inputs
cp *.pdb inputs/ 2>/dev/null || true

python prepare_receptor_for_docking.py
python generate_refinement_report.py

echo "Node 03 completed"