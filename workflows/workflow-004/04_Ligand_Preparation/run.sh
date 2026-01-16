#!/bin/bash
set -e
echo "Starting Node 04: Ligand Preparation"

mkdir -p inputs
cp *.sdf inputs/ 2>/dev/null || true

python prepare_ligands_for_docking.py
python generate_ligand_refinement_report.py

echo "Node 04 completed"