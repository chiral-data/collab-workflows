#!/bin/bash
set -e

echo "Starting Node 05: Generate Report"

cp inputs/vina_screening_scores.csv . 2>/dev/null || true
cp inputs/gnina_screening_scores.csv . 2>/dev/null || true
cp inputs/vina_docking_report.json . 2>/dev/null || true
cp inputs/gnina_docking_report.json . 2>/dev/null || true
cp inputs/vina_screening_poses.pdbqt . 2>/dev/null || true
cp inputs/gnina_screening_poses.sdf . 2>/dev/null || true

python3 generate_report.py

echo "Node 05 completed"
