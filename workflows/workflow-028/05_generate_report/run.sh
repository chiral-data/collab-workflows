#!/bin/bash
set -e

echo "Starting Node 05: Generate Report"

cp inputs/plddt_*.npz            . 2>/dev/null || true
cp inputs/pae_*.npz              . 2>/dev/null || true
cp inputs/affinity_*.json        . 2>/dev/null || true
cp inputs/confidence_*.json      . 2>/dev/null || true
cp inputs/input_summary.json     . 2>/dev/null || true
cp inputs/predictions.csv        . 2>/dev/null || true
cp inputs/selected_model_id.txt  . 2>/dev/null || true
cp inputs/pocket_qc.json         . 2>/dev/null || true
cp inputs/docked_poses.sdf       . 2>/dev/null || true
cp inputs/docking_summary.json   . 2>/dev/null || true

python3 generate_report.py

echo "Node 05 completed"
