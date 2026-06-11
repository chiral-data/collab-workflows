#!/bin/bash
set -e

echo "Starting Node 04a: AutoDock Vina Screening"

cp inputs/receptor.pdbqt . 2>/dev/null || true
cp inputs/optimized_screening_library.pdbqt . 2>/dev/null || true
cp inputs/pocket_config.txt . 2>/dev/null || true
cp inputs/crystal_ligand.pdb . 2>/dev/null || true

python3 run_vina.py \
    --receptor "${PARAM_RECEPTOR:-receptor.pdbqt}" \
    --library "${PARAM_LIBRARY:-optimized_screening_library.pdbqt}" \
    --pocket-config "${PARAM_POCKET_CONFIG:-pocket_config.txt}" \
    --exhaustiveness "${PARAM_EXHAUSTIVENESS:-8}" \
    --num-modes "${PARAM_NUM_MODES:-9}"

echo "Node 04a completed"
