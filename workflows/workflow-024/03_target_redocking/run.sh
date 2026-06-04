#!/bin/bash
set -e

echo "Starting Node 03: Target Redocking QC"

cp inputs/receptor.pdb . 2>/dev/null || true
cp inputs/receptor.pdbqt . 2>/dev/null || true
cp inputs/native_ligand.pdbqt . 2>/dev/null || true
cp inputs/pocket_config.txt . 2>/dev/null || true
cp inputs/optimized_screening_library.pdbqt . 2>/dev/null || true

python3 redock_target.py \
    --receptor-pdbqt "${PARAM_RECEPTOR_PDBQT:-receptor.pdbqt}" \
    --native-ligand-pdbqt "${PARAM_NATIVE_LIGAND_PDBQT:-native_ligand.pdbqt}" \
    --receptor-pdb "${PARAM_RECEPTOR_PDB:-receptor.pdb}" \
    --pocket-config "${PARAM_POCKET_CONFIG:-pocket_config.txt}" \
    --cnn-score-threshold "${PARAM_CNN_SCORE_THRESHOLD:-0.90}" \
    --rmsd-threshold "${PARAM_RMSD_THRESHOLD:-2.0}"

echo "Node 03 completed"
