#!/bin/bash
set -e

echo "Starting Node 04: Uni-Mol Docking V2"

cp inputs/receptor.pdb  . 2>/dev/null || true
cp inputs/ligand.sdf    . 2>/dev/null || true
cp inputs/grid.json     . 2>/dev/null || true
cp inputs/pocket_qc.json . 2>/dev/null || true

python3 dock.py \
    --num-poses    "${PARAM_NUM_POSES:-10}" \
    --weights-path "${PARAM_WEIGHTS_PATH:-/opt/unimol_weights}"

echo "Node 04 completed"
