#!/bin/bash
set -e

RECEPTOR="/workspace/input/receptor.pdbqt"
POCKET="/workspace/input/pocket.json"
LIGAND="/workspace/input/ligand.pdbqt"
BOX_SIZE="${BOX_SIZE:-80}"

OUTPUT_DIR="/workspace/out"
mkdir -p "$OUTPUT_DIR"

python3 /workspace/run_vina.py \
    "$RECEPTOR" \
    "$LIGAND" \
    "$POCKET" \
    "$BOX_SIZE" \
    "$OUTPUT_DIR"
