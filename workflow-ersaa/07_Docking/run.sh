#!/bin/bash
set -e

RECEPTOR="{{inputs.receptor_pdbqt}}"
LIGANDS="{{inputs.ligands_prepared}}"
GRIDS="{{inputs.grid_boxes}}"
BOX_SIZE="${BOX_SIZE:-80}"

OUTPUT_DIR="/workspace/out"
mkdir -p "$OUTPUT_DIR"

python3 /workspace/run_vina.py \
    "$RECEPTOR" \
    "$LIGAND" \
    "$POCKET" \
    "$BOX_SIZE" \
    "$OUTPUT_DIR"
