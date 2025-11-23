#!/bin/bash
set -e

PROTEIN_PDB="/workspace/input/protein_prepared.pdb"
OUTPUT_DIR="/workspace/out"

mkdir -p "$OUTPUT_DIR"

python3 /workspace/run_p2rank.py "$PROTEIN_PDB" "$OUTPUT_DIR"
