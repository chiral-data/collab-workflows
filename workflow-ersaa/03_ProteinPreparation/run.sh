#!/bin/bash
set -e

INPUT_PDB="{{inputs.receptor_pdb}}"
OUTPUT_PDBQT="/workspace/out/receptor.pdbqt"

mkdir -p /workspace/out

python3 /workspace/prepare_receptor.py "$INPUT_PDB" "$OUTPUT_PDBQT"
