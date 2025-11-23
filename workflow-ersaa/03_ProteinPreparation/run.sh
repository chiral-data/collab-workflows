#!/bin/bash
set -e

INPUT_PDB="/workspace/input/protein.pdb"
OUTPUT_PDBQT="/workspace/out/protein.pdbqt"

mkdir -p /workspace/out

python3 /workspace/prepare_receptor.py "$INPUT_PDB" "$OUTPUT_PDBQT"
