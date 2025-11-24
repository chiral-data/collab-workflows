#!/bin/bash
set -e

INPUT_PDB="{{inputs.receptor_pdb}}"
OUTPUT_JSON="/workspace/out/pockets.json"

mkdir -p /workspace/out

python3 /workspace/run_p2rank.py "$INPUT_PDB" "$OUTPUT_JSON"
