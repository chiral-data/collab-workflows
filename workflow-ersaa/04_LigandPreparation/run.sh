#!/bin/bash
set -e

INPUT_DIR="{{inputs.ligands_folder}}"
OUTPUT_DIR="/workspace/out/ligands_prepared"

mkdir -p "$OUTPUT_DIR"

python3 /workspace/prepare_ligands.py "$INPUT_DIR" "$OUTPUT_DIR"
