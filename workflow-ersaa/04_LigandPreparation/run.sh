#!/bin/bash
set -e

INPUT_DIR="/workspace/input/ligands"
OUTPUT_DIR="/workspace/out"

mkdir -p "$OUTPUT_DIR"

python3 /workspace/prepare_ligands.py "$INPUT_DIR" "$OUTPUT_DIR"
