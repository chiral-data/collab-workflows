#!/bin/bash
set -e

# Inputs from Silva
PDB_ID="${PDB_ID}"
OUTPUT_DIR="/workspace/out"

mkdir -p "$OUTPUT_DIR"

python3 /workspace/download_receptor.py "$PDB_ID" "$OUTPUT_DIR"
