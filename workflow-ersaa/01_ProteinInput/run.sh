#!/bin/bash
set -e

PDB_ID="{{inputs.protein_id}}"

OUTPUT_DIR="/workspace/out"
mkdir -p "$OUTPUT_DIR"

python3 /workspace/download_receptor.py "$PDB_ID" "$OUTPUT_DIR"
