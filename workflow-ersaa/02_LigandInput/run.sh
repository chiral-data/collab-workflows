#!/bin/bash
set -e

CIDS="${CIDS}"
RECORD_TYPE="${RECORD_TYPE:-3d}"

OUTPUT_DIR="/workspace/out"
mkdir -p "$OUTPUT_DIR"

python3 /workspace/download_ligands.py "$CIDS" "$OUTPUT_DIR" "$RECORD_TYPE"
