#!/bin/bash
set -e

P2RANK_DIR="/workspace/input/p2rank_output"
OUTPUT_JSON="/workspace/out/pocket.json"

mkdir -p /workspace/out

python3 /workspace/generate_grids.py "$P2RANK_DIR" "$OUTPUT_JSON"
