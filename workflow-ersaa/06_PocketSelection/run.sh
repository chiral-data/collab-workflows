#!/bin/bash
set -e

POCKETS="{{inputs.pockets_json}}"
OUTPUT_JSON="/workspace/out/grids.json"

mkdir -p /workspace/out

python3 /workspace/generate_grids.py "$POCKETS" "$OUTPUT_JSON"
