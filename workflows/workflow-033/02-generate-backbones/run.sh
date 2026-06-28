#!/bin/bash
set -euo pipefail
echo "=== 02 Generate Backbones ==="
mkdir -p outputs/backbones
python3 generate_backbones.py
echo "Done — outputs: backbones/bb*.pdb  backbone_list.json"
