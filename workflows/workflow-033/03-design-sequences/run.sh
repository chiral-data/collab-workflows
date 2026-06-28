#!/bin/bash
set -euo pipefail
echo "=== 03 Design Sequences ==="
mkdir -p outputs/sequences
python3 design_sequences.py
echo "Done — outputs: sequences/*.json  sequence_manifest.json"
