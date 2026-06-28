#!/bin/bash
set -euo pipefail
echo "=== 04 Co-fold and Score ==="
mkdir -p outputs/complexes
python3 cofold_score.py
echo "Done — outputs: complexes/*.cif  scores.json  manifest.json"
