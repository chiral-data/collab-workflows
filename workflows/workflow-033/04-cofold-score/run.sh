#!/bin/bash
set -euo pipefail
export NGC_API_KEY="${NGC_API_KEY:-}"
echo "=== 04 Co-fold and Score ==="
mkdir -p outputs/complexes
python3 cofold_score.py
echo "Done — outputs: complexes/*.cif  scores.json  manifest.json"
