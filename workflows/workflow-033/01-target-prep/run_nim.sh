#!/bin/bash
set -euo pipefail
export NGC_API_KEY="${NGC_API_KEY:-}"
echo "=== 01 Target Prep ==="
mkdir -p outputs
python3 target_prep.py
echo "Done — outputs: target.pdb  chain_seq.txt  hotspots.json  target_a3m.txt  prep_report.json"
