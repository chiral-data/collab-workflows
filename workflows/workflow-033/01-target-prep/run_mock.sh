#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-033/output_files/01-target-prep"
mkdir -p outputs
for f in target.pdb chain_seq.txt hotspots.json target_a3m.txt prep_report.json; do
    echo "[01-mock] downloading $f"
    curl -fsSL "$BASE/$f" -o "outputs/$f"
done
echo "[01-mock] done"
