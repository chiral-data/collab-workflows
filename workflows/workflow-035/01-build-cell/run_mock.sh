#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-035/output_files/01-build-cell"
mkdir -p outputs
for f in system.gro topol.top cell.pdb build_report.json; do
    echo "[01-mock] downloading $f"
    curl -fsSL "$BASE/$f" -o "outputs/$f"
done
echo "[01-mock] done"
