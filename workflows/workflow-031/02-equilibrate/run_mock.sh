#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-031/output_files/02-equilibrate"
mkdir -p outputs
for f in equilibrated.gro equil.xtc density.xvg equil_report.json; do
    echo "[02-mock] downloading $f"
    curl -fsSL "$BASE/$f" -o "outputs/$f"
done
echo "[02-mock] done"
