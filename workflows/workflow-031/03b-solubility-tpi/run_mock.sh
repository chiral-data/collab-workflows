#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-031/output_files/03b-solubility-tpi"
mkdir -p outputs
for f in solubility_report.json tpi.xvg; do
    echo "[03b-mock] downloading $f"
    curl -fsSL "$BASE/$f" -o "outputs/$f"
done
echo "[03b-mock] done"
