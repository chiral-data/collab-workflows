#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-032/output_files/02-cooling-series"
mkdir -p outputs
for f in cooling_series.json equilibrated.gro density_melt.xvg density_q200.xvg density_q150.xvg density_q80.xvg density_q25.xvg; do
    echo "[02-mock] downloading $f"
    curl -fsSL "$BASE/$f" -o "outputs/$f"
done
echo "[02-mock] done"
