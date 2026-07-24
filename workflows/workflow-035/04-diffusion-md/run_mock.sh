#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-035/output_files/04-diffusion-md"
mkdir -p outputs
for f in msd.xvg diffusion_report.json; do
    echo "[04-mock] downloading $f"
    curl -fsSL "$BASE/$f" -o "outputs/$f"
done
echo "[04-mock] done"
