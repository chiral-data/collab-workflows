#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-031/output_files/04-diffusion-md"
mkdir -p outputs
for f in msd.xvg diffusion_report.json; do
    echo "[04-mock] downloading $f"
    python3 -c "import urllib.request; urllib.request.urlretrieve('$BASE/$f', 'outputs/$f')"
done
echo "[04-mock] done"
