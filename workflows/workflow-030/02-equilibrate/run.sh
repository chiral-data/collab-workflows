#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-030/output_files/02-equilibrate"
mkdir -p outputs
for f in equilibrated.gro density.xvg equil_report.json; do
    echo "[02-mock] downloading $f"
    python3 -c "import urllib.request; urllib.request.urlretrieve('$BASE/$f', 'outputs/$f')"
done
echo "[02-mock] done"
