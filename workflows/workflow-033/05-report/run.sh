#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/feat/workflow-033/workflows/workflow-033/output_files/05-report"
mkdir -p outputs
for f in candidates.csv summary.json report.html; do
    echo "[05-mock] downloading $f"
    python3 -c "import urllib.request; urllib.request.urlretrieve('$BASE/$f', 'outputs/$f')"
done
echo "[05-mock] done"
