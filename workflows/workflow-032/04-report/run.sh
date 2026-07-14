#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-032/output_files/04-report"
mkdir -p outputs
for f in report.html summary.json; do
    echo "[04-mock] downloading $f"
    python3 -c "import urllib.request; urllib.request.urlretrieve('$BASE/$f', 'outputs/$f')"
done
echo "[04-mock] done"
