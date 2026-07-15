#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-031/output_files/03-insert-penetrant"
mkdir -p outputs
for f in system_with_penetrant.gro topol_penetrant.top penetrant_report.json; do
    echo "[03-mock] downloading $f"
    python3 -c "import urllib.request; urllib.request.urlretrieve('$BASE/$f', 'outputs/$f')"
done
echo "[03-mock] done"
