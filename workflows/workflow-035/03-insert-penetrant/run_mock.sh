#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-035/output_files/03-insert-penetrant"
mkdir -p outputs
for f in system_with_penetrant.gro topol_penetrant.top penetrant_report.json; do
    echo "[03-mock] downloading $f"
    curl -fsSL "$BASE/$f" -o "outputs/$f"
done
echo "[03-mock] done"
