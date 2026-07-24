#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-036/output_files/03-measure-tg"
mkdir -p outputs
for f in tg_report.json; do
    echo "[03-mock] downloading $f"
    curl -fsSL "$BASE/$f" -o "outputs/$f"
done
echo "[03-mock] done"
