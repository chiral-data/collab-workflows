#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-036/output_files/04-report"
mkdir -p outputs
for f in report.html summary.json; do
    echo "[04-mock] downloading $f"
    curl -fsSL "$BASE/$f" -o "outputs/$f"
done
echo "[04-mock] done"
