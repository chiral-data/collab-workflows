#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-030/output_files/04-apply-gf-correction"
mkdir -p outputs
for f in corrected_properties.json; do
    echo "[04-mock] downloading $f"
    curl -fsSL "$BASE/$f" -o "outputs/$f"
done
echo "[04-mock] done"
