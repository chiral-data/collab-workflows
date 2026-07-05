#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/feat/workflow-033/workflows/workflow-033/output_files/04-cofold-score"
mkdir -p outputs/complexes
for f in scores.json manifest.json cofold_report.json; do
    echo "[04-mock] downloading $f"
    python3 -c "import urllib.request; urllib.request.urlretrieve('$BASE/$f', 'outputs/$f')"
done
for f in bb_0000_seq000.cif bb_0000_seq001.cif bb_0001_seq000.cif bb_0001_seq001.cif; do
    echo "[04-mock] downloading complexes/$f"
    python3 -c "import urllib.request; urllib.request.urlretrieve('$BASE/complexes/$f', 'outputs/complexes/$f')"
done
echo "[04-mock] done"
