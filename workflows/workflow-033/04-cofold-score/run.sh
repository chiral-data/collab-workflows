#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-033/output_files/04-cofold-score"
mkdir -p outputs/complexes

# A failed download fails the job; silva re-runs it (`retry` in job.toml).
fetch() {
    python3 -c "import sys, urllib.request; urllib.request.urlretrieve(sys.argv[1], sys.argv[2])" "$1" "$2"
    sleep 1
}

for f in scores.json manifest.json cofold_report.json; do
    echo "[04-mock] downloading $f"
    fetch "$BASE/$f" "outputs/$f"
done
for f in bb_0000_seq000.cif bb_0000_seq001.cif bb_0001_seq000.cif bb_0001_seq001.cif; do
    echo "[04-mock] downloading complexes/$f"
    fetch "$BASE/complexes/$f" "outputs/complexes/$f"
done
echo "[04-mock] done"
