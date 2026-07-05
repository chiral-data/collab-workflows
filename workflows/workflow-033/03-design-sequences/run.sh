#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/feat/workflow-033/workflows/workflow-033/output_files/03-design-sequences"
mkdir -p outputs/sequences/bb_0000 outputs/sequences/bb_0001
for f in sequence_manifest.json seq_report.json; do
    echo "[03-mock] downloading $f"
    python3 -c "import urllib.request; urllib.request.urlretrieve('$BASE/$f', 'outputs/$f')"
done
for bb in bb_0000 bb_0001; do
    echo "[03-mock] downloading sequences/$bb/seqs.fa"
    python3 -c "import urllib.request; urllib.request.urlretrieve('$BASE/sequences/$bb/seqs.fa', 'outputs/sequences/$bb/seqs.fa')"
done
echo "[03-mock] done"
