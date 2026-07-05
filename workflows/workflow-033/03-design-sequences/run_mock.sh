#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-033/output_files/03-design-sequences"
mkdir -p outputs/sequences/bb_0000 outputs/sequences/bb_0001
for f in sequence_manifest.json seq_report.json; do
    echo "[03-mock] downloading $f"
    curl -fsSL "$BASE/$f" -o "outputs/$f"
done
for bb in bb_0000 bb_0001; do
    echo "[03-mock] downloading sequences/$bb/seqs.fa"
    curl -fsSL "$BASE/sequences/$bb/seqs.fa" -o "outputs/sequences/$bb/seqs.fa"
done
echo "[03-mock] done"
