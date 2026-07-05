#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-033/output_files/02-generate-backbones"
mkdir -p outputs/backbones
for f in backbone_list.json gen_report.json; do
    echo "[02-mock] downloading $f"
    python3 -c "import urllib.request; urllib.request.urlretrieve('$BASE/$f', 'outputs/$f')"
done
for f in bb_0000.pdb bb_0001.pdb; do
    echo "[02-mock] downloading backbones/$f"
    python3 -c "import urllib.request; urllib.request.urlretrieve('$BASE/backbones/$f', 'outputs/backbones/$f')"
done
echo "[02-mock] done"
