#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-030/output_files/03-measure-properties"
mkdir -p outputs
for f in properties.json stress_strain.xvg density_prod.xvg; do
    echo "[03-mock] downloading $f"
    python3 -c "import urllib.request; urllib.request.urlretrieve('$BASE/$f', 'outputs/$f')"
done
echo "[03-mock] done"
