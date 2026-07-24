#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-032/output_files/02-cooling-series"
mkdir -p outputs

fetch() {
    python3 -c "
import sys, time, urllib.error, urllib.request
url, out = sys.argv[1], sys.argv[2]
for attempt in range(1, 6):
    try:
        urllib.request.urlretrieve(url, out)
        break
    except urllib.error.HTTPError as e:
        if attempt == 5 or e.code not in (429, 500, 502, 503, 504):
            raise
        wait = 2 ** attempt
        print(f'  HTTP {e.code}; retry {attempt}/5 in {wait}s', flush=True)
        time.sleep(wait)
" "$1" "$2"
    sleep 1
}

for f in cooling_series.json equilibrated.gro density_melt.xvg density_q200.xvg density_q150.xvg density_q80.xvg density_q25.xvg; do
    echo "[02-mock] downloading $f"
    fetch "$BASE/$f" "outputs/$f"
done
echo "[02-mock] done"
