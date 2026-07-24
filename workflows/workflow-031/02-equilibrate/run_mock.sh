#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-031/output_files/02-equilibrate"
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

for f in equilibrated.gro equil.xtc density.xvg equil_report.json; do
    echo "[02-mock] downloading $f"
    fetch "$BASE/$f" "outputs/$f"
done
echo "[02-mock] done"
