#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-033/output_files/04-cofold-score"
mkdir -p outputs/complexes

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

for f in scores.json manifest.json cofold_report.json; do
    echo "[04-mock] downloading $f"
    fetch "$BASE/$f" "outputs/$f"
done
for f in bb_0000_seq000.cif bb_0000_seq001.cif bb_0001_seq000.cif bb_0001_seq001.cif; do
    echo "[04-mock] downloading complexes/$f"
    fetch "$BASE/complexes/$f" "outputs/complexes/$f"
done
echo "[04-mock] done"
