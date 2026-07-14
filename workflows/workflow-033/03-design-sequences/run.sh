#!/bin/bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-033/output_files/03-design-sequences"
mkdir -p outputs/sequences/bb_0000 outputs/sequences/bb_0001

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

for f in sequence_manifest.json seq_report.json; do
    echo "[03-mock] downloading $f"
    fetch "$BASE/$f" "outputs/$f"
done
for bb in bb_0000 bb_0001; do
    echo "[03-mock] downloading sequences/$bb/seqs.fa"
    fetch "$BASE/sequences/$bb/seqs.fa" "outputs/sequences/$bb/seqs.fa"
done
echo "[03-mock] done"
