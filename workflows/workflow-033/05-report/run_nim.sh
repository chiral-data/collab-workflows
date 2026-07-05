#!/bin/bash
set -euo pipefail
echo "=== 05 Report ==="
mkdir -p outputs
python3 generate_report.py
echo "Done — outputs: report.html  summary.json  ranked_binders.csv"
