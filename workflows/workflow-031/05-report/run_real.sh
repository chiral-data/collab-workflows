#!/bin/bash
set -euo pipefail
echo "=== Generate Barrier Performance Report ==="
mkdir -p outputs
python3 generate_report.py
echo "Done — outputs: report.html  summary.json"
