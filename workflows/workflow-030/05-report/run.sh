#!/bin/bash
set -euo pipefail

echo "=== Generate Property Report ==="
mkdir -p outputs
python3 generate_report.py
echo ""
echo "Done — outputs: report.html  summary.json"
