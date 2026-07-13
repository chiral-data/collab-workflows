#!/bin/bash
set -e
mkdir -p outputs
python3 generate_report.py
echo ""
echo "Done — outputs: report.html  summary.json"
