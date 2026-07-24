#!/bin/bash
set -e
mkdir -p outputs
python3 measure_tg.py
echo ""
echo "Done — outputs: tg_report.json"
