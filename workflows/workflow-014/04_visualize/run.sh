#!/bin/bash
set -e
echo "Starting Node 04: Generate Report"
mkdir -p outputs

python generate_report.py

echo "Node 04 completed"
