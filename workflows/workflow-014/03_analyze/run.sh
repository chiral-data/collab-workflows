#!/bin/bash
set -e
echo "Starting Node 03: Filter and Rank"
mkdir -p outputs

python analyze.py

echo "Node 03 completed"
