#!/bin/bash
set -e
echo "Starting Node 01: Validate Inputs"
mkdir -p outputs

python validate.py

echo "Node 01 completed"
