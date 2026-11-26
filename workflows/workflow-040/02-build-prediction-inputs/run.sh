#!/bin/bash
set -e
echo "Starting Node 02: Build Prediction Inputs"
python script.py "$@"
echo "Node 02 completed"
