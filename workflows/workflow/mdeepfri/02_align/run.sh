#!/bin/bash
set -e
echo "Starting Node 02: Download Model Weights"
python align.py --model-version "${PARAM_MODEL_VERSION:-1.1}"
echo "Node 02 completed"
