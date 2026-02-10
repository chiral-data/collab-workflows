#!/bin/bash
set -e
echo "Starting Node 02: Model Train"
mkdir -p outputs

python 3_model_train.py
python 4_model_train_HTML.py

echo "Node 02 completed"
