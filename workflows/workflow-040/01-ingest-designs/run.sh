#!/bin/bash
set -e
echo "Starting Node 01: Ingest Designs"
python script.py "$@"
echo "Node 01 completed"
