#!/bin/bash
set -e
echo "Starting Node 01: Receptor Acquisition"

# Dependencies are assumed to be in the docker image
python download_receptor_from_pdb.py
python generate_receptor_report.py

echo "Node 01 completed"