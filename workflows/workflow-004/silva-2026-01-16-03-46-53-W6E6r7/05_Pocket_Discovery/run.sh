#!/bin/bash
set -e
echo "Starting Node 05: Pocket Discovery"

# P2Rank setup (if not in image) or execution
# Assuming P2Rank is installed in the image at /opt/p2rank/p2rank_2.5.1/prank
export P2RANK_HOME=/opt/p2rank/p2rank_2.5.1

python predict_binding_pockets.py
python generate_pocket_analysis.py

echo "Node 05 completed"