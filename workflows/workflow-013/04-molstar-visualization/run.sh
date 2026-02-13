#!/bin/bash
set -e

echo "Starting Node 04: Mol* 3D Visualization"

# Select top designs for Mol* viewer
python select_top_designs.py \
    --input-dir "." \
    --top-n 5 \
    --output-dir "."

echo "Selected top designs for visualization:"
ls -la top_design_*.cif 2>/dev/null || echo "Warning: No top design files generated"

echo "Node 04 completed"
