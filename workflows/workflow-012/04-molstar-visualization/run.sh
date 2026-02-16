#!/bin/bash
set -e

echo "Starting Node 04: Mol* 3D Visualization"

# Copy input files from silva's inputs/ directory to working directory
cp inputs/* . 2>/dev/null || true

# Select top models by confidence score for Mol* viewer
python select_top_models.py \
    --input-dir "." \
    --top-n 5 \
    --output-dir "."

echo "Selected top models for visualization:"
ls -la top_model_*.pdb 2>/dev/null || echo "Warning: No top model files generated"

echo "Node 04 completed"
