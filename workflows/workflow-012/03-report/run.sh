#!/bin/bash
set -e

echo "Starting Node 03: Report Generation"

# Copy input files from silva's inputs/ directory to working directory
cp inputs/* . 2>/dev/null || true

# Install dependencies if needed
pip install -q -r requirements.txt 2>/dev/null || true

# Validate that we have result files from the prediction step
if ! ls confidence_*.json 1>/dev/null 2>&1; then
    echo "Error: No confidence JSON files found from prediction step"
    exit 1
fi

echo "Found result files:"
ls -la confidence_*.json *.pdb *.npz 2>/dev/null || true

# Generate the dashboard
python boltz_dashboard.py .

echo "Generated dashboard files:"
ls -la boltz_dashboard_*.html 2>/dev/null || echo "Warning: No dashboard files generated"

echo "Node 03 completed"
