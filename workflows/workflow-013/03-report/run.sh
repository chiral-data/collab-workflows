#!/bin/bash
set -e

echo "Starting Node 03: BoltzGen Report Generation"

# Install dependencies if needed
pip install -q -r requirements.txt 2>/dev/null || true

echo "Found result files:"
ls -la *.cif *.csv *.pdf 2>/dev/null || true

# Generate the dashboard
python boltzgen_dashboard.py .

echo "Generated dashboard files:"
ls -la boltzgen_dashboard_*.html 2>/dev/null || echo "Warning: No dashboard files generated"

echo "Node 03 completed"
