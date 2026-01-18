#!/bin/bash
set -e
echo "Starting Node 04: Prediction"
mkdir -p outputs

# Run batch predictions
python run_prediction.py

# Generate static report (for backup/export)
python generate_prediction_report.py

# Try to start Flask web application (optional - skip if Flask not installed)
echo ""
echo "========================================="
echo "Checking for Flask web server..."
echo "========================================="

if python -c "import flask" 2>/dev/null; then
    echo "✓ Flask found - Starting web application..."
    python app.py
else
    echo "⚠️  Flask not installed - Skipping web server"
    echo ""
    echo "To enable the interactive web application:"
    echo "  1. Add 'flask' and 'flask-cors' to Docker image"
    echo "  2. Rebuild image"
    echo "  3. Run workflow again"
    echo ""
    echo "For now, view the static dashboard:"
    echo "  outputs/report.html"
    echo ""
    echo "Node 04 completed (batch predictions only)"
fi
