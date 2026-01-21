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
    echo "✓ Flask found - Preparing Web Application..."
    
    # Explicitly ensure outputs directory exists
    mkdir -p outputs

    # Copy artifacts locally with ERROR CHECKING
    echo "Copying model files for web app..."
    
    # Define paths
    MODEL_SRC="../03-model-training/outputs/model.h5"
    SCALER_SRC="../02-feature-engineering/outputs/scaler.pkl"
    AD_SRC="../02-feature-engineering/outputs/ad_stats.json"

    # Copy with verbose output
    if [ -f "$MODEL_SRC" ]; then
        cp "$MODEL_SRC" outputs/
        echo "✓ Copied model.h5"
    else
        echo "⚠️  WARNING: Could not find $MODEL_SRC"
        # Debug: Look for it elsewhere
        find .. -name "model.h5"
    fi

    if [ -f "$SCALER_SRC" ]; then
        cp "$SCALER_SRC" outputs/
        echo "✓ Copied scaler.pkl"
    else
        echo "⚠️  WARNING: Could not find $SCALER_SRC"
    fi

    if [ -f "$AD_SRC" ]; then
        cp "$AD_SRC" outputs/
        echo "✓ Copied ad_stats.json"
    else
        echo "⚠️  WARNING: Could not find $AD_SRC"
    fi

    echo "Verifying local artifacts in outputs/:"
    ls -l outputs/

    # DISABLED: Flask server blocks workflow execution
    # echo "Starting web application..."
    # python app.py
    echo "Web server disabled for batch workflow execution."
else
    echo "⚠️  Flask not installed - Skipping web server"
    
    # Copy results to 'output' folder if needed for Silva artifacts
    # (Usually Silva picks up from CWD/outputs)
    echo "Node 04 execution complete."
fi
