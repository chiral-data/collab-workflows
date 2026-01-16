#!/bin/bash
#
# Docking Report Job Script Template
# This script reads configuration from a JSON file
#

# assume job_config.json always exist
CONFIG_FILE="job_config.json"

# Check if config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "Error: Config file '$CONFIG_FILE' not found"
    exit 1
fi

# Parse JSON config
JOB_ID=$(jq -r '.job_id' "$CONFIG_FILE")
OUTPUT_DIR=$(jq -r '.output.directory // "outputs"' "$CONFIG_FILE")

echo "Starting Docking Report job: $JOB_ID"
echo "Output directory: $OUTPUT_DIR"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Verify all input files exist
echo "Verifying input files..."
jq -c '.inputs[]' "$CONFIG_FILE" | while read -r input; do
    FILENAME=$(echo "$input" | jq -r '.filename')
    
    if [ ! -f "$FILENAME" ]; then
        echo "Error: Input file '$FILENAME' not found"
        exit 1
    fi
    echo "Found: $FILENAME"
done

# Parse analysis type parameter
ANALYSIS_TYPE=$(jq -r '.parameters.analysis_type // "autodock_vina"' "$CONFIG_FILE")

echo "=== Docking Report Parameters ==="
echo "Analysis type: $ANALYSIS_TYPE"

# Run the appropriate analysis script
if [ "$ANALYSIS_TYPE" = "diffdock_pp" ]; then
    echo "=== Running DiffDock-PP Analysis ==="
    python3 /workspace/diffdock_pp_dashboard.py
elif [ "$ANALYSIS_TYPE" = "diffdock" ]; then
    echo "=== Running DiffDock Analysis ==="
    python3 /workspace/diffdock_dashboard.py
elif [ "$ANALYSIS_TYPE" = "autodock_vina" ]; then
    echo "=== Running AutoDock Vina Analysis ==="
    python3 /workspace/autodock_vina_dashboard.py
else
    echo "Error: Unknown analysis type: $ANALYSIS_TYPE"
    echo "Supported types: diffdock_pp, diffdock, autodock_vina"
    exit 1
fi

# Check if the analysis was successful
if [ $? -eq 0 ]; then
    echo "=== Analysis completed successfully ==="
    
    # Copy all generated files to output directory
    # The Python script generates files in the current directory
    echo "Copying results to output directory..."
    cp *.html "$OUTPUT_DIR/" 2>/dev/null || true
    cp *.png "$OUTPUT_DIR/" 2>/dev/null || true
    cp *.json "$OUTPUT_DIR/" 2>/dev/null || true
    
    # List generated files
    echo "Generated files:"
    ls -la "$OUTPUT_DIR/"
else
    echo "Error: Analysis failed"
    exit 1
fi

echo "Docking report generation completed. Results are available in: $OUTPUT_DIR"