#!/bin/bash
#
# Boltz-2 Job Script Template
# This script reads configuration from a JSON file
#

# # Check if config file is provided
# if [ "$#" -ne 1 ]; then
#     echo "Usage: $0 <config.json>"
#     exit 1
# fi
#
# CONFIG_FILE="$1"

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

echo "Starting Boltz job: $JOB_ID"
echo "Output directory: $OUTPUT_DIR"

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Verify all input files exist (files should be uploaded already)
echo "Verifying input files..."
jq -c '.inputs[]' "$CONFIG_FILE" | while read -r input; do
    FILENAME=$(echo "$input" | jq -r '.filename')

    if [ ! -f "$FILENAME" ]; then
        echo "Error: Input file '$FILENAME' not found"
        exit 1
    fi
    echo "Found: $FILENAME"
done

# Build Boltz command with parameters
build_boltz_command() {
    local input_file="$1"
    local cmd="python3 -m boltz.main predict $input_file"

    # Add parameters from config with new defaults
    cmd="$cmd --use_msa_server"

    local output_format=$(jq -r '.parameters.output_format // "pdb"' "$CONFIG_FILE")
    cmd="$cmd --output_format $output_format"

    local recycling_steps=$(jq -r '.parameters.recycling_steps // 3' "$CONFIG_FILE")
    cmd="$cmd --recycling_steps $recycling_steps"

    local sampling_steps=$(jq -r '.parameters.sampling_steps // 200' "$CONFIG_FILE")
    cmd="$cmd --sampling_steps $sampling_steps"

    local diffusion_samples=$(jq -r '.parameters.diffusion_samples // 1' "$CONFIG_FILE")
    cmd="$cmd --diffusion_samples $diffusion_samples"

    local step_scale=$(jq -r '.parameters.step_scale // 1.638' "$CONFIG_FILE")
    cmd="$cmd --step_scale $step_scale"

    local devices=$(jq -r '.parameters.devices // 1' "$CONFIG_FILE")
    cmd="$cmd --devices $devices"

    local accelerator=$(jq -r '.parameters.accelerator // "gpu"' "$CONFIG_FILE")
    cmd="$cmd --accelerator $accelerator"

    echo "$cmd"
}

# Run Boltz predictions
echo "Running Boltz predictions..."
jq -c '.inputs[]' "$CONFIG_FILE" | while read -r input; do
    FILENAME=$(echo "$input" | jq -r '.filename')
    # Extract base name without extension for result directory name
    BASENAME=$(basename "$FILENAME" | sed 's/\.[^.]*$//')

    echo "Processing $FILENAME..."
    cmd=$(build_boltz_command "$FILENAME")
    echo "Command: $cmd"
    eval "$cmd"

    # Copy results to output directory
    RESULTS_DIR="boltz_results_${BASENAME}/predictions/${BASENAME}"
    if [ -d "$RESULTS_DIR" ]; then
        echo "Copying results for $BASENAME..."
        cp "$RESULTS_DIR"/* "$OUTPUT_DIR/" 2>/dev/null || true
    else
        echo "Warning: Results directory not found: $RESULTS_DIR"
    fi
done

# Clean up intermediate files
echo "Cleaning up intermediate files..."
jq -c '.inputs[]' "$CONFIG_FILE" | while read -r input; do
    FILENAME=$(echo "$input" | jq -r '.filename')
    BASENAME=$(basename "$FILENAME" | sed 's/\.[^.]*$//')
    rm -rf "boltz_results_${BASENAME}"
done

echo "Prediction completed. Results are available in: $OUTPUT_DIR"
ls -la "$OUTPUT_DIR/"
