#!/bin/bash
set -e

# Workaround: container runs as host UID which may not exist in /etc/passwd,
# causing getpass.getuser() to fail inside PyTorch/Lightning cache init.
export LOGNAME="${LOGNAME:-user}"
export NUMBA_CACHE_DIR="/tmp/numba_cache"

echo "Starting Node 02: Boltz-2 Structure Prediction"

# Copy input files from silva's inputs/ directory to working directory
cp inputs/* . 2>/dev/null || true

# Find validated YAML input from previous node
INPUT_FILE=$(ls validated_input.yaml *.yaml 2>/dev/null | head -1)

if [ -z "$INPUT_FILE" ]; then
    echo "Error: No input YAML file found from previous step"
    exit 1
fi

echo "Input file: $INPUT_FILE"

python run_prediction.py \
    --input "$INPUT_FILE" \
    --diffusion-samples "${PARAM_DIFFUSION_SAMPLES:-10}" \
    --recycling-steps "${PARAM_RECYCLING_STEPS:-5}" \
    --use-msa-server "${PARAM_USE_MSA_SERVER:-true}" \
    --output-dir "."

echo "Node 02 completed"
