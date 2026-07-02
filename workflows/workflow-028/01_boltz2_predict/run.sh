#!/bin/bash
set -e

# Workaround: container may run as host UID absent from /etc/passwd,
# causing getpass.getuser() to fail inside PyTorch/Lightning cache init.
export LOGNAME="${LOGNAME:-user}"
export NUMBA_CACHE_DIR="$(pwd)/.numba_cache"
mkdir -p "$NUMBA_CACHE_DIR"

echo "Starting Node 01: Boltz-2 Structure Prediction"

INPUT_FILE=$(ls inputs/*.yaml inputs/*.yml 2>/dev/null | head -1)
if [ -z "$INPUT_FILE" ]; then
    echo "ERROR: No input YAML file found in inputs/"
    exit 1
fi

echo "Input file: $INPUT_FILE"
echo "Accelerator: ${PARAM_ACCELERATOR:-gpu}"

python3 predict.py \
    --input "$INPUT_FILE" \
    --diffusion-samples "${PARAM_DIFFUSION_SAMPLES:-2}" \
    --recycling-steps "${PARAM_RECYCLING_STEPS:-3}" \
    --use-msa-server "${PARAM_USE_MSA_SERVER:-true}" \
    --accelerator "${PARAM_ACCELERATOR:-gpu}"

echo "Node 01 completed"
