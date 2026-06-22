#!/bin/bash
set -e

# Workaround: container may run as host UID absent from /etc/passwd,
# causing getpass.getuser() to fail inside PyTorch/Lightning cache init.
export LOGNAME="${LOGNAME:-user}"
export NUMBA_CACHE_DIR="/tmp/numba_cache"

echo "Starting Node 01: Boltz-2 Structure Prediction"

cp inputs/* . 2>/dev/null || true

INPUT_FILE=$(ls *.yaml *.yml 2>/dev/null | head -1)
if [ -z "$INPUT_FILE" ]; then
    echo "ERROR: No input YAML file found in inputs/"
    exit 1
fi

echo "Input file: $INPUT_FILE"

python3 predict.py \
    --input "$INPUT_FILE" \
    --diffusion-samples "${PARAM_DIFFUSION_SAMPLES:-10}" \
    --recycling-steps "${PARAM_RECYCLING_STEPS:-5}" \
    --use-msa-server "${PARAM_USE_MSA_SERVER:-true}"

echo "Node 01 completed"
