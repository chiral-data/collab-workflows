#!/bin/bash
set -e

# Workaround: container runs as host UID which may not exist in /etc/passwd,
# causing getpass.getuser() to fail inside PyTorch/Lightning cache init.
export LOGNAME="${LOGNAME:-user}"
export NUMBA_CACHE_DIR="/tmp/numba_cache"

echo "Starting Node 03: Boltz-2 Structure Prediction"

# Copy preprocessed inputs from node 02
cp ../02-preprocessing/boltz_input.yaml . 2>/dev/null || true

# Expect boltz_input.yaml specifically from node 02
INPUT_FILE=$(ls boltz_input.yaml *.yaml 2>/dev/null | head -1)

if [ -z "$INPUT_FILE" ]; then
    echo "Error: No YAML input file found — did node 02 preprocessing run?"
    exit 1
fi

echo "Input file: $INPUT_FILE"

python3 predict_boltz.py \
    --input "$INPUT_FILE" \
    --diffusion-samples "${PARAM_DIFFUSION_SAMPLES:-2}" \
    --recycling-steps "${PARAM_RECYCLING_STEPS:-3}" \
    --use-msa-server "${PARAM_USE_MSA_SERVER:-true}" \
    --output-dir "."

echo "Node 03 completed"
