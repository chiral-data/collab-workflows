#!/bin/bash
set -e

export PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True
export MPLCONFIGDIR=/tmp/matplotlib

echo "Starting Node 04: Chai-1 Structure Prediction"

cp inputs/* . 2>/dev/null || true

FASTA_FILE=$(ls chai_input.fasta *.fasta *.fa 2>/dev/null | head -1)

if [ -z "$FASTA_FILE" ]; then
    echo "Error: No FASTA input file found — did node 02 preprocessing run?"
    exit 1
fi

echo "Input file: $FASTA_FILE"

python3 predict_chai.py \
    --input "$FASTA_FILE" \
    --num-trunk-recycles "${PARAM_NUM_TRUNK_RECYCLES:-3}" \
    --num-diffusion-timesteps "${PARAM_NUM_DIFFUSION_TIMESTEPS:-50}" \
    --output-dir "./chai_output"

echo "Node 04 completed"
