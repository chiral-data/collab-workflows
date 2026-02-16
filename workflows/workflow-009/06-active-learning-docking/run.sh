#!/bin/bash
set -e

# Find chemical space from previous step
CHEMSPACE_FILE=$(ls inputs/chemspace.pkl 2>/dev/null | head -1)

if [ -z "$CHEMSPACE_FILE" ]; then
    echo "Error: chemspace.pkl not found from previous step"
    exit 1
fi

workflow-run python run_active_learning.py \
    --chemspace "${CHEMSPACE_FILE}" \
    --initial-molecules "${PARAM_INITIAL_MOLECULES:-10}" \
    --num-cycles "${PARAM_NUM_CYCLES:-3}" \
    --molecules-per-cycle "${PARAM_MOLECULES_PER_CYCLE:-50}" \
    --model-type "${PARAM_MODEL_TYPE:-gaussian_process}" \
    --query-type "${PARAM_QUERY_TYPE:-UCB}" \
    --output outputs/chemspace_evaluated.sdf
