#!/bin/bash
set -e

# Find evaluated chemical space from previous step
INPUT_FILE=$(ls chemspace_evaluated.sdf 2>/dev/null | head -1)

if [ -z "$INPUT_FILE" ]; then
    echo "Error: chemspace_evaluated.sdf not found from previous step"
    exit 1
fi

workflow-run python generate_report.py \
    --input "${INPUT_FILE}" \
    --top-n "${PARAM_TOP_N:-10}" \
    --output-dir "."
