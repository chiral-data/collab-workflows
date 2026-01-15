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

# Generate visualization HTML
workflow-run python visualize.py top_compounds.sdf top_compounds_report.csv summary.txt report_viz.html
