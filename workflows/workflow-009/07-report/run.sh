#!/bin/bash
set -e

# Find evaluated chemical space from previous step
INPUT_FILE=$(ls inputs/chemspace_evaluated.sdf 2>/dev/null | head -1)

if [ -z "$INPUT_FILE" ]; then
    echo "Error: chemspace_evaluated.sdf not found from previous step"
    exit 1
fi

workflow-run python generate_report.py \
    --input "${INPUT_FILE}" \
    --top-n "${PARAM_TOP_N:-10}" \
    --output-dir "outputs"

# Generate visualization HTML
workflow-run python visualize.py outputs/top_compounds.sdf outputs/top_compounds_report.csv outputs/summary.txt outputs/report_viz.html
