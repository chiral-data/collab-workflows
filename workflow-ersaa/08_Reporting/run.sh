#!/bin/bash
set -e

RESULTS="{{inputs.vina_results}}"
OUTDIR="/workspace/out"

mkdir -p "$OUTDIR"

python3 /workspace/rank_vina.py "$RESULTS"

# Move output to Silva output directory
mv binding_affinities.xlsx "$OUTDIR/results.xlsx"
