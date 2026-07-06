#!/bin/bash
set -e

echo "Starting Node 03: Pocket QC + Grid Preparation"

cp inputs/* . 2>/dev/null || true

python3 pocket_qc_grid.py \
    --plddt-threshold "${PARAM_PLDDT_THRESHOLD:-70.0}" \
    --box-size        "${PARAM_BOX_SIZE:-22.5}" \
    --pocket-rank     "${PARAM_POCKET_RANK:-1}"

echo "Node 03 completed"
