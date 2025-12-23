#!/bin/bash
set -e

workflow-run python validate_ligand.py \
    --input "${PARAM_LIGAND_FILE}" \
    --output ligand.sdf
