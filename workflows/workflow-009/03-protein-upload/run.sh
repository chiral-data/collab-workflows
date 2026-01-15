#!/bin/bash
set -e

workflow-run python validate_protein.py \
    --input "${PARAM_PROTEIN_FILE}" \
    --output protein.pdb
