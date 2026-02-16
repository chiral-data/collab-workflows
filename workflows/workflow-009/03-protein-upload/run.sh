#!/bin/bash
set -e

workflow-run python validate_protein.py \
    --input "inputs/${PARAM_PROTEIN_FILE}" \
    --output outputs/protein.pdb
