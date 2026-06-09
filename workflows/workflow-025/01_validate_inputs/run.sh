#!/bin/bash
set -e

echo "Starting Node 01: Validate Inputs"

python3 validate_inputs.py \
    --pdb-id "${PARAM_TARGET_PDB_ID:-1OKL}" \
    --ligand-input "${PARAM_LIGAND_INPUT:-inputs/ligands.smiles}" \
    --resolution-cutoff "${PARAM_RESOLUTION_CUTOFF:-3.0}"

echo "Node 01 completed"
