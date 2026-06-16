#!/bin/bash
set -e

echo "Starting Node 02: Prepare Structures"

# Flatten inputs from upstream dependency outputs
cp inputs/receptor.pdb . 2>/dev/null || true
cp inputs/native_ligand.sdf . 2>/dev/null || true
cp inputs/validated_ligands.smiles . 2>/dev/null || true
cp inputs/validation_report.json . 2>/dev/null || true

python3 prepare_structures.py \
    --receptor "${PARAM_RECEPTOR:-receptor.pdb}" \
    --native-ligand "${PARAM_NATIVE_LIGAND:-native_ligand.sdf}" \
    --ligands "${PARAM_LIGAND_INPUT:-validated_ligands.smiles}" \
    --box-padding "${PARAM_BOX_PADDING:-5.0}"

echo "Node 02 completed"
