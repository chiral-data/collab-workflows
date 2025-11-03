#!/bin/bash
set -e

echo "=== Starting Protein and Ligand Preparation ==="

# Fix protein
conda run -n prep_env python prepare_protein.py input/protein.pdb output/protein_fixed.pdbqt

# Convert ligands
conda run -n prep_env python prepare_ligands.py input/ligands output/ligands_pdbqt

echo "=== Preparation Complete ==="