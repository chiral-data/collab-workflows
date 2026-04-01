#!/bin/bash
set -e
echo "Starting Node 02: Ligand Collection"

python download_ligands_from_pubchem.py
python generate_ligand_report.py

echo "Node 02 completed"