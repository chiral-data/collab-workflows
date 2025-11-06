#!/bin/bash
set -e

#python download.py
python3 protein_preparation.py
python3 extract_2TK_from_pdb.py
python3 generate_ligand_mods.py
