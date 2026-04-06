#!/bin/bash
set -e
echo "=== Node 1: Complex splitting ==="

ROOT_DIR="$(dirname "$PWD")"
INPUT_PDB_PATH="inputs/${INPUT_PDB}"

if [ ! -f "$INPUT_PDB_PATH" ]; then
	if [ -f "$ROOT_DIR/input_files/${INPUT_PDB}" ]; then
		INPUT_PDB_PATH="$ROOT_DIR/input_files/${INPUT_PDB}"
	elif [ -f "$ROOT_DIR/input_files/5B8C.pdb" ]; then
		INPUT_PDB_PATH="$ROOT_DIR/input_files/5B8C.pdb"
	else
		echo "Error: could not locate input PDB. Expected inputs/${INPUT_PDB}"
		exit 1
	fi
fi

python 1_complex_split_science.py --input_pdb "$INPUT_PDB_PATH" --antibody_chains "$ANTIBODY_CHAINS" --output_dir "outputs"
python 2_complex_split_html.py --data_json "outputs/data.json" --output_html "outputs/report.html"
echo "Node 1 finished successfully"