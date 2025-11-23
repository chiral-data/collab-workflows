#!/bin/bash
set -e

# Silva gives a JSON-style list, convert to comma-separated string
LIGAND_LIST="{{inputs.ligand_ids}}"
RECORD_TYPE="{{inputs.record_type}}"

# Remove brackets and quotes: ["1","2"] → 1,2
CIDS=$(echo "$LIGAND_LIST" | tr -d '[]" ' | tr ',' ',')

OUTPUT_DIR="/workspace/out/ligands"
mkdir -p "$OUTPUT_DIR"

python3 /workspace/download_ligands.py "$CIDS" "$OUTPUT_DIR" "$RECORD_TYPE"
