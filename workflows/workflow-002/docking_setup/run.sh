#!/bin/bash
set -e

# ============================================================================
# Load global parameters from global_params.json
# ============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GLOBAL_PARAMS_FILE="${SCRIPT_DIR}/../global_params.json"

if [ -f "$GLOBAL_PARAMS_FILE" ]; then
    # Use Python to parse JSON and set environment variables
    export PARAM_PDB_ID=$(python3 -c "import json; print(json.load(open('$GLOBAL_PARAMS_FILE'))['pdb_id'])")
    export PARAM_LIGAND_NAME=$(python3 -c "import json; print(json.load(open('$GLOBAL_PARAMS_FILE'))['ligand_name'])")
    # Also set for compatibility
    export PDB_ID="${PARAM_PDB_ID}"
    export LIGAND_NAME="${PARAM_LIGAND_NAME}"
    echo "✅ Loaded parameters from global_params.json:"
    echo "   PARAM_PDB_ID=${PARAM_PDB_ID}"
    echo "   PARAM_LIGAND_NAME=${PARAM_LIGAND_NAME}"
else
    echo "⚠️  Warning: global_params.json not found at $GLOBAL_PARAMS_FILE"
    echo "   Please set PARAM_PDB_ID and PARAM_LIGAND_NAME environment variables manually."
fi

# ============================================================================
# Create input and outputs directories
# ============================================================================
INPUT_DIR="${SCRIPT_DIR}/input"
OUTPUT_DIR="${SCRIPT_DIR}/outputs"
mkdir -p "${INPUT_DIR}"
mkdir -p "${OUTPUT_DIR}"

# ============================================================================
# Copy PDB files from protein_input/outputs/ to input/ directory
# ============================================================================
PROTEIN_INPUT_DIR="${SCRIPT_DIR}/../protein_input/outputs"
if [ -d "${PROTEIN_INPUT_DIR}" ]; then
    for pdb_file in "${PROTEIN_INPUT_DIR}"/*.pdb; do
        if [ -f "${pdb_file}" ]; then
            cp "${pdb_file}" "${INPUT_DIR}/"
            echo "✓ Copied $(basename "${pdb_file}") from protein_input/outputs/"
        fi
    done
fi

# ============================================================================
# Move files from Silva mounts to input/ directory (fallback)
# ============================================================================
# Move PDB files from root directory (Silva mounts) to input/
if [ -d "${SCRIPT_DIR}" ]; then
    for pdb_file in "${SCRIPT_DIR}"/*.pdb; do
        if [ -f "${pdb_file}" ]; then
            mv "${pdb_file}" "${INPUT_DIR}/"
            echo "✓ Moved $(basename "${pdb_file}") to input/"
        fi
    done
fi

# ============================================================================
# Run nodes
# ============================================================================
# Change to script directory to ensure relative paths work correctly
cd "${SCRIPT_DIR}"
python3 "${SCRIPT_DIR}/ligand_loading.py"
python3 "${SCRIPT_DIR}/extract_chains.py"
python3 "${SCRIPT_DIR}/ligand_center_identification.py"

