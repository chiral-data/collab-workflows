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
# Copy files from source directories to input/ directory
# ============================================================================
# Copy constructed_library directory from ligand_unpack/outputs/constructed_library/ to input/
LIGAND_UNPACK_DIR="${SCRIPT_DIR}/../ligand_unpack/outputs/constructed_library"
if [ -d "${LIGAND_UNPACK_DIR}" ]; then
    # Remove existing directory if it exists (to avoid "Directory not empty" error)
    if [ -d "${INPUT_DIR}/constructed_library" ]; then
        rm -rf "${INPUT_DIR}/constructed_library"
    fi
    cp -r "${LIGAND_UNPACK_DIR}" "${INPUT_DIR}/"
    echo "✓ Copied constructed_library directory from ligand_unpack/outputs/"
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

# Move real_ligand.sdf from root directory (Silva mounts) to input/
if [ -f "${SCRIPT_DIR}/real_ligand.sdf" ]; then
    mv "${SCRIPT_DIR}/real_ligand.sdf" "${INPUT_DIR}/"
    echo "✓ Moved real_ligand.sdf to input/"
fi

# Copy constructed_library directory from root directory (Silva mounts) to input/
if [ -d "${SCRIPT_DIR}/constructed_library" ]; then
    # Remove existing directory if it exists (to avoid "Directory not empty" error)
    if [ -d "${INPUT_DIR}/constructed_library" ]; then
        rm -rf "${INPUT_DIR}/constructed_library"
    fi
    cp -r "${SCRIPT_DIR}/constructed_library" "${INPUT_DIR}/"
    echo "✓ Copied constructed_library directory to input/"
fi

# ============================================================================
# Run nodes
# ============================================================================
# Change to script directory to ensure relative paths work correctly
cd "${SCRIPT_DIR}"
python3 "${SCRIPT_DIR}/ligand_selection.py"

