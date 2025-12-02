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
# Copy selected_compounds directory from ligand_selection/outputs/selected_compounds/ to input/
LIGAND_SELECTION_DIR="${SCRIPT_DIR}/../ligand_selection/outputs/selected_compounds"
if [ -d "${LIGAND_SELECTION_DIR}" ]; then
    # Remove existing directory if it exists (to avoid "Directory not empty" error)
    if [ -d "${INPUT_DIR}/selected_compounds" ]; then
        rm -rf "${INPUT_DIR}/selected_compounds"
    fi
    cp -r "${LIGAND_SELECTION_DIR}" "${INPUT_DIR}/"
    echo "✓ Copied selected_compounds directory from ligand_selection/outputs/"
fi

# Copy true_ligand.sdf from true_ligand_addition/outputs/ to input/
TRUE_LIGAND_FILE="${SCRIPT_DIR}/../true_ligand_addition/outputs/true_ligand.sdf"
if [ -f "${TRUE_LIGAND_FILE}" ]; then
    cp "${TRUE_LIGAND_FILE}" "${INPUT_DIR}/"
    echo "✓ Copied true_ligand.sdf from true_ligand_addition/outputs/"
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

# Remove selected_compounds from root directory if Silva mounted it there
# (Silva may mount outputs/selected_compounds as root/selected_compounds)
# We don't use root/selected_compounds - we only use input/selected_compounds from run.sh copy
if [ -d "${SCRIPT_DIR}/selected_compounds" ]; then
    rm -rf "${SCRIPT_DIR}/selected_compounds"
    echo "✓ Removed selected_compounds from root directory (Silva mount)"
fi

# Remove true_ligand.sdf from root directory if Silva mounted it there
# (Silva may mount outputs/selected_compounds/true_ligand.sdf as root/true_ligand.sdf due to *.sdf pattern in outputs)
# We only use input/true_ligand.sdf - root directory should not have it
if [ -f "${SCRIPT_DIR}/true_ligand.sdf" ]; then
    rm -f "${SCRIPT_DIR}/true_ligand.sdf"
    echo "✓ Removed true_ligand.sdf from root directory (Silva mount)"
fi

# ============================================================================
# Run nodes
# ============================================================================
# Change to script directory to ensure relative paths work correctly
cd "${SCRIPT_DIR}"
python3 "${SCRIPT_DIR}/prepare_ligands.py"
python3 "${SCRIPT_DIR}/ligand_view.py"

# ============================================================================
# Post-run cleanup: Remove true_ligand.sdf from root directory if it was created
# ============================================================================
# After processing, if true_ligand.sdf exists in root directory, remove it
# (input/true_ligand.sdf should remain, only root directory version is removed)
if [ -f "${SCRIPT_DIR}/true_ligand.sdf" ]; then
    rm -f "${SCRIPT_DIR}/true_ligand.sdf"
    echo "✓ Removed true_ligand.sdf from root directory (post-run cleanup)"
fi

