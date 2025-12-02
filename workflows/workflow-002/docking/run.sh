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

# Clean up any stray files in outputs/ directory
# Remove any .sdf or .txt files directly in outputs/
# Also remove old docking_results subdirectory if it exists (for backward compatibility)
if [ -d "${OUTPUT_DIR}" ]; then
    find "${OUTPUT_DIR}" -maxdepth 1 -type f \( -name "*.sdf" -o -name "*_log.txt" -o -name "*.txt" \) -delete 2>/dev/null || true
    # Remove old docking_results subdirectory if it exists
    if [ -d "${OUTPUT_DIR}/docking_results" ]; then
        rm -rf "${OUTPUT_DIR}/docking_results"
        echo "Removed old docking_results subdirectory"
    fi
    echo "Cleaned up stray files in ${OUTPUT_DIR}"
fi

# ============================================================================
# Copy files from source directories to input/ directory
# ============================================================================
# Copy PDB files from protein_extraction/outputs/ to input/
PROTEIN_EXTRACTION_DIR="${SCRIPT_DIR}/../protein_extraction/outputs"
if [ -d "${PROTEIN_EXTRACTION_DIR}" ]; then
    for pdb_file in "${PROTEIN_EXTRACTION_DIR}"/*.pdb; do
        if [ -f "${pdb_file}" ]; then
            cp "${pdb_file}" "${INPUT_DIR}/"
            echo "✓ Copied $(basename "${pdb_file}") from protein_extraction/outputs/"
        fi
    done
fi

# Copy config.txt from docking_setup/outputs/ to input/
DOCKING_SETUP_DIR="${SCRIPT_DIR}/../docking_setup/outputs"
if [ -f "${DOCKING_SETUP_DIR}/config.txt" ]; then
    cp "${DOCKING_SETUP_DIR}/config.txt" "${INPUT_DIR}/"
    echo "✓ Copied config.txt from docking_setup/outputs/"
fi

# Copy real_ligand.sdf from docking_setup/outputs/ to input/ (optional)
if [ -f "${DOCKING_SETUP_DIR}/real_ligand.sdf" ]; then
    cp "${DOCKING_SETUP_DIR}/real_ligand.sdf" "${INPUT_DIR}/"
    echo "✓ Copied real_ligand.sdf from docking_setup/outputs/"
fi

# Copy selected_compounds directory from library_construction/outputs/selected_compounds/ to input/
LIBRARY_CONSTRUCTION_DIR="${SCRIPT_DIR}/../library_construction/outputs/selected_compounds"
if [ -d "${LIBRARY_CONSTRUCTION_DIR}" ]; then
    # Remove existing directory if it exists (to avoid "Directory not empty" error)
    if [ -d "${INPUT_DIR}/selected_compounds" ]; then
        rm -rf "${INPUT_DIR}/selected_compounds"
    fi
    cp -r "${LIBRARY_CONSTRUCTION_DIR}" "${INPUT_DIR}/"
    echo "✓ Copied selected_compounds directory from library_construction/outputs/"
fi

# ============================================================================
# Move files from Silva mounts to input/ directory (fallback for other files only)
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

# Move config.txt from root directory (Silva mounts) to input/
if [ -f "${SCRIPT_DIR}/config.txt" ]; then
    mv "${SCRIPT_DIR}/config.txt" "${INPUT_DIR}/"
    echo "✓ Moved config.txt to input/"
fi

# Move real_ligand.sdf from root directory (Silva mounts) to input/
if [ -f "${SCRIPT_DIR}/real_ligand.sdf" ]; then
    mv "${SCRIPT_DIR}/real_ligand.sdf" "${INPUT_DIR}/"
    echo "✓ Moved real_ligand.sdf to input/"
fi

# Remove selected_compounds from root directory if Silva mounted it there
# (Silva may mount outputs/selected_compounds as root/selected_compounds)
if [ -d "${SCRIPT_DIR}/selected_compounds" ]; then
    rm -rf "${SCRIPT_DIR}/selected_compounds"
    echo "✓ Removed selected_compounds from root directory (Silva mount)"
fi

# ============================================================================
# Run nodes
# ============================================================================
# Change to script directory to ensure relative paths work correctly
cd "${SCRIPT_DIR}"
python3 "${SCRIPT_DIR}/smina_screening.py"

