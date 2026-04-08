#!/bin/bash
set -e
echo "=== Node 4: DiffDock-PP inference ==="

ROOT_DIR="$(dirname "$PWD")"

RECEPTOR_PDB="$(find inputs -type f -name processed_antibody.pdb | head -n 1)"
LIGAND_PDB="$(find inputs -type f -name processed_antigen.pdb | head -n 1)"
RECEPTOR_FEATURES="$(find inputs -type f -name antibody_features.pt | head -n 1)"
LIGAND_FEATURES="$(find inputs -type f -name antigen_features.pt | head -n 1)"

if [ -z "$RECEPTOR_PDB" ] && [ -f "$ROOT_DIR/2_diffdock_prep/outputs/processed_antibody.pdb" ]; then
    RECEPTOR_PDB="$ROOT_DIR/2_diffdock_prep/outputs/processed_antibody.pdb"
fi
if [ -z "$LIGAND_PDB" ] && [ -f "$ROOT_DIR/2_diffdock_prep/outputs/processed_antigen.pdb" ]; then
    LIGAND_PDB="$ROOT_DIR/2_diffdock_prep/outputs/processed_antigen.pdb"
fi
if [ -z "$RECEPTOR_FEATURES" ] && [ -f "$ROOT_DIR/3_diffdock_features/outputs/antibody_features.pt" ]; then
    RECEPTOR_FEATURES="$ROOT_DIR/3_diffdock_features/outputs/antibody_features.pt"
fi
if [ -z "$LIGAND_FEATURES" ] && [ -f "$ROOT_DIR/3_diffdock_features/outputs/antigen_features.pt" ]; then
    LIGAND_FEATURES="$ROOT_DIR/3_diffdock_features/outputs/antigen_features.pt"
fi

if [ -z "$RECEPTOR_PDB" ] || [ -z "$LIGAND_PDB" ] || [ -z "$RECEPTOR_FEATURES" ] || [ -z "$LIGAND_FEATURES" ]; then
    echo "Error: missing required inputs for Node 4"
    exit 1
fi

MOCK_FLAG=()
if [[ "${ALLOW_MOCK_FALLBACK:-0}" == "1" ]]; then
    MOCK_FLAG+=(--allow_mock_fallback)
fi

GPU_FLAG=()
if [[ "${USE_GPU:-1}" == "1" || "${USE_GPU:-1}" == "true" ]]; then
    GPU_FLAG+=(--use_gpu)
elif [[ "${USE_GPU:-1}" == "auto" ]]; then
    if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi >/dev/null 2>&1; then
        GPU_FLAG+=(--use_gpu)
    fi
fi

python 7_diffdock_inference_science.py \
    --receptor_pdb "$RECEPTOR_PDB" \
    --ligand_pdb "$LIGAND_PDB" \
    --receptor_features "$RECEPTOR_FEATURES" \
    --ligand_features "$LIGAND_FEATURES" \
    --diffdock_path "${DIFFDOCK_PP_PATH:-}" \
    --num_samples "${PARAM_NUM_SAMPLES:-10}" \
    --inference_steps "${PARAM_INFERENCE_STEPS:-20}" \
    --output_dir "outputs" \
    "${GPU_FLAG[@]}" \
    "${MOCK_FLAG[@]}"
python 8_diffdock_inference_html.py --data_json "outputs/data.json" --output_html "outputs/report.html"
echo "Node 4 finished"