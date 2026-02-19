#!/usr/bin/env bash
# run_all.sh – Full ABB3 Workflow
# =================================
# Runs all four nodes in sequence:
#
#   Node 1  Read & validate FASTA → paired .pt files
#   Node 2  (Optional) ProtT5 PLM embeddings
#   Node 3  ABB3 structure prediction → PDB files
#   Node 4  HTML visualization report
#
# ---------------------------------------------------------------
# REQUIRED environment variables:
#   HEAVY_FASTA        path to heavy-chain FASTA
#   LIGHT_FASTA        path to light-chain FASTA
#   CHECKPOINT_PATH    path to ABB3 model checkpoint (.ckpt)
#
# OPTIONAL environment variables:
#   USE_PLM            set to "1" to enable ABB3-LM via Node 2
#   PRECOMPUTED_DIR    directory with pre-computed PLM .pt files
#   DEVICE             cpu | cuda | auto  (default: auto)
#   REPORT_TITLE       title for the HTML report
#   RESULTS_ROOT       root directory for all results (default: results)
# ---------------------------------------------------------------
#
# Quick start:
#   HEAVY_FASTA=data/heavy.fasta \
#   LIGHT_FASTA=data/light.fasta \
#   CHECKPOINT_PATH=/path/to/best_second_stage.ckpt \
#   ./run_all.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

# ---- Config ---------------------------------------------------------------
HEAVY_FASTA="${HEAVY_FASTA:-data/heavy.fasta}"
LIGHT_FASTA="${LIGHT_FASTA:-data/light.fasta}"
CHECKPOINT_PATH="${CHECKPOINT_PATH:-}"
USE_PLM="${USE_PLM:-0}"
PRECOMPUTED_DIR="${PRECOMPUTED_DIR:-}"
DEVICE="${DEVICE:-auto}"
REPORT_TITLE="${REPORT_TITLE:-ABB3 Structure Predictions}"
RESULTS_ROOT="${RESULTS_ROOT:-results}"

NODE1_OUT="${RESULTS_ROOT}/node1"
NODE2_OUT="${RESULTS_ROOT}/node2"
NODE3_OUT="${RESULTS_ROOT}/node3"
NODE4_OUT="${RESULTS_ROOT}/node4"

echo "======================================================="
echo "  ABB3 Workflow"
echo "======================================================="
echo "  HEAVY_FASTA      = ${HEAVY_FASTA}"
echo "  LIGHT_FASTA      = ${LIGHT_FASTA}"
echo "  CHECKPOINT_PATH  = ${CHECKPOINT_PATH:-<not set>}"
echo "  USE_PLM          = ${USE_PLM}"
echo "  DEVICE           = ${DEVICE}"
echo "  RESULTS_ROOT     = ${RESULTS_ROOT}"
echo "  Started          : $(date)"
echo "======================================================="

if [[ -z "${CHECKPOINT_PATH}" ]]; then
  echo "ERROR: CHECKPOINT_PATH is required."
  echo "       Set it via: export CHECKPOINT_PATH=/path/to/model.ckpt"
  exit 2
fi

# ---- Node 1 ---------------------------------------------------------------
echo ""
echo "--- Node 1: Input Preparation ---"
HEAVY_FASTA="${HEAVY_FASTA}" \
LIGHT_FASTA="${LIGHT_FASTA}" \
OUT_DIR="${NODE1_OUT}" \
bash run1.sh

# ---- Node 2 (optional) ----------------------------------------------------
if [[ "${USE_PLM}" == "1" ]]; then
  echo ""
  echo "--- Node 2: PLM Embedding Generation ---"
  EXTRA=""
  if [[ -n "${PRECOMPUTED_DIR}" ]]; then
    EXTRA="PRECOMPUTED_DIR=${PRECOMPUTED_DIR}"
  fi
  INPUTS_DIR="${NODE1_OUT}" \
  OUTPUTS_DIR="${NODE2_OUT}" \
  DEVICE="${DEVICE}" \
  ${EXTRA} \
  bash run2.sh
  NODE3_INPUT="${NODE2_OUT}"
else
  echo ""
  echo "--- Node 2: Skipped (USE_PLM != 1) ---"
  NODE3_INPUT="${NODE1_OUT}"
fi

# ---- Node 3 ---------------------------------------------------------------
echo ""
echo "--- Node 3: Structure Prediction ---"
INPUTS_DIR="${NODE3_INPUT}" \
OUTPUTS_DIR="${NODE3_OUT}" \
CHECKPOINT_PATH="${CHECKPOINT_PATH}" \
DEVICE="${DEVICE}" \
bash run3.sh

# ---- Node 4 ---------------------------------------------------------------
echo ""
echo "--- Node 4: Visualization ---"
INPUTS_DIR="${NODE3_OUT}" \
OUTPUTS_DIR="${NODE4_OUT}" \
REPORT_TITLE="${REPORT_TITLE}" \
bash run4.sh

# ---- Summary ---------------------------------------------------------------
echo ""
echo "======================================================="
echo "  Workflow complete: $(date)"
echo "  Results:"
echo "    Node 1 inputs  : ${NODE1_OUT}/"
[[ "${USE_PLM}" == "1" ]] && echo "    Node 2 embeddings: ${NODE2_OUT}/"
echo "    Node 3 PDBs    : ${NODE3_OUT}/"
echo "    Node 4 report  : ${NODE4_OUT}/report.html"
echo "======================================================="