#!/usr/bin/env bash
# run3.sh – Node 3: ABB3 Structure Prediction
# =============================================
# Runs ABB3 (or ABB3-LM if inputs contain PLM embeddings from Node 2)
# and writes one PDB file per antibody pair.
#
# Environment variables (all optional – defaults shown):
#   INPUTS_DIR       Node 1 or Node 2 output dir   (default: results/node1)
#   OUTPUTS_DIR      output directory               (default: results/node3)
#   CHECKPOINT_PATH  path to model .ckpt file       (REQUIRED)
#   DEVICE           cpu | cuda | auto              (default: auto)
#
# Usage (plain ABB3, reading from Node 1):
#   CHECKPOINT_PATH=/path/to/best_second_stage.ckpt ./run3.sh
#
# Usage (ABB3-LM, reading from Node 2):
#   INPUTS_DIR=results/node2 CHECKPOINT_PATH=/path/to/lm.ckpt ./run3.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${WORKFLOW_ROOT}"

INPUTS_DIR="${INPUTS_DIR:-results/node1}"
OUTPUTS_DIR="${OUTPUTS_DIR:-results/node3}"
CHECKPOINT_PATH="${CHECKPOINT_PATH:-${1:-}}"
DEVICE="${DEVICE:-auto}"

echo "[Node 3] ============================================"
echo "[Node 3] Starting: $(date)"
echo "[Node 3] INPUTS_DIR      = ${INPUTS_DIR}"
echo "[Node 3] OUTPUTS_DIR     = ${OUTPUTS_DIR}"
echo "[Node 3] CHECKPOINT_PATH = ${CHECKPOINT_PATH}"
echo "[Node 3] DEVICE          = ${DEVICE}"
echo "[Node 3] python          = $(command -v python || echo 'NOT FOUND')"
python --version || true
echo "[Node 3] ============================================"

if [[ -z "${CHECKPOINT_PATH}" ]]; then
  echo "[Node 3] ERROR: CHECKPOINT_PATH is required."
  echo "         Set it via env var or pass as first argument:"
  echo "         CHECKPOINT_PATH=/path/to/model.ckpt ./run3.sh"
  exit 2
fi

if [[ ! -f "${CHECKPOINT_PATH}" ]]; then
  echo "[Node 3] ERROR: checkpoint file not found: ${CHECKPOINT_PATH}"
  exit 2
fi

if [[ ! -d "${INPUTS_DIR}" ]]; then
  echo "[Node 3] ERROR: inputs directory does not exist: ${INPUTS_DIR}"
  exit 2
fi

echo "[Node 3] Input .pt files found: $(ls -1 "${INPUTS_DIR}"/*.pt 2>/dev/null | wc -l)"
ls -1 "${INPUTS_DIR}"/*.pt 2>/dev/null | head -n 10 || true

python -u node3/node3.py \
  --inputs     "${INPUTS_DIR}"     \
  --checkpoint "${CHECKPOINT_PATH}" \
  --outputs    "${OUTPUTS_DIR}"    \
  --device     "${DEVICE}"

echo "[Node 3] Completed: $(date)"
echo "[Node 3] Output PDB files:"
ls -1 "${OUTPUTS_DIR}"/*.pdb 2>/dev/null | head -n 20 || true