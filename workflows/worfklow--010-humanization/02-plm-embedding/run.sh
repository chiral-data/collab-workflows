#!/usr/bin/env bash
# run2.sh – Node 2: PLM Embedding Generation (OPTIONAL)
# =======================================================
# Generates ProtT5 embeddings for each antibody pair from Node 1.
# Skip this node entirely if you only want plain ABB3 (no language model).
#
# Environment variables (all optional – defaults shown):
#   INPUTS_DIR       Node 1 output directory       (default: results/node1)
#   OUTPUTS_DIR      output directory               (default: results/node2)
#   PRECOMPUTED_DIR  directory with cached .pt files (default: unset)
#   DEVICE           cpu | cuda | auto              (default: auto)
#
# Usage:
#   ./run2.sh
#   PRECOMPUTED_DIR=data/embeddings ./run2.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${WORKFLOW_ROOT}"

INPUTS_DIR="${INPUTS_DIR:-results/node1}"
OUTPUTS_DIR="${OUTPUTS_DIR:-results/node2}"
PRECOMPUTED_DIR="${PRECOMPUTED_DIR:-}"
DEVICE="${DEVICE:-auto}"

echo "[Node 2] ============================================"
echo "[Node 2] Starting: $(date)"
echo "[Node 2] INPUTS_DIR      = ${INPUTS_DIR}"
echo "[Node 2] OUTPUTS_DIR     = ${OUTPUTS_DIR}"
echo "[Node 2] PRECOMPUTED_DIR = ${PRECOMPUTED_DIR:-<not set, will compute live>}"
echo "[Node 2] DEVICE          = ${DEVICE}"
echo "[Node 2] python          = $(command -v python || echo 'NOT FOUND')"
python --version || true
echo "[Node 2] ============================================"

if [[ ! -d "${INPUTS_DIR}" ]]; then
  echo "[Node 2] ERROR: inputs directory does not exist: ${INPUTS_DIR}"
  exit 2
fi

EXTRA_ARGS=""
if [[ -n "${PRECOMPUTED_DIR}" ]]; then
  EXTRA_ARGS="--precomputed_dir ${PRECOMPUTED_DIR}"
fi

python -u node2/node2.py \
  --inputs  "${INPUTS_DIR}"  \
  --outputs "${OUTPUTS_DIR}" \
  --device  "${DEVICE}"      \
  ${EXTRA_ARGS}

echo "[Node 2] Completed: $(date)"
echo "[Node 2] Output files:"
ls -1 "${OUTPUTS_DIR}"/*.pt 2>/dev/null | head -n 20 || true