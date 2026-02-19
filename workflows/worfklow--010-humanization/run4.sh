#!/usr/bin/env bash
# run4.sh – Node 4: PDB Export & HTML Visualization
# ===================================================
# Copies predicted PDB files from Node 3 and generates
# a self-contained interactive HTML report with py3Dmol viewers.
#
# Environment variables (all optional – defaults shown):
#   INPUTS_DIR     Node 3 output directory     (default: results/node3)
#   OUTPUTS_DIR    output directory             (default: results/node4)
#   REPORT_TITLE   title shown in the report    (default: ABB3 Structure Predictions)
#
# Usage:
#   ./run4.sh
#   REPORT_TITLE="My Antibody Run" ./run4.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKFLOW_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${WORKFLOW_ROOT}"

INPUTS_DIR="${INPUTS_DIR:-results/node3}"
OUTPUTS_DIR="${OUTPUTS_DIR:-results/node4}"
REPORT_TITLE="${REPORT_TITLE:-ABB3 Structure Predictions}"

echo "[Node 4] ============================================"
echo "[Node 4] Starting: $(date)"
echo "[Node 4] INPUTS_DIR    = ${INPUTS_DIR}"
echo "[Node 4] OUTPUTS_DIR   = ${OUTPUTS_DIR}"
echo "[Node 4] REPORT_TITLE  = ${REPORT_TITLE}"
echo "[Node 4] python        = $(command -v python || echo 'NOT FOUND')"
python --version || true
echo "[Node 4] ============================================"

if [[ ! -d "${INPUTS_DIR}" ]]; then
  echo "[Node 4] ERROR: inputs directory does not exist: ${INPUTS_DIR}"
  exit 2
fi

echo "[Node 4] PDB files found: $(ls -1 "${INPUTS_DIR}"/*.pdb 2>/dev/null | wc -l)"
ls -1 "${INPUTS_DIR}"/*.pdb 2>/dev/null | head -n 10 || true

python -u node4/node4.py \
  --inputs  "${INPUTS_DIR}"   \
  --outputs "${OUTPUTS_DIR}"  \
  --title   "${REPORT_TITLE}"

echo "[Node 4] Completed: $(date)"
echo "[Node 4] Outputs:"
ls -1 "${OUTPUTS_DIR}" | head -n 30 || true