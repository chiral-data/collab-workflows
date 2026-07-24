#!/bin/bash
set -euo pipefail

echo "=== Apply Glass-Fiber Correction (Halpin-Tsai) ==="
mkdir -p outputs
python3 apply_correction.py
echo ""
echo "Done — outputs: corrected_properties.json"
