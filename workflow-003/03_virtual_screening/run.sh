#!/usr/bin/env bash
set -euo pipefail

# Input and output directories
INDIR=${INDIR:-/work/in}
OUTDIR=${OUTDIR:-/work/out}

mkdir -p "$INDIR" "$OUTDIR"

echo "Inputs: $INDIR"
echo "Outputs: $OUTDIR"

# Run docking pipeline
python /work/03_virtual_screening.py

# Optional: automatically run ranking step
echo "Running ranking step (rank_vina.py)..."
python /work/rank_vina.py

echo "=== Step 2: Running docking ==="
python vina_screen.py --center "$BOX_CENTER" --size "$BOX_SIZE"

echo "=== Step 3: Ranking ligands ==="
python rank_vina.py

echo "=== All steps completed successfully! ==="