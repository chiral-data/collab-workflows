#!/bin/bash
set -e

echo "Starting Node 02: P2Rank Pocket Detection"

cp inputs/* . 2>/dev/null || true

# Select the Boltz-2 model with the highest confidence score.
# confidence_model_N.json -> model_N -> *_model_N.cif
BEST_MODEL_ID=$(python3 - <<'EOF'
import json, glob, sys

jsons = sorted(glob.glob("confidence_*.json"))
if not jsons:
    print("model_0")
    sys.exit(0)

best_score = -1
best_id = "model_0"
for jf in jsons:
    try:
        with open(jf) as fh:
            d = json.load(fh)
        score = float(d.get("confidence_score", 0))
        if score > best_score:
            best_score = score
            # confidence_model_0.json -> model_0
            best_id = jf.replace("confidence_", "").replace(".json", "")
    except Exception:
        pass

print(best_id)
EOF
)

BEST_CIF=$(ls *_${BEST_MODEL_ID}.cif 2>/dev/null | head -1)
if [ -z "$BEST_CIF" ]; then
    BEST_CIF=$(ls *.cif 2>/dev/null | head -1)
    BEST_MODEL_ID=$(echo "$BEST_CIF" | grep -oP 'model_\d+' || echo "model_0")
fi
if [ -z "$BEST_CIF" ]; then
    echo "ERROR: No .cif structure file found in inputs/"
    exit 1
fi

echo "Selected model: $BEST_MODEL_ID ($BEST_CIF)"

# Run P2Rank with AlphaFold pLDDT profile.
# -c alphafold: alphafold config drops B-factor from the feature set
mkdir -p ./p2rank_out
prank predict -f "$BEST_CIF" -c alphafold -o ./p2rank_out

# P2Rank writes <input_filename>_predictions.csv and <input_filename>_residues.csv
# into the output dir. Rename to stable names for downstream nodes.
PRED_FILE=$(ls ./p2rank_out/*_predictions.csv 2>/dev/null | head -1)
RESID_FILE=$(ls ./p2rank_out/*_residues.csv 2>/dev/null | head -1)

if [ -z "$PRED_FILE" ]; then
    echo "ERROR: P2Rank produced no predictions.csv — check prank output above"
    exit 1
fi

mkdir -p ./outputs
cp "$PRED_FILE"  ./outputs/predictions.csv
cp "$RESID_FILE" ./outputs/residues.csv    2>/dev/null || true
cp "$BEST_CIF"   ./outputs/selected_structure.cif
echo "$BEST_MODEL_ID" > ./outputs/selected_model_id.txt

NPOCKETS=$(tail -n +2 ./outputs/predictions.csv | wc -l)
echo "P2Rank identified $NPOCKETS pocket(s) in $BEST_CIF"

echo "Node 02 completed"
