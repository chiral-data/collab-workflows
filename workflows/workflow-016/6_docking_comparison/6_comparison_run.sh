#!/bin/bash
set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo -e "${GREEN}=== Node 6: Comparison with wet‑lab structure ===${NC}"

ROOT_DIR="$(dirname "$PWD")"

# Define paths
CHAIN_INFO_FILE="$(find inputs -type f -name chain_info.json | head -n 1)"
CHAIN_DATA_FILE="$(find inputs -type f -name data.json | grep -m1 '1_complex_splitting\|1_complex_split\|node1' || true)"
ORIG_COMPLEX="$(find inputs -type f -name original_complex.pdb | head -n 1)"
PRED_POSE="$(find inputs -type f -name rank1.pdb | head -n 1)"
OUTPUT_DIR="outputs"

if [ -z "$CHAIN_INFO_FILE" ] && [ -f "$ROOT_DIR/1_complex_splitting/outputs/chain_info.json" ]; then
    CHAIN_INFO_FILE="$ROOT_DIR/1_complex_splitting/outputs/chain_info.json"
fi
if [ -z "$CHAIN_DATA_FILE" ] && [ -f "$ROOT_DIR/1_complex_splitting/outputs/data.json" ]; then
    CHAIN_DATA_FILE="$ROOT_DIR/1_complex_splitting/outputs/data.json"
fi
if [ -z "$ORIG_COMPLEX" ] && [ -f "$ROOT_DIR/1_complex_splitting/outputs/original_complex.pdb" ]; then
    ORIG_COMPLEX="$ROOT_DIR/1_complex_splitting/outputs/original_complex.pdb"
fi
if [ -z "$PRED_POSE" ] && [ -f "$ROOT_DIR/4_diffdock_inference/outputs/rank1.pdb" ]; then
    PRED_POSE="$ROOT_DIR/4_diffdock_inference/outputs/rank1.pdb"
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Validate input files exist
if [ ! -f "$CHAIN_INFO_FILE" ]; then
    echo -e "${RED}Error: Chain info not found at $CHAIN_INFO_FILE${NC}"
    echo "Please run Node 1 first to generate chain information."
    exit 1
fi

if [ ! -f "$ORIG_COMPLEX" ]; then
    echo -e "${YELLOW}Warning: Original complex not found at $ORIG_COMPLEX${NC}"
    echo "Comparison cannot be performed. Skipping Node 6."
    exit 0
fi

if [ ! -f "$PRED_POSE" ]; then
    echo -e "${RED}Error: Predicted pose not found at $PRED_POSE${NC}"
    echo "Please run Node 4 first to generate docking poses."
    exit 1
fi

# Extract chain information
echo -e "${YELLOW}Reading chain information${NC}"
export CHAIN_DATA_FILE
export CHAIN_INFO_FILE
AB_CHAINS=$(python3 - <<'PY'
import json
import os
from pathlib import Path

data_path = Path(os.environ.get('CHAIN_DATA_FILE', ''))
chain_path = Path(os.environ.get('CHAIN_INFO_FILE', ''))

ab = []
if data_path.is_file():
    d = json.load(open(data_path))
    ab = d.get('antibody', {}).get('chains', [])

if not ab and chain_path.is_file():
    d = json.load(open(chain_path))
    info = d.get('chain_info', {})
    ab = [cid for cid, meta in info.items() if meta.get('type') == 'antibody']

print(','.join(ab))
PY
)
AG_CHAINS=$(python3 - <<'PY'
import json
import os
from pathlib import Path

data_path = Path(os.environ.get('CHAIN_DATA_FILE', ''))
chain_path = Path(os.environ.get('CHAIN_INFO_FILE', ''))

ag = []
if data_path.is_file():
    d = json.load(open(data_path))
    ag = d.get('antigen', {}).get('chains', [])

if not ag and chain_path.is_file():
    d = json.load(open(chain_path))
    info = d.get('chain_info', {})
    ag = [cid for cid, meta in info.items() if meta.get('type') == 'antigen']

print(','.join(ag))
PY
)

if [ -z "$AB_CHAINS" ] || [ -z "$AG_CHAINS" ]; then
    echo -e "${RED}Error: Could not extract chain information${NC}"
    exit 1
fi

echo -e "${GREEN}Antibody chains: $AB_CHAINS${NC}"
echo -e "${GREEN}Antigen chains: $AG_CHAINS${NC}"

# Run comparison science script
echo -e "${YELLOW}Running comparison analysis...${NC}"
python3 11_comparison_science.py \
    --original_complex "$ORIG_COMPLEX" \
    --pred_pose "$PRED_POSE" \
    --ab_chains "$AB_CHAINS" \
    --ag_chains "$AG_CHAINS" \
    --output_dir "$OUTPUT_DIR"

# Generate HTML report
echo -e "${YELLOW}Generating HTML report...${NC}"
python3 12_comparison_html.py \
    --data_json "$OUTPUT_DIR/data.json" \
    --output_html "$OUTPUT_DIR/report.html"

echo -e "${GREEN}✓ Node 6 finished successfully${NC}"
echo -e "${GREEN}Outputs written to: $OUTPUT_DIR/${NC}"