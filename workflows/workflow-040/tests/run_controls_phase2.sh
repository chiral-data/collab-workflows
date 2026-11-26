#!/bin/bash
# Tier 2, Phase 2: the controls that can falsify INTERFACE scoring rather than the plumbing.
#
#   C8  1BRS barnase on barstar, native sequence     must score well   (positive)
#   C9  the same complex, binder sequence scrambled  must score badly  (negative)
#
# One run, one target, both binders. The separation between them is read off a single
# report rather than compared across runs, for the same reason Phase 1 does it that way:
# two runs differ by more than the thing under test.
#
# What this establishes, and what it does not
# -------------------------------------------
# It establishes SEPARATION: whether 040's interface metrics tell a real binder from the
# same backbone carrying a shuffled sequence. That is the question Gate B Phase 2 asks,
# and §8.6 is why it has to be asked at all -- iptm is a number whatever went in.
#
# It does NOT establish calibration. The target alignment here is SYNTHETIC: homologs built
# from the target's own sequence by point substitution (tests/make_target_msa.py), because
# this run is offline and a real target MSA needs api.colabfold.com. A tiny synthetic
# alignment depresses the absolute numbers -- §11.1 measured that directly, where a
# synthetic a3m scored WORSE than no alignment at all. So no threshold may be set from the
# absolute values this run produces. Only the gap between C8 and C9 means anything, and
# the report says so in its own verdict rather than leaving it to be remembered.
#
# C7 is the control that needs the real alignment, and it is the reason §4's decision is
# still only half measured. It is NOT run here: it needs one fetch from api.colabfold.com,
# which is outward-facing and is the user's call, not this script's. Run it with
# --allow-egress once that call is made.
#
# Why no single_sequence arm
# --------------------------
# Node 03 refuses msa_mode=single_sequence on the complex path by name, because colabfold
# clamps max_seq to 1 and discards the alignment while still writing every artifact that
# looks like proof it was used (§11.1, tests/probe_msa_clamp.sh). That refusal is the
# behaviour under test, so this script shows it refusing rather than working around it.
#
# Needs: docker, the colabfold:2026_10_08 image, a GPU, curl, python3.
#
# Usage:  bash tests/run_controls_phase2.sh [OUT_DIR] [DOWNLOADS_DIR] [--allow-egress]
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WF="$(dirname "$HERE")"
WORKFLOWS="$(dirname "$WF")"
IMAGE="colabfold:2026_10_08"

# WHICH PAIR. Default is C8/C9, the pair this script was written for, so every command
# recorded in the docs and in the evidence directories keeps working unchanged.
#
# --pair c13 selects the designed-binder control. The difference that matters is not the
# entry, it is that C13's BINDER is a synthetic construct rather than a native chain, so
# both chains can be modelled under the flags section 4 mandates: the target from its
# alignment, the binder from nothing. C8/C9's binder is native barnase and cannot be,
# which is why they came back NOT INTERPRETABLE twice.
#
# Parameterised rather than copied into a second script: the staging, the egress scoping,
# the n3-refuse step and the free-space floor are identical, and a 160-line copy would be
# two things to keep in step (the docs-with-code rule, and 8.15).
PAIR="c8"
ALLOW_EGRESS=0
ARGS=()
while [ "$#" -gt 0 ]; do
  case "$1" in
    --allow-egress) ALLOW_EGRESS=1 ;;
    --pair) PAIR="${2:?--pair needs a value: c8 or c13}"; shift ;;
    --pair=*) PAIR="${1#--pair=}" ;;
    *) ARGS+=("$1") ;;
  esac
  shift
done

case "$PAIR" in
  c8)
    ENTRY="1BRS"; BINDER_CHAIN="A"; TARGET_CHAIN="D"
    POS="C8"; NEG="C9"
    PAIR_NOTE="barnase (auth A) on barstar (auth D). The binder is NATIVE, so it cannot be modelled without an alignment and section 4 gives it none: this pair is NOT INTERPRETABLE by construction and is kept as the measured record of that."
    ;;
  c13)
    ENTRY="7XYQ"; BINDER_CHAIN="B"; TARGET_CHAIN="A"
    POS="C13"; NEG="C14"
    PAIR_NOTE="DBL1_03 (auth B, a synthetic construct, 124 res) on human PD-L1 (auth A, 203 res, 4-residue gap at 165->170). The binder is DESIGNED, so both chains can be modelled under section 4's flags. Pre-registered in tests/reference/tier2-phase2-controls.md BEFORE this code existed."
    ;;
  *)
    echo "run_controls_phase2: --pair must be c8 or c13, got '$PAIR'." >&2
    echo "Guessing would mean guessing which experiment you meant." >&2
    exit 1
    ;;
esac
OUT="${ARGS[0]:-$HOME/wf040-controls-phase2/$(date +%Y%m%d-%H%M%S)}"
DL="${ARGS[1]:-$HOME/wf040-downloads}"

# The GPU host's own ~/.bashrc exports a DOCKER_HOST that does not exist, so a non-login
# shell sees no daemon at all rather than an error (README, "Working on the GPU host").
if [ -S "/run/user/$(id -u)/docker.sock" ]; then
  export DOCKER_HOST="unix:///run/user/$(id -u)/docker.sock"
fi

say() { printf '\n=== %s\n' "$*"; }
free_gb() { df -BG --output=avail "$HOME" | tail -1 | tr -dc '0-9'; }

BEFORE_GB="$(free_gb)"
say "free space before: ${BEFORE_GB}G"
# Never delete anything on this host, so the only lever is refusing to start. 8G is the
# Phase 1 floor and this run writes less than that one does.
if [ "$BEFORE_GB" -lt 8 ]; then
  echo "REFUSING: ${BEFORE_GB}G free, under the 8G floor. Nothing was written." >&2
  exit 1
fi

mkdir -p "$OUT" "$DL"
say "out: $OUT"

# ---------------------------------------------------------------- the entry
# Its own directory, nothing else in it, and every tool that reads it runs with -I: a file
# fetched from the PDB is untrusted input and must not be able to put a module on the
# import path. (NOT -I for this repo's own node scripts -- it implies -s and hides the
# user site-packages where biopython lives.)
if [ ! -s "$DL/$ENTRY.cif" ]; then
  say "fetching $ENTRY"
  curl -fsSL "https://files.rcsb.org/download/$ENTRY.cif" -o "$DL/$ENTRY.cif" || {
    echo "REFUSING: could not fetch $ENTRY." >&2; exit 1; }
fi
sha256sum "$DL/$ENTRY.cif" | tee "$OUT/entry-sha256.txt"

# Probed, then chosen from the probe. 037 lost five of 1BRS's six chains to "chain A"
# (53d1d460), so the chain ids are read off the probe rather than assumed.
say "probing $ENTRY"
python3 -I "$HERE/make_controls.py" --list "$DL/$ENTRY.cif" | tee "$OUT/probe-$ENTRY.txt"

# ---------------------------------------------------------------- build the two
# Barnase (auth A) as the binder, barstar (auth D) as the target: a real, tight,
# well-characterised pair, and the same chains Phase 1's C2 used so the provenance is
# already probed. C9 is the SAME backbone with the SAME residues in a shuffled order, so
# composition, length, numbering and every coordinate are identical between the two and
# the only thing that differs is the sequence. That is what makes the comparison a control
# rather than two experiments.
say "pair $PAIR: $PAIR_NOTE"
say "building $POS (native) and $NEG (scrambled) on one target"
mkdir -p "$OUT/g"
printf '%s\n' "$PAIR" "$ENTRY" "$BINDER_CHAIN" "$TARGET_CHAIN" "$POS" "$NEG" \
  > "$OUT/pair.txt"
python3 -I "$HERE/make_controls.py" "$DL/$ENTRY.cif" "$OUT/g" \
  --name "$POS" --binder "$BINDER_CHAIN" --target "$TARGET_CHAIN" \
  --sequence native || exit 1
python3 -I "$HERE/make_controls.py" "$DL/$ENTRY.cif" "$OUT/g" \
  --name "$NEG" --binder "$BINDER_CHAIN" --target "$TARGET_CHAIN" \
  --sequence scrambled --seed 0 || exit 1
python3 -I "$HERE/probe_design.py" "$OUT/g" | tee "$OUT/probe-controls.txt"

# ---------------------------------------------------------------- run the nodes
run_node() {           # run_node <step-dir> <script> [extra docker args...]
  local step="$1" script="$2"; shift 2
  mkdir -p "$OUT/$step/inputs" "$OUT/$step/outputs"
  docker run --rm "$NETWORK" \
    -e PARAM_SCORE_INTERFACE=true \
    -e PARAM_MSA_SOURCE="$MSA_SOURCE" \
    -v "$WORKFLOWS:$WORKFLOWS:ro" -v "$OUT:$OUT" -w "$OUT/$step" \
    "$@" "$IMAGE" python "$script" 2>&1 | tee "$OUT/$step.log"
  return "${PIPESTATUS[0]}"
}

NETWORK="--network=none"
MSA_SOURCE="file"
if [ "$ALLOW_EGRESS" = "1" ]; then
  MSA_SOURCE="server"
  say "EGRESS ALLOWED: node 01b will align the TARGET against api.colabfold.com"
  echo "  One sequence goes out: the fused TARGET chain, public protein from the PDB."
  echo "  The designed binders ARE in 01b's input (they ride in design_inventory.json);"
  echo "  assert_no_designed_sequence() refuses the POST if one is in the payload, and"
  echo "  fixture case egress-refuses-a-designed-sequence shows it refusing."
fi

say "037 node 04: the real design_metrics.json producer"
mkdir -p "$OUT/m/inputs"
cp "$OUT"/g/*_generated.cif "$OUT/g/generation_manifest.json" "$OUT/m/inputs/"
run_node m "$WORKFLOWS/workflow-037/04-analyze/script.py"

say "040 node 01: ingest"
mkdir -p "$OUT/n1/inputs"
cp "$OUT"/g/*_generated.cif "$OUT/m/outputs/design_metrics.json" "$OUT/n1/inputs/"
run_node n1 "$WF/01-ingest-designs/script.py" || exit 1

# The target alignment. Built from the inventory node 01 just wrote, not hard-coded: 01b
# validates the query row against the target sequence it derives from the inventory, so a
# hard-coded a3m would be refused and the run would test the refusal instead of the path.
say "040 node 01b: align the target"
mkdir -p "$OUT/n1b/inputs"
cp "$OUT/n1/outputs/design_inventory.json" "$OUT/n1b/inputs/"
if [ "$MSA_SOURCE" = "file" ]; then
  python3 -I "$HERE/make_target_msa.py" \
    "$OUT/n1b/inputs/design_inventory.json" "$OUT/n1b/inputs/target_msa.a3m" 8 \
    | tee "$OUT/target-msa.txt"
fi
# Egress, when it is allowed, is allowed for THIS node only and for nothing else in the
# run. Every other node above and below stays on --network=none.
if [ "$MSA_SOURCE" = "server" ]; then NETWORK="--network=bridge"; fi
run_node n1b "$WF/01b-align-target/script.py" || exit 1
NETWORK="--network=none"

say "040 node 02: build BOTH prediction inputs"
mkdir -p "$OUT/n2/inputs"
cp "$OUT"/n1/outputs/* "$OUT"/n1b/outputs/* "$OUT/n2/inputs/"
run_node n2 "$WF/02-build-prediction-inputs/script.py" || exit 1
ls "$OUT/n2/outputs" | sed 's/^/  /'

# The refusal, shown refusing, before the real run spends any GPU time on it. A guard that
# cannot be shown refusing is not a guard (§8.13), and this one is the whole reason the
# complex path has its own msa_mode.
say "node 03 must REFUSE single_sequence on the complex path"
mkdir -p "$OUT/n3-refuse/inputs"
cp "$OUT"/n2/outputs/* "$OUT/n3-refuse/inputs/"
docker run --rm --network=none \
  -e PARAM_SCORE_INTERFACE=true -e PARAM_MSA_MODE_CPLX=single_sequence \
  -v "$WORKFLOWS:$WORKFLOWS:ro" -v "$OUT:$OUT" -w "$OUT/n3-refuse" \
  "$IMAGE" python "$WF/03-predict/script.py" --dry-run \
  > "$OUT/n3-refuse.log" 2>&1
echo "  exit=$? (non-zero is the pass)" | tee -a "$OUT/n3-refuse.log"
tail -3 "$OUT/n3-refuse.log" | sed 's/^/  /'

say "040 node 03: REAL multimer prediction on the GPU, offline"
mkdir -p "$OUT/n3/inputs"
cp "$OUT"/n2/outputs/* "$OUT/n3/inputs/"
nvidia-smi --query-gpu=name,memory.used,memory.total --format=csv,noheader \
  | tee "$OUT/gpu-before.txt"
run_node n3 "$WF/03-predict/script.py" --gpus all
nvidia-smi --query-gpu=name,memory.used,memory.total --format=csv,noheader \
  | tee "$OUT/gpu-after.txt"

say "040 node 04: score, with the interface overlaid"
mkdir -p "$OUT/n4/inputs"
cp "$OUT"/n1/outputs/* "$OUT"/n3/outputs/* "$OUT/n4/inputs/"
run_node n4 "$WF/04-score/script.py" || exit 1

say "040 node 05: filter and report"
mkdir -p "$OUT/n5/inputs"
cp "$OUT"/n4/outputs/* "$OUT/n5/inputs/"
run_node n5 "$WF/05-filter-and-report/script.py"

# ---------------------------------------------------------------- the verdict
say "free space after: $(free_gb)G (was ${BEFORE_GB}G)"
python3 -I "$HERE/report_controls_phase2.py" "$OUT" | tee "$OUT/VERDICT.txt"
