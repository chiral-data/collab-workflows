#!/bin/bash
# C10: what does ColabFold ACTUALLY do with a complex input whose ':' is missing?
#
# Node 03 refuses this input by name, and the Tier 1 case `cplx-input-missing-colon` shows
# it refusing before the predictor starts. So the guard is tested. What the guard protects
# against was not measured until this control ran.
#
# MEASURED, and it corrected the rationale this script was written with. The header used to
# say the hazard was a colon-less complex coming back as one fused chain WHILE STILL
# PRODUCING AN IPTM -- a plausible interface number from an input declaring no interface
# (§8.4, §8.6). Half of that is true and half of it was mine, asserted rather than
# measured, which is §8.17's error in the rationale of the very control sent to settle it.
#
#   true        the two sequences are silently fused into one 195-residue chain, and the
#               predictor does not complain at all
#   NOT true    no iptm is reported. ColabFold omits it for a single-chain prediction,
#               which is the sane behaviour: there is no interface to score
#
# So the guard is still necessary -- a chimera of two proteins is not the complex anyone
# asked for, and the failure is SILENT -- but the specific §8.6 trap is absent here. Node
# 03's own refusal message never made the iptm claim; it says only that the colon count is
# wrong and the file must have changed since node 02 verified it, which remains exactly
# right. The overstatement was in this file, and it is corrected rather than deleted.
#
# It cannot be measured through node 03, because node 03 correctly refuses to run it. So
# this control calls colabfold_batch DIRECTLY. That is the one place in this suite where
# bypassing a node is the right thing to do, and it is worth being explicit about why: the
# subject here is the predictor's behaviour, not the workflow's, and the workflow is
# deliberately built so that this input never reaches the predictor.
#
# TWO ARMS, ONE VARIABLE
#
#   A  the real cplx_ fasta, target:binder          expected 2 chains, 87 + 108
#   B  byte-identical except the ':' is deleted     the question
#
# Same model, same msa_mode, same seed handling, same container. The colon is the only
# difference, which is what makes this a control rather than two observations.
#
# --msa-mode single_sequence with NO a3m is honest here, unlike on the complex path where
# node 03 refuses it: §11.1's clamp matters only when an alignment is handed in and then
# discarded. Nothing is handed in. It also keeps both arms fast, and the chain count does
# not depend on MSA depth.
#
# Needs: docker, the colabfold:2026_10_08 image, a GPU, python3, and a completed Phase 2
# control run to take the real cplx_ fasta from.
#
# Usage:  bash tests/run_c10.sh <PHASE2_RUN_DIR> [OUT_DIR]
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
SRC="${1:?usage: run_c10.sh <PHASE2_RUN_DIR> [OUT_DIR]}"
OUT="${2:-$HOME/wf040-c10/$(date +%Y%m%d-%H%M%S)}"
IMAGE="colabfold:2026_10_08"
DESIGN="C8"

# The GPU host's own ~/.bashrc exports a DOCKER_HOST that does not exist, so a non-login
# shell sees no daemon at all rather than an error.
if [ -z "${DOCKER_HOST:-}" ] && [ -S "/run/user/$(id -u)/docker.sock" ]; then
  export DOCKER_HOST="unix:///run/user/$(id -u)/docker.sock"
fi

say() { printf '\n=== %s\n' "$*"; }
free_gb() { df --output=avail -BG / | tail -1 | tr -dc '0-9'; }

BEFORE_GB="$(free_gb)"
say "free space before: ${BEFORE_GB}G"
if [ "$BEFORE_GB" -lt 8 ]; then
  echo "run_c10: only ${BEFORE_GB}G free. Refusing to start: this host runs near full and" >&2
  echo "a prediction run that fills the disk is worse than one that never ran." >&2
  exit 1
fi

CPLX="$SRC/n2/outputs/cplx_$DESIGN.fasta"
[ -f "$CPLX" ] || { echo "run_c10: no $CPLX. Run tests/run_controls_phase2.sh first." >&2; exit 1; }

docker image inspect "$IMAGE" >/dev/null 2>&1 || {
  echo "run_c10: image $IMAGE not found. Build it from workflow-026/Dockerfile." >&2
  exit 1
}

mkdir -p "$OUT/A/in" "$OUT/A/out" "$OUT/B/in" "$OUT/B/out"
say "out: $OUT"

# ---------------------------------------------------------------- the two inputs
cp "$CPLX" "$OUT/A/in/$DESIGN.fasta"
# Arm B differs by exactly one deleted byte. Built with tr rather than retyped, so the
# two arms cannot drift in any other character (§8.14: a bundled example must be what it
# says it is).
tr -d ':' < "$CPLX" > "$OUT/B/in/$DESIGN.fasta"

python3 -I - "$OUT/A/in/$DESIGN.fasta" "$OUT/B/in/$DESIGN.fasta" \
  > "$OUT/inputs.json" <<'PY'
import json, sys
out = {}
for arm, path in zip(("A", "B"), sys.argv[1:3]):
    seq = "".join(l.strip() for l in open(path) if not l.startswith(">"))
    out[arm] = {
        "path": path, "colons": seq.count(":"),
        "residues_excluding_colon": len(seq.replace(":", "")),
        "declared_chain_lengths": [len(p) for p in seq.split(":")],
    }
# The one assertion this staging step makes about itself: the arms must differ ONLY in the
# colon. If they differ in residues too, every number below is about two inputs rather
# than about the colon, and the control is worthless (§8.16, suspect the harness first).
a, b = out["A"], out["B"]
assert a["residues_excluding_colon"] == b["residues_excluding_colon"], (a, b)
assert a["colons"] == 1 and b["colons"] == 0, (a, b)
json.dump(out, open("/dev/stdout", "w"), indent=2)
PY
cat "$OUT/inputs.json"

# ---------------------------------------------------------------- predict both
nvidia-smi --query-gpu=name,memory.used,memory.total --format=csv,noheader \
  | tee "$OUT/gpu-before.txt"

for ARM in A B; do
  say "arm $ARM: colabfold_batch, offline, directly (NOT through node 03)"
  # --network none: the weights are in the image and no --host-url is passed, so a
  # prediction that needed the network fails here rather than quietly reaching
  # api.colabfold.com.
  docker run --rm --network none --gpus all \
    -v "$OUT:$OUT" -w "$OUT" "$IMAGE" \
    colabfold_batch "$OUT/$ARM/in" "$OUT/$ARM/out" \
      --model-type alphafold2_multimer_v3 \
      --msa-mode single_sequence \
      --num-models 1 \
      --num-recycle 3 \
    2>&1 | tee "$OUT/$ARM.log" || echo "arm $ARM exited non-zero: $?" | tee -a "$OUT/$ARM.log"
done

nvidia-smi --query-gpu=name,memory.used,memory.total --format=csv,noheader \
  | tee "$OUT/gpu-after.txt"

say "free space after: $(free_gb)G (was ${BEFORE_GB}G)"
python3 -I "$HERE/report_c10.py" "$OUT" | tee "$OUT/VERDICT.txt"
