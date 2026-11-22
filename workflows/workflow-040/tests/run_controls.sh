#!/bin/bash
# Tier 2, Phase 1: the controls that can falsify the scoring rather than the plumbing.
#
# C12 an idealised designed helix with a designed sequence   must PASS   (positive)
# C3  a native backbone with the sequence shuffled           must FAIL   (negative)
# C4  the same backbone with a poly-alanine sequence         must FAIL   (negative)
# C1  4ZQK PD-1, native, its own sequence on its own backbone            (diagnostic)
# C2  1BRS barnase, native, numbered from 3                              (diagnostic)
#
# C1 and C2 were the positive controls and are not any more. They were run, and they
# failed: 49.33 and 43.09 pLDDT. The control was wrong rather than the scorer -- spec §4
# predicts with msa_mode = single_sequence, which is the mode a native protein cannot be
# folded in, and C12 scoring 0.93 A at 95.05 in the same run is what tells those two
# explanations apart. They are kept as the single_sequence baseline for §4's Phase 2 MSA
# arms. See tests/reference/tier2-phase1-controls.md and spec §8.4.
#
# All five go through one run, so the separation margin between the worst positive and the
# best negative is read off a single report instead of compared across runs. If they do
# not separate on sc_rmsd_A, 040 is wrong and does not merge -- that is the whole point of
# the section, and it is the first thing in either workflow that tests the SCORE rather
# than the wiring. plddt_binder is measured but does not gate: C4 is 115 alanines and
# scores 95.81, above the real design (§8.6).
#
# This is deliberately not part of run_fixtures.sh. It needs a GPU, a network fetch from
# RCSB and several minutes of real prediction, so it is run by hand and its numbers are
# recorded. Everything it does is logged under $OUT so the numbers can be read again.
#
# Needs: docker, the colabfold:2026_10_08 image, a GPU, curl, python3.
#
# Usage:  bash tests/run_controls.sh [OUT_DIR] [DOWNLOADS_DIR]
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WF="$(dirname "$HERE")"
WORKFLOWS="$(dirname "$WF")"
OUT="${1:-$HOME/wf040-controls/$(date +%Y%m%d-%H%M%S)}"
DL="${2:-$HOME/wf040-downloads}"
IMAGE="colabfold:2026_10_08"

# The GPU host's own ~/.bashrc exports a DOCKER_HOST that does not exist, so any
# non-login shell sees no daemon at all rather than an error. Set it from the socket that
# is actually there (README, "Working on the GPU host").
if [ -z "${DOCKER_HOST:-}" ] && [ -S "/run/user/$(id -u)/docker.sock" ]; then
  export DOCKER_HOST="unix:///run/user/$(id -u)/docker.sock"
fi

say() { printf '\n=== %s\n' "$*"; }

free_gb() { df --output=avail -BG / | tail -1 | tr -dc '0-9'; }

BEFORE_GB="$(free_gb)"
say "free space before: ${BEFORE_GB}G"
if [ "$BEFORE_GB" -lt 8 ]; then
  echo "run_controls: only ${BEFORE_GB}G free. Refusing to start: this host runs at 94%" >&2
  echo "and a prediction run that fills the disk is worse than one that never ran." >&2
  exit 1
fi

docker image inspect "$IMAGE" >/dev/null 2>&1 || {
  echo "run_controls: image $IMAGE not found. Build it from workflow-026/Dockerfile." >&2
  exit 1
}

mkdir -p "$OUT" "$DL"
say "out: $OUT"

# ---------------------------------------------------------------- the entries
# Their own directory, nothing else in it, and every tool that reads them runs with -I:
# a file fetched from the PDB is untrusted input and must not be able to put a module on
# the import path.
for id in 4ZQK 1BRS; do
  if [ ! -f "$DL/$id.cif" ]; then
    say "downloading $id"
    curl -fsSL "https://files.rcsb.org/download/$id.cif" -o "$DL/$id.cif"
  fi
done
sha256sum "$DL"/4ZQK.cif "$DL"/1BRS.cif | tee "$OUT/downloads.sha256"

# Probed, then chosen from the probe. "chain A" is not a sufficient description of an
# entry: 037 lost five of 1BRS's six chains to that (53d1d460), and 4ZQK's chain B turns
# out to carry an 8-residue internal gap, which would make it a positive control that
# fails for a reason unrelated to the score.
say "probing the entries"
python3 -I "$HERE/make_controls.py" --list "$DL/4ZQK.cif" | tee "$OUT/probe-4ZQK.txt"
python3 -I "$HERE/make_controls.py" --list "$DL/1BRS.cif" | tee "$OUT/probe-1BRS.txt"

# ---------------------------------------------------------------- build the four
say "building the controls"
mkdir -p "$OUT/g"
python3 -I "$HERE/make_controls.py" "$DL/4ZQK.cif" "$OUT/g" \
  --name C1 --binder A --target B --sequence native
python3 -I "$HERE/make_controls.py" "$DL/1BRS.cif" "$OUT/g" \
  --name C2 --binder A --target D --sequence native
python3 -I "$HERE/make_controls.py" "$DL/4ZQK.cif" "$OUT/g" \
  --name C3 --binder A --target B --sequence scrambled --seed 0
python3 -I "$HERE/make_controls.py" "$DL/4ZQK.cif" "$OUT/g" \
  --name C4 --binder A --target B --sequence poly-alanine

# C12: an idealised DESIGNED sequence on an idealised backbone -- a textbook alpha helix
# (CA radius 2.3 A, rise 1.5 A, 100 deg per residue, cross-checked against 2ZTA) carrying
# make_designs.py's ALEKIAELKA cycle.
#
# Added after the first run of this script, and it is the control that tells two very
# different explanations apart. C1 and C2 are NATIVE proteins, and §4 says why this
# workflow predicts with msa_mode = single_sequence: a de novo binder has no homologs, and
# designed proteins are idealised. A native immunoglobulin domain is the opposite of that
# -- it depends on the alignment it is not being given -- so a native positive control at
# single_sequence is the HARD case, not the easy one. If C12 scores well where C1 and C2
# do not, the fault is in the choice of positive control; if C12 scores badly too, it is
# in the scorer. Neither answer can be had without it.
say "building C12: an idealised designed sequence on an idealised backbone"
mkdir -p "$OUT/g12"
python3 -I "$WORKFLOWS/workflow-037/tests/make_designs.py" "$OUT/g12" 1 helix 60
cp "$OUT/g12/design_0_generated.cif" "$OUT/g/C12_generated.cif"
# make_designs.py writes the design as chain B and the target as chain A, which is the
# same convention make_controls.py uses, so the manifest already written is correct for it.
# The provenance is MEASURED off the file just written, not asserted -- see the header of
# make_c12_provenance.py for what happened when it was asserted.
python3 -I "$HERE/make_c12_provenance.py" \
  "$OUT/g/C12_generated.cif" "$OUT/g/C12_provenance.json"

python3 -I "$HERE/probe_design.py" "$OUT/g" | tee "$OUT/probe-controls.txt"

# ---------------------------------------------------------------- run the nodes
# One function, so every node runs the same way and the only difference between the CPU
# nodes and the GPU one is the flag that says so.
run_node() {           # run_node <step-dir> <script> [extra docker args...]
  local step="$1" script="$2"; shift 2
  mkdir -p "$OUT/$step/inputs" "$OUT/$step/outputs"
  docker run --rm --network none \
    -v "$WORKFLOWS:$WORKFLOWS:ro" -v "$OUT:$OUT" -w "$OUT/$step" \
    "$@" "$IMAGE" python "$script" 2>&1 | tee "$OUT/$step.log"
}

say "037 node 04: the real design_metrics.json producer"
mkdir -p "$OUT/m/inputs"
cp "$OUT"/g/*_generated.cif "$OUT/g/generation_manifest.json" "$OUT/m/inputs/"
run_node m "$WORKFLOWS/workflow-037/04-analyze/script.py"

say "040 node 01: ingest"
mkdir -p "$OUT/n1/inputs"
cp "$OUT"/g/*_generated.cif "$OUT/m/outputs/design_metrics.json" "$OUT/n1/inputs/"
run_node n1 "$WF/01-ingest-designs/script.py"

# 01b runs here too, although this is the Phase 1 path and it has nothing to align. It is
# a dependency of node 02 in .chiral/workflow.toml, so skipping it would mean this script
# ran a DAG the workflow does not have -- and the no-op branch, the one that writes a
# target_alignment.json saying it was a no-op so that "did 01b run?" has an answer at all
# (§8.3), would never be exercised on the GPU host.
#
# The recorded Gate B Phase 1 evidence predates this step. It is unaffected: with
# score_interface off, node 02 never calls load_target, so 01b contributes nothing to any
# number in that run. See tests/reference/tier2-phase1-controls.md.
say "040 node 01b: align target (a no-op on the Phase 1 path, which is the point)"
mkdir -p "$OUT/n1b/inputs"
cp "$OUT/n1/outputs/design_inventory.json" "$OUT/n1b/inputs/"
run_node n1b "$WF/01b-align-target/script.py"

say "040 node 02: prediction inputs"
mkdir -p "$OUT/n2/inputs"
cp "$OUT"/n1/outputs/* "$OUT"/n1b/outputs/* "$OUT/n2/inputs/"
run_node n2 "$WF/02-build-prediction-inputs/script.py"

say "040 node 03: REAL prediction on the GPU, offline"
mkdir -p "$OUT/n3/inputs"
cp "$OUT"/n2/outputs/* "$OUT/n3/inputs/"
nvidia-smi --query-gpu=name,memory.used,memory.total --format=csv,noheader \
  | tee "$OUT/gpu-before.txt"
# --network none stays. The weights are baked into the image and node 03 never passes
# --host-url, so a prediction that needed the network would fail here rather than quietly
# reaching api.colabfold.com -- which is what 026 does (02_predict/script.py:58).
run_node n3 "$WF/03-predict/script.py" --gpus all
nvidia-smi --query-gpu=name,memory.used,memory.total --format=csv,noheader \
  | tee "$OUT/gpu-after.txt"

say "040 node 04: score"
mkdir -p "$OUT/n4/inputs"
cp "$OUT"/n1/outputs/* "$OUT"/n3/outputs/* "$OUT/n4/inputs/"
run_node n4 "$WF/04-score/script.py"

say "040 node 05: filter and report"
mkdir -p "$OUT/n5/inputs"
cp "$OUT"/n4/outputs/* "$OUT/n5/inputs/"
run_node n5 "$WF/05-filter-and-report/script.py"

# ---------------------------------------------------------------- the verdict
say "free space after: $(free_gb)G (was ${BEFORE_GB}G)"
python3 -I "$HERE/report_controls.py" "$OUT" | tee "$OUT/VERDICT.txt"
