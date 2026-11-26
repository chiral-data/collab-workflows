#!/bin/bash
# C6: a REAL prediction killed at a short budget. The last Phase 1 control.
#
# What is being tested is node 03's kill_tree(): Popen in its own process group, SIGTERM
# then SIGKILL, because subprocess.run's timeout signals the direct child only while the
# predictor's workers go on holding GPU memory. Workflow-037 hit exactly that -- its
# generation spawned workers that survived the parent and had to be freed by hand
# (755d4af2) -- and 040 copied the fix. This is the control that shows the fix working
# against a real colabfold_batch rather than against a shell stand-in.
#
# WHY THIS DOES NOT JUST USE run_controls.sh's docker run
#
# Tier 1 covers the timeout path with a stand-in predictor, and the suite reports a GAP
# for the worker-survival claim in the image, because a pid recorded inside a container
# is meaningless on the host. The same problem wrecks the naive version of this control:
# `docker run` a node that times out and the container exits, and a container's exit kills
# every process in it and frees the GPU no matter how badly kill_tree() is written. The
# check would pass for the wrong reason, which is worse than no check (§8.4).
#
# So the container is started FIRST and outlives the node:
#
#   1. docker run -d ... sleep infinity          (the container survives the node)
#   2. docker exec  -> node 03, short budget     (the node times out and kills its group)
#   3. docker exec  -> ps inside that container  (did anything survive? answerable now)
#   4. nvidia-smi on the host                    (was the device actually freed?)
#
# And step 3 is only worth anything if workers existed in the first place, so a sampler
# runs alongside and records what was alive mid-prediction. "Nothing survived" is vacuous
# if nothing was ever there.
#
# Needs: docker, the colabfold:2026_10_08 image, a GPU, and an existing control run to
# take node 02's outputs from (its real mono prediction inputs).
#
# Usage:  bash tests/run_c6.sh <control-run-dir> [OUT_DIR] [TIMEOUT_S]
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WF="$(dirname "$HERE")"
WORKFLOWS="$(dirname "$WF")"
SRC="${1:?usage: run_c6.sh <control-run-dir> [OUT_DIR] [TIMEOUT_S]}"
OUT="${2:-$HOME/wf040-c6/$(date +%Y%m%d-%H%M%S)}"
# 35s, not 20s. At 20s the kill lands before colabfold_batch has printed anything, so the
# recorded reason carries an empty "Last lines:" -- true, and useless to whoever reads it.
# 35s is still well inside the 45-67s a real prediction of these inputs takes, so the kill
# is certain, and by then the predictor has said something worth retaining.
BUDGET="${3:-35}"
IMAGE="colabfold:2026_10_08"
NAME="wf040-c6-$$"

if [ -z "${DOCKER_HOST:-}" ] && [ -S "/run/user/$(id -u)/docker.sock" ]; then
  export DOCKER_HOST="unix:///run/user/$(id -u)/docker.sock"
fi

say() { printf '\n=== %s\n' "$*"; }
free_gb() { df --output=avail -BG / | tail -1 | tr -dc '0-9'; }

BEFORE_GB="$(free_gb)"
say "free space before: ${BEFORE_GB}G"
[ "$BEFORE_GB" -lt 8 ] && { echo "run_c6: only ${BEFORE_GB}G free, refusing" >&2; exit 1; }

[ -d "$SRC/n2/outputs" ] || { echo "run_c6: no $SRC/n2/outputs" >&2; exit 2; }
docker image inspect "$IMAGE" >/dev/null 2>&1 || {
  echo "run_c6: image $IMAGE not found" >&2; exit 1; }

mkdir -p "$OUT/n3/inputs" "$OUT/n3/outputs"
cp "$SRC"/n2/outputs/* "$OUT/n3/inputs/"
say "out: $OUT  (budget ${BUDGET}s, inputs from $SRC)"
ls "$OUT/n3/inputs" | sed 's/^/  /'

cleanup() {
  docker rm -f "$NAME" >/dev/null 2>&1
}
trap cleanup EXIT

# ------------------------------------------------------------------ 1. the container
# --network none stays: the weights are in the image and node 03 never passes --host-url,
# so a prediction that needed the network would fail here rather than quietly reaching
# api.colabfold.com.
say "starting a container that will OUTLIVE the node"
docker run -d --rm --name "$NAME" --gpus all --network none \
  -v "$WORKFLOWS:$WORKFLOWS:ro" -v "$OUT:$OUT" -w "$OUT/n3" \
  "$IMAGE" sleep infinity > "$OUT/container-id.txt"
echo "  container $NAME: $(cut -c1-12 < "$OUT/container-id.txt")"

nvidia-smi --query-gpu=name,memory.used,memory.total --format=csv,noheader \
  | tee "$OUT/gpu-before.txt" | sed 's/^/  /'

# ------------------------------------------------------------------ 2. the sampler
# Records what is alive inside the container while the prediction runs, so that step 3's
# "nothing survived" is a measurement and not a tautology.
# `docker top`, not `docker exec ps`. The image has no ps -- the first version of this
# script used one and the OCI "executable file not found" error landed in the after-listing,
# where the judge read a file with no process lines in it as "nothing survived" and passed.
# A check that passes for the wrong reason is worse than no check (§8.4), and this one was
# in the checking code rather than in the workflow.
#
# docker top runs the HOST's ps against the container's pids, so it needs nothing inside
# the image, and it prints a header row the judge can require as proof it really listed.
procs() {
  docker top "$NAME" -eo pid,ppid,pgid,etime,rss,args 2>&1
}

sample() {
  local n=0
  while [ "$n" -lt 400 ]; do
    {
      echo "--- t+${n}s"
      procs
      nvidia-smi --query-gpu=memory.used --format=csv,noheader 2>/dev/null \
        | sed 's/^/gpu: /'
    } >> "$OUT/sampler.txt"
    sleep 2
    n=$((n + 2))
  done
}
sample &
SAMPLER=$!

# ------------------------------------------------------------------ 3. the node
say "node 03 inside that container, PARAM_PREDICT_TIMEOUT_S=$BUDGET"
START="$(date +%s)"
docker exec \
  -e PARAM_PREDICT_TIMEOUT_S="$BUDGET" \
  -e PARAM_MSA_MODE=single_sequence \
  -e PARAM_NUM_MODELS=5 \
  -e PARAM_NUM_RECYCLE=3 \
  -w "$OUT/n3" "$NAME" \
  python "$WF/03-predict/script.py" > "$OUT/n3.log" 2>&1
NODE_RC=$?
ELAPSED=$(( $(date +%s) - START ))
echo "  node 03 exited $NODE_RC after ${ELAPSED}s"

kill "$SAMPLER" 2>/dev/null
wait "$SAMPLER" 2>/dev/null

# ------------------------------------------------------------------ 4. what survived
# The whole point of the exercise, and only askable because the container is still up.
say "what is still alive INSIDE the container, after the node returned"
sleep 3
procs | tee "$OUT/ps-after.txt" | sed 's/^/  /'

nvidia-smi --query-gpu=name,memory.used,memory.total --format=csv,noheader \
  | tee "$OUT/gpu-after.txt" | sed 's/^/  /'

printf '%s\n' "$NODE_RC" > "$OUT/node-rc.txt"
printf '%s\n' "$BUDGET" > "$OUT/budget.txt"
printf '%s\n' "$ELAPSED" > "$OUT/elapsed.txt"

say "free space after: $(free_gb)G (was ${BEFORE_GB}G)"
say "verdict"
python3 -I "$HERE/report_c6.py" "$OUT" | tee "$OUT/VERDICT.txt"
