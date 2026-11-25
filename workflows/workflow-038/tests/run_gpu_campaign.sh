#!/usr/bin/env bash
# Run the whole pipeline, for real, over several literature systems on a GPU host.
#
#   ./tests/run_gpu_campaign.sh [workdir] [system ...]      default workdir: ./.gpu-campaign
#
# This is the counterpart to run_fixtures.sh, which is offline and covers the plumbing.
# Everything here needs what that one deliberately does not: a GPU, the public ColabFold
# server, and minutes per system. It is how the claims that cannot be checked without a
# model get checked -- whether the ensemble is physically sensible, whether a heteromer's
# paired alignment survives the round trip, whether a run fits in the card's memory.
#
# Needs: the built image (see ../Dockerfile), docker with --gpus, curl, python3.
# Structures are fetched from RCSB once and cached. They are untrusted input: they stay in
# their own directory, nothing is executed from it, and the measuring script runs with -I.
set -uo pipefail

TESTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WF="$(dirname "$TESTS")"
WORK="${1:-$PWD/.gpu-campaign}"
shift || true
SYSTEMS=("$@")
if [ ${#SYSTEMS[@]} -eq 0 ]; then
  # 1UBQ   control: a small, rigid, exhaustively characterised fold. Rg is known.
  # 1CLL   four structural Ca2+ ions, every one an atom named CA. The selection regression.
  # 5PTI   three disulfides, so bad_bonds_sulfur and the system's declared bonds are live.
  # 2ZTA   homodimer: two chains, one sequence, so no paired alignment is due.
  # 1A3N   heteromer: four chains over two sequences, so a paired alignment IS due.
  # 1AKE   upstream's own example and a real two-state enzyme: the one case where the
  #        ensemble is expected to be interesting rather than merely well formed.
  # 1BNA   no protein chain at all, so no alignment is due and none must be demanded.
  SYSTEMS=(1UBQ 1CLL 5PTI 2ZTA 1A3N 1AKE 1BNA)
fi

DL="$WORK/downloads"
OUT="$WORK/results"
IMG="${BIOEMU2_IMAGE:-bioemu2:2026_11_20}"
SAMPLES="${NUM_SAMPLES:-8}"
# 0 is BioEmu2's own default: denoise every sample in one batch. Overridable because that is
# the knob a system too large for the card needs -- measured on an 8 GB RTX 3060 Ti, 542 tokens
# runs out of memory at batch 8 while upstream's own estimate said 33 would fit.
BATCH="${BATCH_SIZE:-0}"

mkdir -p "$DL" "$OUT" || exit 1

if ! docker image inspect "$IMG" >/dev/null 2>&1; then
  echo "ERROR: image '$IMG' not found. Build it from $WF/Dockerfile, or set BIOEMU2_IMAGE." >&2
  exit 1
fi
if ! docker run --rm --gpus all "$IMG" true >/dev/null 2>&1; then
  echo "ERROR: 'docker run --gpus all' does not work here. This script needs a GPU." >&2
  exit 1
fi

echo "== image and device"
docker run --rm --gpus all "$IMG" python -c "
import torch
print('  torch', torch.__version__, '| CUDA', torch.version.cuda)
print('  architectures compiled in:', ' '.join(torch.cuda.get_arch_list()))
" 2>/dev/null
docker run --rm --gpus all "$IMG" \
  nvidia-smi --query-gpu=name,compute_cap,driver_version,memory.total --format=csv,noheader \
  2>/dev/null | sed 's/^/  /'

echo
echo "== fetching structures into $DL"
for ID in "${SYSTEMS[@]}"; do
  if [ ! -s "$DL/$ID.cif" ]; then
    curl -fsSL "https://files.rcsb.org/download/$ID.cif" -o "$DL/$ID.cif" \
      || { echo "ERROR: could not fetch $ID.cif" >&2; exit 1; }
    echo "   fetched $ID.cif"
  fi
done

echo
echo "== measured shape of each input (tests/probe_entry.py)"
python3 -I "$TESTS/probe_entry.py" "$DL"

# One node, in its own directory, with its inputs copied from the nodes it depends on.
# --network none everywhere except node 02, which is the only node that may reach the
# internet -- so a node that quietly wants the network fails here rather than succeeding
# for a reason nobody records.
stage () {  # $1 = case dir, $2 = node dir name, $3.. = source node dirs
  local D="$1" NODE="$2"; shift 2
  local TAG="${NODE%%-*}"
  mkdir -p "$D/$TAG/inputs" "$D/$TAG/outputs"
  cp "$WF/$NODE/script.py" "$WF/$NODE/run.sh" "$D/$TAG/"
  for SRC in "$@"; do
    cp -a "$D/$SRC/outputs/." "$D/$TAG/inputs/" 2>/dev/null
  done
}

run_node () {  # $1 = case dir, $2 = node dir name, $3 = net|none, $4 = gpu|cpu, $5.. = env
  local D="$1" NODE="$2" NET="$3" GPU="$4"; shift 4
  local TAG="${NODE%%-*}"
  local ARGS=()
  [ "$NET" = "none" ] && ARGS+=(--network none)
  [ "$GPU" = "gpu" ] && ARGS+=(--gpus all)
  # bioemu2 caches each chain's ColabFold result under BIOEMU2_MSA_CACHE, which the image
  # points inside the container. Pointing it at a host directory makes a second campaign skip
  # the server entirely -- which matters because the MSA fetches, not the GPU, are most of the
  # wall time, and because the public server rate-limits.
  if [ "$NET" != "none" ]; then
    mkdir -p "$WORK/msa-cache"
    ARGS+=(-v "$WORK/msa-cache":/msa-cache -e BIOEMU2_MSA_CACHE=/msa-cache)
  fi
  for kv in "$@"; do ARGS+=(-e "$kv"); done
  docker run --rm "${ARGS[@]}" -v "$D/$TAG":/workspace -w /workspace "$IMG" \
    bash run.sh > "$D/$TAG.log" 2>&1
  local RC=$?
  # The status goes next to the log so check_campaign.py can tell a node that ran and
  # produced nothing from a node that never ran -- and can require that a node which
  # stopped said why before it did.
  echo "$RC" > "$D/$TAG.rc"
  return $RC
}

echo
echo "== full pipeline, $SAMPLES sample(s) per system"
printf '%-8s %-5s %-5s %-5s %-5s %-5s %-5s  %s\n' system n01 n02 n03 n04 n05 n06 notes
printf -- '--------------------------------------------------------------------------\n'

for ID in "${SYSTEMS[@]}"; do
  D="$OUT/$ID"; rm -rf "$D"; mkdir -p "$D"

  stage "$D" 01-prepare-system
  cp "$DL/$ID.cif" "$D/01/inputs/"
  run_node "$D" 01-prepare-system none cpu; RC1=$?

  stage "$D" 02-fetch-msa 01
  run_node "$D" 02-fetch-msa net cpu; RC2=$?

  stage "$D" 03-sample 02
  run_node "$D" 03-sample none gpu \
    "PARAM_MODE=equilibrium" "PARAM_NUM_SAMPLES=$SAMPLES" "PARAM_SEED=42" \
    "PARAM_BATCH_SIZE=$BATCH"; RC3=$?

  stage "$D" 04-propagate 03
  run_node "$D" 04-propagate none gpu; RC4=$?

  stage "$D" 05-analyze 03 04
  run_node "$D" 05-analyze none cpu; RC5=$?

  stage "$D" 06-report 03 04 05
  run_node "$D" 06-report none cpu; RC6=$?

  NOTE=$(python3 -I -c "
import json, os
path = '$D/05/outputs/ensemble_summary.json'
if not os.path.exists(path):
    print('no summary'); raise SystemExit
d = json.load(open(path))
s = d['ensembles']['sampled']
pw = s.get('pairwise_rmsd') or {}
rg = s.get('radius_of_gyration_nm') or {}
print(f\"Rg {rg.get('mean')} nm, max RMSD {pw.get('max_nm')} nm, \"
      f\"{s.get('n_clusters')} cluster(s), CA set {len(s.get('residue_rmsf') or [])}\")
" 2>/dev/null)
  printf '%-8s %-5s %-5s %-5s %-5s %-5s %-5s  %s\n' \
    "$ID" "$RC1" "$RC2" "$RC3" "$RC4" "$RC5" "$RC6" "$NOTE"
done

# The time-lagged path, run once for real. Nothing else in either suite executes
# propagate_with_embeddings: the offline suite only reaches node 04's identity pass-through,
# and the loop above leaves `propagate` at its default of false. A feature the workflow
# advertises and no test has ever run is a feature nobody knows works.
PROP_SRC="1UBQ"
for ID in "${SYSTEMS[@]}"; do [ "$ID" = "$PROP_SRC" ] && FOUND=1; done
[ "${FOUND:-}" = "1" ] || PROP_SRC="${SYSTEMS[0]}"

if [ -s "$OUT/$PROP_SRC/03/outputs/samples.xtc" ]; then
  echo
  echo "== time-lagged propagation, for real, on $PROP_SRC"
  D="$OUT/$PROP_SRC-propagated"; rm -rf "$D"; mkdir -p "$D"
  # Node 04 takes node 03's outputs, exactly as the DAG gives them to it.
  mkdir -p "$D/03/outputs"
  cp -a "$OUT/$PROP_SRC/03/outputs/." "$D/03/outputs/"
  cp -a "$OUT/$PROP_SRC/01/outputs/system_inventory.json" "$D/03/outputs/" 2>/dev/null

  stage "$D" 04-propagate 03
  run_node "$D" 04-propagate none gpu \
    "PARAM_PROPAGATE=true" "PARAM_LAG_TIME_NS=1.0" "PARAM_SEED=42"; RC4=$?

  stage "$D" 05-analyze 03 04
  run_node "$D" 05-analyze none cpu; RC5=$?

  stage "$D" 06-report 03 04 05
  run_node "$D" 06-report none cpu; RC6=$?

  printf '%-8s %-5s %-5s %-5s %-5s %-5s %-5s  %s\n' \
    "$PROP_SRC+p" "-" "-" "-" "$RC4" "$RC5" "$RC6" \
    "$(python3 -I -c "
import json, os
path = '$D/05/outputs/ensemble_summary.json'
if not os.path.exists(path):
    print('no summary'); raise SystemExit
d = json.load(open(path))
labels = sorted((d.get('ensembles') or {}).keys())
out = []
for label in labels:
    e = d['ensembles'][label]
    rg = (e.get('radius_of_gyration_nm') or {}).get('mean')
    pw = (e.get('pairwise_rmsd') or {}).get('max_nm')
    out.append(f\"{label}: Rg {rg} nm, max RMSD {pw} nm\")
print(' | '.join(out))
" 2>/dev/null)"
fi

# The cases this run actually drove, recorded for the checker to read rather than infer.
# check_campaign.py used to take every directory under $OUT as a case, which made a scratch
# directory left in the tree look like a failed system, and -- worse -- made a system this
# driver was asked for but never created a directory for simply absent from the report.
: > "$OUT/campaign_cases.txt"
for ID in "${SYSTEMS[@]}"; do echo "$ID" >> "$OUT/campaign_cases.txt"; done
[ -d "$OUT/$PROP_SRC-propagated" ] && echo "$PROP_SRC-propagated" >> "$OUT/campaign_cases.txt"

echo
echo "== checking the measured ensembles against the measured inputs"
python3 -I "$TESTS/check_campaign.py" "$OUT" "$DL"
