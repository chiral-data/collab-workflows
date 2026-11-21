#!/usr/bin/env bash
# The offline suite: everything that can be checked without a GPU or the MSA server.
#
#   ./tests/run_fixtures.sh [workdir]            default workdir: ./.fixture-run
#   ./tests/run_fixtures.sh [workdir] <slug>...  only these cases
#
# Four sections, all driven from tests/fixtures.json so the run and the expectations cannot
# drift apart:
#
#   unit         the three checks that need only the node code and mdtraj.
#   fixtures     nodes 01 and 02 over real entries, then node 03 --dry-run, with the
#                output checked against a measurement of the input.
#   rejections   inputs the pipeline must refuse, plus the counter-cases it must accept.
#   downstream   nodes 04, 05 and 06 over ensembles built for the purpose.
#
# Node 03's sampling is out of scope here: it needs a GPU and a 1.9 GiB checkpoint. The
# dry run still exercises the staging, the parameter resolution and the alignment guard,
# which is where the plumbing defects were. tests/run_gpu_campaign.sh covers the rest.
#
# Needs docker, curl, python3 and the image from ../Dockerfile (override with
# BIOEMU2_IMAGE). Entries are fetched from RCSB once and cached. They are untrusted input:
# they stay in their own directory, nothing is executed from it, and every python3 here
# runs with -I so a file dropped beside a script cannot be imported for a stdlib module.
set -uo pipefail

TESTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WF="$(dirname "$TESTS")"
WORK="${1:-$PWD/.fixture-run}"
shift || true
ONLY=("$@")

DL="$WORK/downloads"
OUT="$WORK/results"
IMG="${BIOEMU2_IMAGE:-bioemu2:2026_11_20}"

mkdir -p "$DL" "$OUT" || exit 1

if ! docker image inspect "$IMG" >/dev/null 2>&1; then
  echo "ERROR: image '$IMG' not found. Build it from $WF/Dockerfile, or set BIOEMU2_IMAGE." >&2
  exit 1
fi

wanted () {  # no slugs given means every case
  [ ${#ONLY[@]} -eq 0 ] && return 0
  local slug
  for slug in "${ONLY[@]}"; do [ "$slug" = "$1" ] && return 0; done
  return 1
}

# One node in its own directory. --network none everywhere: nothing in this suite may reach
# the internet, so a node that quietly wants it fails here instead of succeeding for a
# reason nobody records. The exit code goes to <tag>.rc and the output to <tag>.log, both
# read by check_fixtures.py.
run_node () {  # $1 = case dir, $2 = node dir, $3 = "" | --dry-run, $4.. = env
  local D="$1" NODE="$2" FLAG="$3"; shift 3
  local TAG="${NODE%%-*}"
  local ARGS=(--network none)
  local kv
  for kv in "$@"; do ARGS+=(-e "$kv"); done
  mkdir -p "$D/$TAG/outputs"
  cp "$WF/$NODE/script.py" "$D/$TAG/"
  docker run --rm "${ARGS[@]}" -v "$D/$TAG":/workspace -w /workspace "$IMG" \
    python script.py $FLAG > "$D/$TAG.log" 2>&1
  local RC=$?
  echo "$RC" > "$D/$TAG.rc"
  return $RC
}

stage () {  # $1 = case dir, $2 = node tag, $3.. = source node tags
  local D="$1" TAG="$2"; shift 2
  local SRC
  mkdir -p "$D/$TAG/inputs"
  for SRC in "$@"; do cp -a "$D/$SRC/outputs/." "$D/$TAG/inputs/" 2>/dev/null; done
}

# The workflow at /wf and the run directory at /work, so the run directory does not have
# to live inside the workflow tree.
in_image () {
  docker run --rm --network none -v "$WF":/wf -v "$WORK":/work -w /wf "$IMG" "$@"
}

echo "== fetching entries into $DL"
for ID in $(python3 -I "$TESTS/plan.py" entries); do
  if [ ! -s "$DL/$ID.cif" ]; then
    curl -fsSL "https://files.rcsb.org/download/$ID.cif" -o "$DL/$ID.cif" \
      || { echo "ERROR: could not fetch $ID.cif" >&2; exit 1; }
    echo "   fetched $ID.cif"
  fi
done

echo
echo "== measured shape of each entry (tests/probe_entry.py)"
python3 -I "$TESTS/probe_entry.py" "$DL"

echo
echo "== unit"
mkdir -p "$OUT/unit"
in_image python -I tests/unit_checks.py /work/results/unit/unit_results.json \
  2>&1 | sed 's/^/  /'

echo
echo "== fixtures"
while IFS=$'\t' read -r SLUG ID; do
  wanted "$SLUG" || continue
  D="$OUT/fixtures/$SLUG"; rm -rf "$D"; mkdir -p "$D/01/inputs"
  cp "$DL/$ID.cif" "$D/01/inputs/"
  run_node "$D" 01-prepare-system ""; RC1=$?

  stage "$D" 02 01
  # The uploaded-alignment path, so no server is contacted. For a system with no protein
  # chain make_a3m.py writes nothing and node 02 takes its no-alignment path.
  if [ "$RC1" = "0" ]; then
    python3 -I "$TESTS/make_a3m.py" "$DL/$ID.cif" \
      "$D/01/outputs/system_inventory.json" "$D/02/inputs" >/dev/null
  fi
  run_node "$D" 02-fetch-msa ""; RC2=$?

  stage "$D" 03 02
  run_node "$D" 03-sample --dry-run "PARAM_NUM_SAMPLES=4"; RC3=$?
  printf '  %-6s node01=%s node02=%s node03(dry)=%s\n' "$ID" "$RC1" "$RC2" "$RC3"
done < <(python3 -I "$TESTS/plan.py" fixtures)

# The ensembles the downstream and node-04 cases are driven from. Built once per (source,
# frames, shape, physicality) and reused, because make_ensembles.py is deterministic.
base_for () {  # $1 = source id, $2 = n_frames, $3 = shape, $4 = physicality
  local ID="$1" N="$2" SHAPE="$3" PHYS="$4"
  local KEY="$ID-$N-$SHAPE-$PHYS"
  local B="$OUT/bases/$KEY"
  if [ ! -s "$B/samples.xtc" ]; then
    # A fixture case's slug is its entry id, so the fixture run's system.yaml is here.
    local SYS="$OUT/fixtures/$ID/01/outputs/system.yaml"
    if [ ! -s "$SYS" ]; then
      echo "   ERROR: no system.yaml for $ID; its fixture case must run first" >&2
      return 1
    fi
    mkdir -p "$B"
    cp "$SYS" "$B/system.yaml"
    in_image python -I tests/make_ensembles.py \
      "/work/results/bases/$KEY" "/work/results/bases/$KEY/system.yaml" "$N" "$SHAPE" "$PHYS" \
      > "$B/make.log" 2>&1 \
      || { echo "   ERROR: make_ensembles failed for $KEY; see $B/make.log" >&2; return 1; }
  fi
  echo "$B"
}

echo
echo "== rejections"
while IFS=$'\t' read -r SLUG ID STAGE SETUP ENVS; do
  wanted "$SLUG" || continue
  D="$OUT/rejections/$SLUG"; rm -rf "$D"; mkdir -p "$D"
  read -r -a ENV_ARR <<< "${ENVS:-}"

  if [ "$STAGE" = "node03_dryrun" ]; then
    FIX="$OUT/fixtures/$ID"
    mkdir -p "$D/03/inputs"
    cp -a "$FIX/02/outputs/." "$D/03/inputs/" 2>/dev/null
    case "$SETUP" in
      drop_msa)
        # Exactly the staging failure this suite exists for: node 02 named protein chains
        # and the msa directory never arrived.
        rm -rf "$D/03/inputs/msa" ;;
      corrupt_msa_query)
        for A3M in "$D/03/inputs/msa"/*.a3m; do
          [ -e "$A3M" ] || continue
          python3 -I -c "
import sys
path = sys.argv[1]
with open(path) as handle:
    lines = handle.read().splitlines()
# Replace the query sequence with a different one of the same length.
if len(lines) > 1:
    lines[1] = 'G' * len(lines[1])
with open(path, 'w') as handle:
    handle.write('\n'.join(lines) + '\n')
" "$A3M"
        done ;;
    esac
    run_node "$D" 03-sample --dry-run "PARAM_NUM_SAMPLES=4" "${ENV_ARR[@]}"; RC=$?
  else
    B=$(base_for "$ID" 6 spread passed) || continue
    mkdir -p "$D/04/inputs"
    cp -a "$B/." "$D/04/inputs/" 2>/dev/null
    rm -f "$D/04/inputs/make.log"
    run_node "$D" 04-propagate --dry-run "${ENV_ARR[@]}"; RC=$?
  fi
  printf '  %-62s exit=%s\n' "${SLUG:0:62}" "$RC"
done < <(python3 -I "$TESTS/plan.py" rejections)

echo
echo "== downstream"
while IFS=$'\t' read -r SLUG SOURCE N SHAPE PHYS SETUP ENVS; do
  wanted "$SLUG" || continue
  D="$OUT/downstream/$SLUG"; rm -rf "$D"; mkdir -p "$D"
  read -r -a ENV_ARR <<< "${ENVS:-}"

  B=$(base_for "$SOURCE" "$N" "$SHAPE" "$PHYS") || continue
  mkdir -p "$D/base"
  cp -a "$B/." "$D/base/" 2>/dev/null
  rm -f "$D/base/make.log"

  mkdir -p "$D/04/inputs"
  cp -a "$D/base/." "$D/04/inputs/"
  run_node "$D" 04-propagate "" "${ENV_ARR[@]}"; RC4=$?

  stage "$D" 05 04
  cp -a "$D/base/." "$D/05/inputs/"
  cp "$D/base/system.yaml" "$D/05/inputs/" 2>/dev/null
  cp "$OUT/fixtures/$SOURCE/01/outputs/system_inventory.json" "$D/05/inputs/" 2>/dev/null
  case "$SETUP" in
    drop_propagation_manifest)
      # Node 04 ran and its manifest did not arrive. Not the same thing as propagation off.
      rm -f "$D/05/inputs/propagation_manifest.json" ;;
    lose_propagated_ensemble)
      # A manifest saying propagation ran, with the ensemble it produced missing.
      rm -f "$D/05/inputs/propagated_samples.xtc"
      python3 -I -c "
import json, sys
with open(sys.argv[1], 'w') as handle:
    json.dump({'enabled': True, 'reason': 'fixture: manifest says propagation ran',
               'n_frames_in': 6, 'n_frames_out': 6}, handle, indent=2)
" "$D/05/inputs/propagation_manifest.json" ;;
  esac
  run_node "$D" 05-analyze "" "${ENV_ARR[@]}"; RC5=$?

  if [ "$RC5" = "0" ]; then
    stage "$D" 06 05
    cp "$D/base/sampling_manifest.json" "$D/06/inputs/" 2>/dev/null
    cp "$D/04/outputs/propagation_manifest.json" "$D/06/inputs/" 2>/dev/null
    run_node "$D" 06-report "" "${ENV_ARR[@]}"; RC6=$?
  else
    RC6="skipped"
  fi
  printf '  %-56s node04=%s node05=%s node06=%s\n' "${SLUG:0:56}" "$RC4" "$RC5" "$RC6"
done < <(python3 -I "$TESTS/plan.py" downstream)

echo
echo "== checking what each node reported against what the input holds"
python3 -I "$TESTS/check_fixtures.py" "$OUT" "$DL"
