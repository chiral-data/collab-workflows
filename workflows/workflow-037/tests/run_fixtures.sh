#!/usr/bin/env bash
# Run the fixture suite: nodes 01 and 02 over real structures, node 02's refusals, and
# node 04 over structures built for the purpose. Then check it all against the inputs.
#
#   ./tests/run_fixtures.sh [workdir]        default workdir: ./.fixture-run
#
# Needs: the built image (see ../Dockerfile), docker, curl, python3.
# Node 03 is not covered: it needs a GPU and minutes per design. Node 04 is reached
# anyway, with inputs from tests/make_designs.py -- see that file for why that is sound.
#
# Downloaded structures are untrusted data. They go in their own directory, nothing is
# executed from it, and every python3 here runs with -I so a file dropped next to a
# script cannot be imported in place of a stdlib module.
set -uo pipefail

TESTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WF="$(dirname "$TESTS")"
WORK="${1:-$PWD/.fixture-run}"
DL="$WORK/downloads"
OUT="$WORK/results"
IMG="${RFD4_IMAGE:-rfd4-proteina:2026_11_04}"

mkdir -p "$DL" "$OUT" || exit 1

if ! docker image inspect "$IMG" >/dev/null 2>&1; then
  echo "ERROR: image '$IMG' not found. Build it from $WF/Dockerfile, or set RFD4_IMAGE." >&2
  exit 1
fi

# "<fixture id> <source filename>" per line; the id is not always the filename, because
# one entry appears twice -- once as mmCIF and once as PDB.
plan=$(python3 -I -c "
import json
with open('$TESTS/fixtures.json') as f:
    for x in json.load(f)['fixtures']:
        print(x['id'], x.get('source', x['id'] + '.cif'))
") || exit 1

echo "== fetching structures into $DL"
while read -r ID SRC; do
  [ -n "$ID" ] || continue
  if [ ! -s "$DL/$SRC" ]; then
    curl -fsSL "https://files.rcsb.org/download/$SRC" -o "$DL/$SRC" \
      || { echo "ERROR: could not fetch $SRC" >&2; exit 1; }
    echo "   fetched $SRC"
  fi
done <<< "$plan"

echo
echo "== measured shape of each mmCIF input (tests/probe_entry.py)"
python3 -I "$TESTS/probe_entry.py" "$DL"

params_for () {  # $1 = section, $2 = key
  python3 -I -c "
import json
with open('$TESTS/fixtures.json') as f:
    doc = json.load(f)
for x in doc.get('$1', []):
    if x.get('id', x.get('name')) == '''$2''' or x.get('name') == '''$2''':
        for k, v in x.get('params', {}).items():
            print(f'{k}={v}')
        break
"
}

run_nodes () {  # $1 = dir, $2 = source path, $3.. = env
  local D="$1" SRCPATH="$2"; shift 2
  mkdir -p "$D/n1/inputs" "$D/n1/outputs" "$D/n2/inputs" "$D/n2/outputs"
  cp "$WF/01-acquire-target/script.py" "$WF/01-acquire-target/run.sh" "$D/n1/"
  cp "$WF/02-build-spec/script.py" "$WF/02-build-spec/run.sh" "$D/n2/"
  cp "$SRCPATH" "$D/n1/inputs/"
  local ARGS=()
  for kv in "$@"; do ARGS+=(-e "$kv"); done
  docker run --rm --network none "${ARGS[@]}" \
    -v "$D/n1":/workspace -w /workspace "$IMG" bash run.sh > "$D/n1.log" 2>&1
  RC1=$?
  cp -a "$D/n1/outputs/." "$D/n2/inputs/" 2>/dev/null
  docker run --rm --network none "${ARGS[@]}" \
    -v "$D/n2":/workspace -w /workspace "$IMG" bash run.sh > "$D/n2.log" 2>&1
  RC2=$?
}

echo
echo "== nodes 01 and 02 over every fixture"
printf '%-10s %-7s %-7s %s\n' fixture node01 node02 warnings
printf -- '----------------------------------------------------------------\n'
while read -r ID SRC; do
  [ -n "$ID" ] || continue
  D="$OUT/$ID"; rm -rf "$D"
  mapfile -t ENVS < <(params_for fixtures "$ID")
  run_nodes "$D" "$DL/$SRC" "${ENVS[@]}"
  W1=$(grep -c '^WARNING' "$D/n1.log" 2>/dev/null)
  W2=$(grep -c '^WARNING' "$D/n2.log" 2>/dev/null)
  printf '%-10s %-7s %-7s %s\n' "$ID" "$RC1" "$RC2" "${W1:-0}/${W2:-0}"
done <<< "$plan"

echo
echo "== inputs node 02 must refuse"
names=$(python3 -I -c "
import json
with open('$TESTS/fixtures.json') as f:
    for x in json.load(f).get('rejections', []):
        print(x['name'], x['id'], sep='|')
")
while IFS='|' read -r NAME ID; do
  [ -n "$NAME" ] || continue
  SLUG="${NAME// /-}"
  D="$OUT/_rejections/$SLUG"; rm -rf "$D"
  mapfile -t ENVS < <(params_for rejections "$NAME")
  SRC=$(python3 -I -c "
import json
with open('$TESTS/fixtures.json') as f:
    for x in json.load(f)['fixtures']:
        if x['id'] == '$ID':
            print(x.get('source', x['id'] + '.cif')); break
")
  run_nodes "$D/run" "$DL/$SRC" "${ENVS[@]}"
  mkdir -p "$D"
  echo "$RC2" > "$D/rc"
  cp "$D/run/n2.log" "$D/log" 2>/dev/null
  printf '  %-42s node02 rc=%s\n' "$NAME" "$RC2"
done <<< "$names"

echo
echo "== node 04, on structures from tests/make_designs.py"
cases=$(python3 -I -c "
import json
with open('$TESTS/fixtures.json') as f:
    for x in json.load(f).get('node04', []):
        print(x['name'], x['n_designs'], str(x.get('swap_chains', False)), sep='|')
")
while IFS='|' read -r NAME N SWAP; do
  [ -n "$NAME" ] || continue
  D="$OUT/_node04/$NAME"; rm -rf "$D"; mkdir -p "$D/inputs" "$D/outputs"
  cp "$WF/04-analyze/script.py" "$WF/04-analyze/run.sh" "$D/"
  if [ "$SWAP" = "True" ]; then
    python3 -I "$TESTS/make_designs.py" "$D/inputs" "$N" swap > /dev/null
  else
    python3 -I "$TESTS/make_designs.py" "$D/inputs" "$N" > /dev/null
  fi
  docker run --rm --network none -v "$D":/workspace -w /workspace \
    "$IMG" bash run.sh > "$D/log" 2>&1
  printf '  %-42s node04 rc=%s  (%s design(s))\n' "$NAME" "$?" "$N"
done <<< "$cases"

echo
echo "== checking output against measured input"
python3 -I "$TESTS/check_fixtures.py" "$OUT" "$DL"
