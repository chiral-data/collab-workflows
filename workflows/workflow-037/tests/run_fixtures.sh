#!/usr/bin/env bash
# Run nodes 01 and 02 over every fixture in fixtures.json, then check the output
# against a measurement of the input.
#
#   ./tests/run_fixtures.sh [workdir]        default workdir: ./.fixture-run
#
# Needs: the built image (see ../Dockerfile), docker, curl, python3.
# Nodes 03-06 are not covered here: 03 needs a GPU and minutes per design, and the
# input-shape bugs this suite is for all bite in 01 and 02.
#
# Downloaded structures are untrusted data. They go in their own directory, nothing
# is executed from it, and every python3 here runs with -I so a file sitting next to
# a script cannot be imported in place of a stdlib module.
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

ids=$(python3 -I -c "
import json
with open('$TESTS/fixtures.json') as f:
    print(' '.join(x['id'] for x in json.load(f)['fixtures']))
") || exit 1

echo "== fetching structures into $DL"
for ID in $ids; do
  if [ ! -s "$DL/$ID.cif" ]; then
    curl -fsSL "https://files.rcsb.org/download/$ID.cif" -o "$DL/$ID.cif" \
      || { echo "ERROR: could not fetch $ID" >&2; exit 1; }
    echo "   fetched $ID"
  fi
done

echo
echo "== measured shape of each input (tests/probe_entry.py)"
python3 -I "$TESTS/probe_entry.py" "$DL"

echo
echo "== running nodes 01 and 02"
printf '%-8s %-7s %-7s %s\n' fixture node01 node02 warnings
printf -- '------------------------------------------------------------\n'
for ID in $ids; do
  D="$OUT/$ID"
  rm -rf "$D"
  mkdir -p "$D/n1/inputs" "$D/n1/outputs" "$D/n2/inputs" "$D/n2/outputs"
  cp "$WF/01-acquire-target/script.py" "$WF/01-acquire-target/run.sh" "$D/n1/"
  cp "$WF/02-build-spec/script.py" "$WF/02-build-spec/run.sh" "$D/n2/"
  cp "$DL/$ID.cif" "$D/n1/inputs/"

  mapfile -t ENVS < <(python3 -I -c "
import json
with open('$TESTS/fixtures.json') as f:
    for x in json.load(f)['fixtures']:
        if x['id'] == '$ID':
            for k, v in x['params'].items():
                print(f'{k}={v}')
")
  ARGS=()
  for kv in "${ENVS[@]}"; do ARGS+=(-e "$kv"); done

  docker run --rm --network none "${ARGS[@]}" \
    -v "$D/n1":/workspace -w /workspace "$IMG" bash run.sh > "$D/n1.log" 2>&1
  RC1=$?
  cp -a "$D/n1/outputs/." "$D/n2/inputs/" 2>/dev/null
  docker run --rm --network none "${ARGS[@]}" \
    -v "$D/n2":/workspace -w /workspace "$IMG" bash run.sh > "$D/n2.log" 2>&1
  RC2=$?

  W1=$(grep -c '^WARNING' "$D/n1.log" 2>/dev/null)
  W2=$(grep -c '^WARNING' "$D/n2.log" 2>/dev/null)
  printf '%-8s %-7s %-7s %s\n' "$ID" "$RC1" "$RC2" "${W1:-0}/${W2:-0}"
done

echo
echo "== checking output against measured input"
python3 -I "$TESTS/check_fixtures.py" "$OUT" "$DL"
