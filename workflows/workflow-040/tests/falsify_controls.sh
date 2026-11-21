#!/bin/bash
# Show that Gate B can refuse, not just agree.
#
# report_controls.py printed "Gate B Phase 1: MET" for the control run. That is worth
# exactly nothing until the same judge, on the same run, can be made to say NOT MET --
# §8.13, a guard that cannot be shown refusing is not a guard, and §8.4, a check that
# passes for the wrong reason is worse than no check.
#
# Each mutation below changes one number or flag in a COPY of the run's design_scores.json
# and requires a non-zero exit and a NOT MET line. The run itself is never touched.
#
# Usage:  bash tests/falsify_controls.sh <control-run-dir>
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
RUN="${1:?usage: falsify_controls.sh <control-run-dir>}"
WORK="${TMPDIR:-/tmp}/wf040-gateb-falsify.$$"
SCORES="n4/outputs/design_scores.json"

[ -f "$RUN/$SCORES" ] || { echo "no $SCORES under $RUN" >&2; exit 2; }

# The honest run must pass first. If it does not, the mutations below prove nothing,
# because every one of them would "fail" for the pre-existing reason (§8.16).
if ! python3 -I "$HERE/report_controls.py" "$RUN" >/dev/null 2>&1; then
  echo "FAIL: $RUN does not pass Gate B to begin with, so falsifying it is meaningless"
  exit 1
fi
echo "baseline: $RUN passes Gate B"
echo

fails=0
run_mutation() {        # run_mutation <name> <python-expression-file>
  local name="$1" script="$2"
  rm -rf "$WORK"; mkdir -p "$WORK"
  cp -r "$RUN/." "$WORK/"
  python3 -I "$script" "$WORK/$SCORES" "$name" || { echo "  $name: mutation failed"; return 1; }
  local out rc
  out="$(python3 -I "$HERE/report_controls.py" "$WORK" 2>&1)"; rc=$?
  local verdict
  verdict="$(printf '%s\n' "$out" | grep 'Gate B' | head -2 | tr '\n' ';')"
  if [ "$rc" -ne 0 ] && printf '%s\n' "$out" | grep -q 'NOT MET'; then
    echo "  REFUSED   $name"
    echo "            $verdict"
  else
    echo "  NOT CAUGHT  $name  (exit $rc) -- the gate accepted a corrupted run"
    fails=$((fails + 1))
  fi
}

MUT="$WORK.mutate.py"
mkdir -p "$(dirname "$MUT")"
cat > "$MUT" <<'PY'
import json
import sys

path, mode = sys.argv[1], sys.argv[2]
with open(path) as handle:
    doc = json.load(handle)
hit = 0
for row in doc["designs"]:
    name = row["design"]
    if mode == "positive-fails" and name == "C12":
        row["passes"], row["failed_thresholds"] = False, "forced by falsify_controls"
        hit += 1
    elif mode == "negative-passes" and name == "C4":
        row["passes"], row["failed_thresholds"] = True, ""
        hit += 1
    elif mode == "negative-rmsd-as-good-as-positive" and name == "C3":
        row["sc_rmsd_A"] = 0.5
        hit += 1
    elif mode == "positive-rmsd-worse-than-negative" and name == "C12":
        row["sc_rmsd_A"] = 99.0
        hit += 1
    elif mode == "positive-metric-missing" and name == "C12":
        row["sc_rmsd_A"] = None
        hit += 1
if not hit:
    print(f"mutate: {mode} matched no design", file=sys.stderr)
    raise SystemExit(1)
with open(path, "w") as handle:
    json.dump(doc, handle, indent=2)
PY

echo "mutations (each must be REFUSED):"
run_mutation positive-fails "$MUT"
run_mutation negative-passes "$MUT"
run_mutation negative-rmsd-as-good-as-positive "$MUT"
run_mutation positive-rmsd-worse-than-negative "$MUT"
run_mutation positive-metric-missing "$MUT"

# And the one that was found by accident rather than by design: C12's provenance was a
# four-key hand-written dict, the printer raised a KeyError on the residue fields, and the
# crash landed after the table and before the verdict. The gate now refuses a gating
# control with no provenance instead of dying halfway through the report.
rm -rf "$WORK"; mkdir -p "$WORK"; cp -r "$RUN/." "$WORK/"
rm -f "$WORK/g/C12_provenance.json"
out="$(python3 -I "$HERE/report_controls.py" "$WORK" 2>&1)"; rc=$?
if [ "$rc" -ne 0 ] && printf '%s\n' "$out" | grep -q 'NOT MET'; then
  echo "  REFUSED   positive-provenance-removed"
  echo "            $(printf '%s\n' "$out" | grep 'Gate B' | head -1)"
else
  echo "  NOT CAUGHT  positive-provenance-removed (exit $rc)"
  fails=$((fails + 1))
fi

rm -rf "$WORK" "$MUT"
echo
if [ "$fails" -eq 0 ]; then
  echo "Gate B is falsifiable: every mutation was refused."
  exit 0
fi
echo "Gate B accepted $fails corrupted run(s). The MET verdict does not mean anything."
exit 1
