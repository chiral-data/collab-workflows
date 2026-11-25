#!/bin/bash
# Attack report_controls_phase2.py. Every mutation below must be REFUSED.
#
# Why this exists
# ---------------
# A MET verdict from a judge that has never been shown refusing is not evidence. The Phase
# 1 version of this script is the proof: it found two real defects in a judge that had
# already printed MET, and neither was visible by reading the code. One was a crash that
# landed after the table and before the verdict, so the report read as complete with no
# verdict in it; the other was a deleted provenance file passing while a partial one
# refused, because the judge only iterated records that existed.
#
# Each mutation is applied to a COPY of the run, the judge is re-run over the copy, and a
# zero exit from the judge means the mutation got through. Value-level mutations come
# before deletions: a deleted file produces a shorter, louder refusal that would mask a
# subtler one that should also have been caught.
#
# --synthetic, and why it had to exist
# ------------------------------------
# The honest-run check below aborts on every run this repo has: Gate B Phase 2 is NOT MET,
# because the offline target alignment floors the target to ~38 pLDDT and takes iptm with
# it (tests/reference/tier2-phase2-controls.md). That abort is correct -- grading mutations
# against a tree that already fails proves nothing -- but it also means the judge could
# never be shown refusing at all, and 8.13 says a guard that cannot be shown refusing is
# not a guard. So --synthetic fabricates the minimal tree the judge accepts, with a
# separation that was typed rather than measured, purely so the mutations can be graded.
#
# The fabricated numbers are evidence about the JUDGE and about nothing else. They say
# nothing about this workflow's scoring, they are not a Phase 2 result, and the tree itself
# carries a _FABRICATED key saying so -- because a file of plausible scores sitting in a
# temp directory is exactly the thing that gets mistaken for a measurement later.
#
# Usage:  bash tests/falsify_controls_phase2.sh <RUN_DIR>
#         bash tests/falsify_controls_phase2.sh --synthetic
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
JUDGE="$HERE/report_controls_phase2.py"
WORK="${TMPDIR:-/tmp}/wf040-falsify-p2.$$"
SCORES="n4/outputs/design_scores.json"

SYNTHETIC=0
if [ "${1:-}" = "--synthetic" ]; then SYNTHETIC=1; shift; fi

if [ "$SYNTHETIC" = "1" ]; then
  RUN="$WORK/synthetic-run"
  mkdir -p "$RUN/n4/outputs"
  echo "=== building a FABRICATED passing run: this grades the judge, not the workflow"
  python3 -I - "$RUN/$SCORES" <<'FAKE' || exit 2
import json, sys
# Round, obviously-typed numbers. A fabricated tree should not be able to pass for a
# measured one, so nothing here is near a real figure from any recorded run.
json.dump({
    "_FABRICATED": "Typed by falsify_controls_phase2.sh --synthetic to exercise the "
                   "judge. NOT a measurement. Not a Phase 2 result. See the script "
                   "header and tests/reference/tier2-phase2-controls.md.",
    "designs": [
        {"design": "C8", "iptm": 0.90, "pae_interaction": 5.0,
         "plddt_target": 85.0, "plddt_binder_in_complex": 80.0},
        # Both chains above the judge's 70.0 modelled floor for BOTH controls, or the
        # run is NOT INTERPRETABLE and the honest-run check below fails. A scramble can
        # be modelled confidently -- Phase 1 measured 115 alanines at 95.81 -- so a
        # well-modelled negative with a bad interface is the right shape here.
        {"design": "C9", "iptm": 0.10, "pae_interaction": 25.0,
         "plddt_target": 84.0, "plddt_binder_in_complex": 75.0},
    ],
}, open(sys.argv[1], "w"), indent=2)
FAKE
  {
    echo "PARAM_MSA_MODE_CPLX=single_sequence refused"
    echo "exit=1 (non-zero is the pass)"
  } > "$RUN/n3-refuse.log"
else
  RUN="${1:?usage: falsify_controls_phase2.sh <RUN_DIR> | --synthetic}"
fi

[ -f "$RUN/$SCORES" ] || { echo "no $SCORES under $RUN" >&2; exit 2; }

# WHICH CONTROL IS WHICH, read from the run rather than hardcoded. run_controls_phase2.sh
# writes pair.txt with the positive on line 5 and the negative on line 6, and the judge
# reads the same file for the same reason: the run says what it is, so this script cannot
# be pointed at a run and told it is something else. A synthetic tree has no pair.txt and
# is the C8 pair, which is what it was fabricated as.
#
# This is not cosmetic. Before it, every mutation body below named C8 and C9 literally, so
# pointing this script at the C13 run would have raised KeyError inside the mutation,
# which mutate() reports as a non-zero exit, which reads as "refused". A script that cannot
# grade the only run that has ever passed, announcing that it graded it.
POS="C8"; NEG="C9"
if [ -f "$RUN/pair.txt" ]; then
  POS="$(sed -n 5p "$RUN/pair.txt")"
  NEG="$(sed -n 6p "$RUN/pair.txt")"
fi
if [ -z "$POS" ] || [ -z "$NEG" ]; then
  echo "pair.txt under $RUN does not name both controls on lines 5 and 6." >&2
  echo "Refusing to guess which control is the positive." >&2
  exit 2
fi
echo "=== grading the judge on $RUN"
echo "    positive=$POS   negative=$NEG"

# §8.16: check the honest run passes FIRST. If it does not, every refusal below is a
# refusal of something other than the mutation, and the whole pass means nothing.
echo "=== the honest run, which must PASS"
if ! python3 -I "$JUDGE" "$RUN" > "$WORK.honest" 2>&1; then
  echo "ABORTING: the unmutated run does not pass the judge, so nothing below would mean"
  echo "anything. Fix that first (§8.16). Last lines:"
  tail -12 "$WORK.honest"
  exit 1
fi
echo "    passes"

mutate() {            # mutate <name> <python-body-over-'d'>
  local name="$1" body="$2"
  local copy="$WORK/$name"
  rm -rf "$copy"; mkdir -p "$copy"
  cp -r "$RUN"/. "$copy"/ 2>/dev/null
  python3 -I - "$copy/$SCORES" <<PY || return 2
import json, sys
path = sys.argv[1]
with open(path) as handle:
    d = json.load(handle)
rows = {r.get("design"): r for r in d.get("designs") or []}
$body
with open(path, "w") as handle:
    json.dump(d, handle, indent=2)
PY
  if python3 -I "$JUDGE" "$copy" > "$copy/judged.txt" 2>&1; then
    echo "  [GOT THROUGH]  $name"
    echo "      the judge still said MET. Last lines:"
    tail -4 "$copy/judged.txt" | sed 's/^/      /'
    return 1
  fi
  echo "  [refused]      $name -- $(reason "$copy/judged.txt")"
  return 0
}

# The judge refuses in two shapes: an early "NOT MET: <why>" for a tree it cannot read at
# all, and a bulleted list for a run it read and rejected. Both are refusals and both have
# to be quoted here, or a mutation caught by the first shape is recorded with no reason at
# all -- which is a pass that looks like a pass and tells nobody what was caught.
reason() {
  local log="$1"
  grep -m1 -E '^(NOT MET:|  - )' "$log" \
    | sed -e 's/^  - //' -e 's/^NOT MET: //' | cut -c1-110
}

echo
echo "=== mutations, each of which must be REFUSED"
RC=0

# 1. The positive fails. The most direct falsification there is.
mutate positive-iptm-collapses \
  "rows[\"$POS\"][\"iptm\"] = 0.05" || RC=1

# 2. The negative passes. A scorer that cannot see a shuffle is the §8.6 failure.
mutate negative-iptm-as-good-as-positive \
  "rows[\"$NEG\"][\"iptm\"] = rows[\"$POS\"][\"iptm\"]" || RC=1

# 3. The ordering inverts. Separation with the wrong sign is not separation.
mutate separation-inverted \
  "rows[\"$POS\"][\"iptm\"], rows[\"$NEG\"][\"iptm\"] = rows[\"$NEG\"][\"iptm\"], rows[\"$POS\"][\"iptm\"]" \
  || RC=1

# 4. The gating metric goes missing. §8.7: never default a missing metric to a number,
#    and never let an absent one count as separation.
mutate positive-iptm-missing \
  "rows[\"$POS\"].pop(\"iptm\", None)" || RC=1

# 5. A control vanishes entirely. The Phase 1 defect: a record that is absent generates no
#    row to notice.
mutate negative-row-removed \
  "d[\"designs\"] = [r for r in d[\"designs\"] if r.get(\"design\") != \"$NEG\"]" || RC=1

# 6. Two rows for one control. Which one does the judge believe? Neither: it refuses.
mutate duplicate-positive \
  "d[\"designs\"].append(dict(rows[\"$POS\"]))" || RC=1

# 7. A chain is not modelled. Must come back NOT INTERPRETABLE rather than MET -- and
#    rather than NOT MET, which is the conflation this check exists to break. This is
#    also run 2's failure in synthetic form: C8's binder floored at 38.55 is the real
#    measured number, and it is the regression that keeps the judge's narrowing honest
#    -- the POSITIVE's binder must still floor, or the narrowing rescued run 2 after
#    all and was written to do so.
mutate binder-chain-not-modelled \
  "rows[\"$POS\"][\"plddt_binder_in_complex\"] = 38.55" || RC=1

# 8. The POSITIVE's target is not modelled. The chain 8e4772a7 was actually measured on:
#    a receptor at ~42 pLDDT moved the interface PAE floor from 7.6 to 15.4 A.
mutate positive-target-not-modelled \
  "rows[\"$POS\"][\"plddt_target\"] = 38.55" || RC=1

# 9. The NEGATIVE's target is not modelled. This is the exact boundary of the narrowing
#    the judge now makes, so it is the mutation that boundary lives or dies by. The
#    negative's BINDER is exempt, because a scramble is allowed not to fold and that is
#    the intervention. Its TARGET is not exempt, because the target is the same chain with
#    the same alignment as the positive's, so a drop there is the flags' fault. Without
#    this check the narrowing would be untested in the one direction where getting it
#    wrong turns the whole gate into a formality.
mutate negative-target-not-modelled \
  "rows[\"$NEG\"][\"plddt_target\"] = 38.55" || RC=1

# 10. The guard is never shown refusing. Not a mutation of the scores, so it is done here
#    rather than through mutate().
echo "  --- and one that is not about the scores"
COPY="$WORK/refusal-not-shown"
rm -rf "$COPY"; mkdir -p "$COPY"; cp -r "$RUN"/. "$COPY"/ 2>/dev/null
rm -f "$COPY/n3-refuse.log"
if python3 -I "$JUDGE" "$COPY" > "$COPY/judged.txt" 2>&1; then
  echo "  [GOT THROUGH]  refusal-not-shown"
  echo "      the judge said MET without node 03's single_sequence refusal ever being"
  echo "      exercised, which is exactly what §8.13 forbids."
  RC=1
else
  echo "  [refused]      refusal-not-shown -- $(reason "$COPY/judged.txt")"
fi

# ---------------------------------------------------------------------------------------
# THE EXEMPTION, PINNED.
#
# Every check above asks whether the judge still REFUSES. This asks the opposite question,
# and it has to be asked out loud: the judge deliberately does NOT floor the negative's
# binder, and an exemption that nothing tests is indistinguishable from an exemption
# someone will widen back by accident. 8.13 says a guard that cannot be shown refusing is
# not a guard; the mirror of that is that a deliberate permission nobody can show being
# granted is not a decision, it is a gap.
#
# Two things have to hold, not one. The run must still produce a verdict, AND the floored
# negative chain must still be DISCLOSED in the output. An exemption that silently drops
# the observation would give the right verdict while hiding the reason it is weak, which is
# the shape 8.2 is about.
permit() {            # permit <name> <python-body-over-'d'>
  local name="$1" body="$2"
  local copy="$WORK/permit-$name"
  rm -rf "$copy"; mkdir -p "$copy"
  cp -r "$RUN"/. "$copy"/ 2>/dev/null
  python3 -I - "$copy/$SCORES" <<PY || return 2
import json, sys
path = sys.argv[1]
with open(path) as handle:
    d = json.load(handle)
rows = {r.get("design"): r for r in d.get("designs") or []}
$body
with open(path, "w") as handle:
    json.dump(d, handle, indent=2)
PY
  if ! python3 -I "$JUDGE" "$copy" > "$copy/judged.txt" 2>&1; then
    echo "  [WRONGLY REFUSED]  $name"
    echo "      the judge called this uninterpretable. A scrambled negative is allowed"
    echo "      not to fold -- that is the intervention, not an artifact of the flags."
    echo "      Reason given: $(reason "$copy/judged.txt")"
    return 1
  fi
  if ! grep -q "below the 70.0 floor" "$copy/judged.txt"; then
    echo "  [PERMITTED SILENTLY]  $name"
    echo "      the judge returned a verdict without disclosing that a chain was on the"
    echo "      floor. The verdict is then stronger-looking than the run supports."
    return 1
  fi
  echo "  [permitted, and disclosed]  $name"
  return 0
}

echo
echo "=== and one that must be PERMITTED, with its limit disclosed"
permit negative-binder-not-modelled \
  "rows[\"$NEG\"][\"plddt_binder_in_complex\"] = 32.35" || RC=1

echo
echo "=============================================================================="
if [ "$RC" = "0" ]; then
  echo "Every mutation was refused, and the one exemption was granted and disclosed."
  echo "The judge can tell a passing run from a failing one, and can tell a chain that"
  echo "should have modelled from one that was never going to."
  echo "Evidence under $WORK"
else
  echo "AT LEAST ONE MUTATION GOT THROUGH. The Phase 2 verdict means nothing until the"
  echo "judge refuses it. Evidence under $WORK"
fi
echo "=============================================================================="
exit "$RC"
