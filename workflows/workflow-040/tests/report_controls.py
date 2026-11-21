#!/usr/bin/env python3
"""Read a Tier 2 control run and say whether Gate B's Phase 1 condition holds.

Gate B is not "the nodes ran". It is one quantitative claim, and it is the only claim in
either workflow that can falsify the SCORE rather than the plumbing:

    the worst positive must separate from the best negative, on sc_rmsd_A.

C12 is a DESIGNED binder on its own backbone and must pass the thresholds. C3 and C4 are a
native backbone with the sequence shuffled and with the sequence replaced by poly-alanine,
and must fail them. If they do not separate, 040 is wrong and does not merge.

AMENDED after the first run, and the amendment is the point of the file:

  - The positive control used to be C1/C2, NATIVE complexes from the PDB. They failed, at
    49.33 and 43.09 pLDDT. The control was wrong, not the scorer: spec §4 fixes
    msa_mode = single_sequence because a de novo binder has no homologs, and a native
    immunoglobulin domain is the opposite -- it needs the alignment it is not being given.
    A native positive at single_sequence measures AF2, not 040. They are DIAGNOSTIC here
    now, reported with their numbers but not gating, and §4's Phase 2 MSA arms are
    measured against them.

  - plddt_binder is no longer required to separate. C4 is 115 alanines and scores 95.81,
    ABOVE the real design's 95.05, and clears the 80.0 threshold. Confidence in *a* fold
    is not evidence of *the designed* fold. It may reject a design, it may not admit one,
    so requiring it to separate would be requiring something untrue. It is still printed,
    with its margin, because the negative number is the finding.

Both amendments are recorded in tests/reference/tier2-phase1-controls.md with the run that
produced them, and in spec §8.4 (a control that passes review for the wrong reason) and
§8.6 (a gating metric that looks excellent on garbage).

The margin is printed with its sign both ways round, because a negative margin is the
interesting result and rounding it into a verdict would hide it. Everything is read from
design_scores.json, the inventory and the provenance files this run wrote; nothing is
recomputed here.

Run:  python3 -I tests/report_controls.py <control-run-dir>
"""
import json
import os
import sys

# Which control is a positive and which is a negative. Stated here rather than inferred
# from the scores, which would be circular.
POSITIVE = {
    "C12": "an idealised designed helix carrying a designed sequence (DESIGNED)",
}
NEGATIVE = {
    "C3": "C1's backbone, sequence permuted (composition identical)",
    "C4": "C1's backbone, sequence all alanine",
}
# Reported, never gated. See the amendment note above: a native chain predicted at
# single_sequence measures the prediction mode, so it cannot pass or fail 040.
DIAGNOSTIC = {
    "C1": "4ZQK PD-1, native, single_sequence baseline for §4's MSA arms",
    "C2": "1BRS barnase, native, numbered from 3, same baseline",
}
# Lower is better for RMSD, higher is better for pLDDT.
#   gates=True  -> must separate, or Gate B is not met
#   gates=False -> measured and printed, but a negative margin is a known finding
GATING = [
    ("sc_rmsd_A", "max_sc_rmsd_a", "lower", True),
    ("plddt_binder", "min_plddt_binder", "higher", False),
]


def load(path, what):
    if not os.path.exists(path):
        print(f"report_controls: no {what} at {path}", file=sys.stderr)
        return None
    with open(path) as handle:
        return json.load(handle)


def margin(positives, negatives, metric, direction):
    """(margin, worst positive, best negative). Positive margin means they separate."""
    pos = [(d, v) for d, v in positives if v is not None]
    neg = [(d, v) for d, v in negatives if v is not None]
    if not pos or not neg:
        return None, None, None
    if direction == "lower":
        worst = max(pos, key=lambda p: p[1])
        best = min(neg, key=lambda p: p[1])
        return best[1] - worst[1], worst, best
    worst = min(pos, key=lambda p: p[1])
    best = max(neg, key=lambda p: p[1])
    return worst[1] - best[1], worst, best


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    root = sys.argv[1]
    scores = load(os.path.join(root, "n4", "outputs", "design_scores.json"),
                  "design_scores.json")
    if scores is None:
        return 1
    rows = {r["design"]: r for r in scores["designs"]}
    thresholds = scores.get("thresholds") or {}

    width = 78
    print("=" * width)
    print("workflow-040 Tier 2, Phase 1 controls")
    print(f"  run:        {root}")
    predictor = scores.get("predictor") or {}
    print(f"  predictor:  model_type={predictor.get('model_type')} "
          f"msa_mode={predictor.get('msa_mode')} "
          f"num_models={predictor.get('num_models')} "
          f"num_recycle={predictor.get('num_recycle')}")
    print(f"  network:    host_url_passed={predictor.get('host_url_passed')} "
          f"expected={predictor.get('network_expected')!r} "
          f"synthetic={predictor.get('synthetic')}")
    print(f"  thresholds: {thresholds} "
          f"(calibrated: {scores.get('thresholds_are_calibrated')})")
    print("=" * width)

    missing = [d for d in list(POSITIVE) + list(NEGATIVE) + list(DIAGNOSTIC)
               if d not in rows]
    if missing:
        print(f"MISSING from design_scores.json: {missing}")
        print(f"  it holds {sorted(rows)}")

    print(f"\n{'':4} {'role':10} {'len':>4} {'sc_rmsd_A':>10} {'pLDDT bnd':>10} "
          f"{'pLDDT tgt':>10} {'passes':>7}  failed on")
    for design in sorted(rows):
        row = rows[design]
        role = ("positive" if design in POSITIVE
                else "negative" if design in NEGATIVE
                else "diagnostic" if design in DIAGNOSTIC else "-")

        def show(key, places=2):
            value = row.get(key)
            return "-" if value is None else f"{value:.{places}f}"

        print(f"{design:4} {role:10} {row.get('binder_length', '-'):>4} "
              f"{show('sc_rmsd_A', 4):>10} {show('plddt_binder'):>10} "
              f"{show('plddt_target'):>10} "
              f"{('yes' if row.get('passes') else 'no'):>7}  "
              f"{row.get('failed_thresholds') or ''}"
              + (f"  [{(row.get('reason') or '')[:60]}]" if row.get("reason") else ""))

    # Thresholds first: the claim each control makes on its own.
    print("\nper-control verdict against the thresholds")
    verdicts = {}
    for design, blurb in list(POSITIVE.items()) + list(NEGATIVE.items()):
        row = rows.get(design)
        if row is None:
            verdicts[design] = "MISSING"
            print(f"  {design}  MISSING   {blurb}")
            continue
        wanted = "pass" if design in POSITIVE else "fail"
        got = "pass" if row.get("passes") else "fail"
        verdicts[design] = "ok" if got == wanted else "WRONG"
        print(f"  {design}  {verdicts[design]:7}  must {wanted}, did {got}  -- {blurb}")

    # The diagnostics are printed in the same place so nobody has to notice they are
    # absent, but they are deliberately not scored. A control that cannot be wrong must
    # not be allowed to look like a control that passed (§8.3: a silent nothing reads the
    # same as a success).
    if DIAGNOSTIC:
        print("\ndiagnostic, NOT gated -- these cannot pass or fail 040")
        for design, blurb in DIAGNOSTIC.items():
            row = rows.get(design)
            if row is None:
                print(f"  {design}  absent    {blurb}")
                continue
            plddt = row.get("plddt_binder")
            rmsd = row.get("sc_rmsd_A")
            print(f"  {design}  plddt_binder="
                  f"{'-' if plddt is None else format(plddt, '.2f')}"
                  f"  sc_rmsd_A={'-' if rmsd is None else format(rmsd, '.4f')}"
                  f"  -- {blurb}")
        print("  Native chains at msa_mode=single_sequence measure the prediction mode,")
        print("  not the scorer (spec §4, §8.4). They are the baseline Phase 2's")
        print("  target-only and full-MSA arms are compared against.")

    # Then the separation, which is the part a threshold cannot tell you.
    print("\nseparation between the worst positive and the best negative")
    print("  A positive margin means the two groups do not overlap on that metric. The")
    print("  thresholds could be wrong and the separation still real; the separation could")
    print("  be absent while every threshold happens to be met. Both are reported.")
    separations = {}
    for metric, key, direction, gates in GATING:
        gap, worst, best = margin(
            [(d, rows[d].get(metric)) for d in POSITIVE if d in rows],
            [(d, rows[d].get(metric)) for d in NEGATIVE if d in rows],
            metric, direction,
        )
        separations[metric] = (gap, gates)
        if gap is None:
            print(f"  {metric:14} not measurable: a group has no value"
                  + ("" if gates else "   [not gated]"))
            continue
        arrow = "lower is better" if direction == "lower" else "higher is better"
        print(f"  {metric:14} margin {gap:+.4f}  ({arrow}; threshold "
              f"{key}={thresholds.get(key)})"
              + ("" if gates else "   [NOT GATED, see below]"))
        print(f"  {'':14} worst positive {worst[0]}={worst[1]:.4f}, "
              f"best negative {best[0]}={best[1]:.4f}")
        if not gates and gap <= 0:
            print(f"  {'':14} this margin is NEGATIVE and that is the expected, recorded")
            print(f"  {'':14} finding: {best[0]} encodes nothing and the predictor is")
            print(f"  {'':14} confident about it anyway. {metric} may reject a design,")
            print(f"  {'':14} never admit one (spec §8.6).")

    # What the run itself says about how it was produced, because a control set measured
    # with the wrong settings is a measurement of the settings.
    provenance = []
    absent = []
    for design in sorted(set(POSITIVE) | set(NEGATIVE) | set(DIAGNOSTIC)):
        got = load(os.path.join(root, "g", f"{design}_provenance.json"), "provenance")
        if got:
            provenance.append((design, got))
        else:
            # An absent record has to be tracked, not merely skipped. Skipping it was the
            # hole falsify_controls.sh found: a PARTIAL provenance refused the gate while
            # a DELETED one sailed through, because nothing downstream iterated over the
            # controls that were supposed to have one (§8.3 again -- the absence was the
            # thing, and absence produced no row to notice).
            absent.append(design)
    incomplete = []
    if provenance:
        print("\nprovenance")
        for design, got in provenance:
            # Every field is read with .get and a visible marker for what is absent.
            #
            # This printer used to index the keys directly, and C12's provenance -- the
            # one record written by hand in run_controls.sh rather than by
            # make_controls.py -- does not carry the residue fields. The KeyError landed
            # AFTER the table and BEFORE the verdict, so the tee'd VERDICT.txt of the
            # second control run is a file that reads like a complete report and contains
            # no verdict (§8.3: the failure was visible only in the exit code, and the
            # artifact looked fine). Missing provenance is now reported as missing.
            binder = got.get("binder") or {}
            sha = got.get("source_sha256") or "?"
            residues = binder.get("n_residues")
            span = (f"{residues} residues, ids {binder.get('first_id', '?')}.."
                    f"{binder.get('last_id', '?')}, "
                    f"{len(binder.get('gaps') or [])} gap(s)"
                    if residues is not None else "residue span NOT RECORDED")
            print(f"  {design}  {got.get('source_file', '?')} "
                  f"sha256={sha[:16]}  binder auth "
                  f"{binder.get('auth_asym_id', '?')} -> "
                  f"{binder.get('label_asym_id_written', '?')}, {span}")
            print(f"       {got.get('transform_note', '(no transform note)')}")
            if residues is None:
                incomplete.append(design)
        if incomplete:
            print(f"  NOTE: incomplete provenance for {incomplete} -- the controls ran, "
                  f"but what\n        went into them is not fully recorded (§8.14).")
    if absent:
        print(f"\n  NOTE: no provenance record at all for {absent}")

    print("\n" + "=" * width)
    wrong = [d for d, v in verdicts.items() if v != "ok"]
    unseparated = [m for m, (g, gates) in separations.items()
                   if gates and (g is None or g <= 0)]
    # The gate rests on these controls being what they claim, so for a control that
    # gates, provenance that is absent is not a cosmetic problem (§8.14).
    unprovenanced = sorted(d for d in set(incomplete) | set(absent)
                           if d in POSITIVE or d in NEGATIVE)
    if wrong:
        print(f"Gate B Phase 1: NOT MET -- wrong verdict for {wrong}")
    if unseparated:
        print(f"Gate B Phase 1: NOT MET -- no separation on {unseparated}")
    if unprovenanced:
        print(f"Gate B Phase 1: NOT MET -- gating control(s) {unprovenanced} ran without")
        print("recorded provenance, so the right numbers cannot be shown to have come")
        print("from the right input.")
    if not wrong and not unseparated and not unprovenanced:
        print("Gate B Phase 1: MET -- every positive passes, every negative fails, and")
        print("every gating metric separates. This says the scorer is not broken and is")
        print("not merely generous: a designed sequence recovers its own backbone to under")
        print("1 A while a scramble and a poly-alanine of the same backbone do not.")
        print()
        print("It does NOT say the scorer ranks real designs well. The positive is a")
        print("synthetic idealised helix -- the easiest designed positive there is -- and")
        print("the negatives are easy negatives. The hard case is a design that looks fine")
        print("and does not bind, and nothing in this repo can produce one. Nor is any of")
        print("it experimental: C12 recovering its backbone says its sequence encodes its")
        print("structure, and nothing about binding. See tests/README.md, 'What this suite")
        print("does not establish', and tests/reference/tier2-phase1-controls.md.")
    print("=" * width)
    return 0 if not wrong and not unseparated and not unprovenanced else 1


if __name__ == "__main__":
    sys.exit(main())
