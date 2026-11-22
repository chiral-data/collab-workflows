#!/usr/bin/env python3
"""Judge Gate B Phase 2: does 040's interface scoring separate a real binder from a shuffle?

Reads tests/run_controls_phase2.sh's output tree. Recomputes nothing -- the numbers come
out of design_scores.json, which is what the workflow actually produced, for the same
reason report_controls.py does it that way: a judge that re-derived the metric would be
agreeing with its own arithmetic rather than with the node's.

What a MET verdict means here
-----------------------------
Only that C8 and C9 separated, in the stated direction, on the stated metric, in one run.
It is not a calibration: the target alignment in an offline run is synthetic, and §11.1
measured that a small synthetic a3m can score worse than no alignment at all, so the
absolute values are not threshold material. The verdict says this itself rather than
relying on a reader to remember it.

And a MET verdict is worth nothing until this judge has been shown refusing. That is
tests/falsify_controls_phase2.sh, which mutates the scores and requires a refusal for each
mutation -- the same pass that found two defects in the Phase 1 judge, neither of which
was visible by reading it (§8.13, §8.16).

Usage:  report_controls_phase2.py <OUT_DIR>
"""
import json
import os
import sys

# The positive and the negative, with what each is. Named here rather than inferred from
# the scores, because a judge that worked out which control was the positive from which one
# scored better would pass any run (§8.4).
POSITIVE = {"C8": "1BRS barnase on barstar, the deposited sequence (NATIVE)"}
NEGATIVE = {"C9": "the same backbone, same residues, order shuffled (SCRAMBLED)"}

# (column, direction, gates). `gates` records whether this metric is allowed to decide the
# verdict. iptm is the interface quantity the multimer model exists to produce, so it
# gates. pae_interaction is reported and does NOT: §7 ships it as a provisional band, no
# control has calibrated it, and a threshold nobody has measured cannot be the thing that
# fails a run.
METRICS = [
    ("iptm", "higher", True),
    ("pae_interaction", "lower", False),
    ("plddt_target", "higher", False),
    ("plddt_binder_in_complex", "higher", False),
]


def load(path, what):
    if not os.path.exists(path):
        return None, f"{what} is not there ({path})"
    try:
        with open(path) as handle:
            return json.load(handle), None
    except (OSError, ValueError) as error:
        return None, f"{what} could not be read: {error}"


def rows_by_design(scores):
    out = {}
    for row in scores.get("designs") or []:
        name = row.get("design")
        if name in out:
            return None, f"design {name!r} appears twice in design_scores.json"
        out[name] = row
    return out, None


def number(row, column):
    """The value, or None. A missing metric is NOT a number (§8.7)."""
    if row is None:
        return None
    value = row.get(column)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    return float(value)


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    out = os.path.abspath(argv[1])

    print("=" * 78)
    print("workflow-040 Tier 2, Phase 2: does interface scoring separate a shuffle?")
    print(f"  run: {out}")
    print("=" * 78)

    scores, why = load(os.path.join(out, "n4", "outputs", "design_scores.json"),
                       "design_scores.json")
    if scores is None:
        print(f"\nNOT MET: {why}")
        print("Node 04 produced no scores, so there is nothing to judge. The node logs "
              "are in this directory.")
        return 1

    rows, why = rows_by_design(scores)
    if rows is None:
        print(f"\nNOT MET: {why}")
        return 1

    # ---------------------------------------------------------------- provenance
    print("\nprovenance")
    manifest, _ = load(os.path.join(out, "n3", "outputs", "prediction_manifest.json"),
                       "prediction_manifest.json")
    alignment, _ = load(os.path.join(out, "n1b", "outputs", "target_alignment.json"),
                        "target_alignment.json")
    if manifest:
        cplx = manifest.get("cplx") or {}
        print(f"  model          {cplx.get('model_type')}")
        print(f"  msa_mode       {cplx.get('msa_mode')}")
        print(f"  recycles       {cplx.get('num_recycle')}")
        print(f"  binder aligned {cplx.get('binder_aligned')}")
        print(f"  host_url       passed={manifest.get('host_url_passed')}")
    else:
        print("  prediction_manifest.json NOT RECORDED")
    if alignment:
        print(f"  target MSA     {alignment.get('n_sequences')} sequence(s) from "
              f"{alignment.get('source')!r}, network_used="
              f"{alignment.get('network_used')}, sha256="
              f"{str(alignment.get('a3m_sha256'))[:16]}")
        for line in alignment.get("target_fusion") or []:
            print(f"  FUSED          {line}")
    else:
        print("  target_alignment.json NOT RECORDED")

    # ---------------------------------------------------------------- the table
    print("\nmeasured")
    header = f"  {'control':<6} {'role':<10}" + "".join(
        f"{column:>26}" for column, _, _ in METRICS)
    print(header)
    absent = []
    for name in list(POSITIVE) + list(NEGATIVE):
        row = rows.get(name)
        if row is None:
            absent.append(name)
        role = "POSITIVE" if name in POSITIVE else "NEGATIVE"
        cells = ""
        for column, _, _ in METRICS:
            value = number(row, column)
            cells += f"{'ABSENT' if value is None else format(value, '.4f'):>26}"
        print(f"  {name:<6} {role:<10}{cells}")
        if row is not None and row.get("null_reasons"):
            for column, reason in (row["null_reasons"] or {}).items():
                print(f"         null {column}: {reason}")

    if absent:
        print(f"\nNOT MET: {absent} produced no row in design_scores.json. A control that "
              "is not there generates no row to notice, which is how the Phase 1 judge "
              "passed a deleted provenance file while refusing a partial one (§8.13).")
        return 1

    # ------------------------------------------------- is this run interpretable at all?
    # 8e4772a7: the interface can be no better than the WORSE-modelled of the two chains.
    # That rule is SYMMETRIC, and until this check existed nothing enforced it on either
    # chain. Both Phase 2 runs to date were capped by a chain nobody had required to be
    # modelled:
    #
    #   run 1 (synthetic alignment)  the TARGET at 38.28 -- §4's single_sequence figure
    #   run 2 (real 3600-seq MSA)    the BINDER at 38.55 -- native barnase, given no
    #                                alignment by §4's design, which is correct for a
    #                                DESIGNED binder and wrong for this control
    #
    # Run 2's pre-registration got this half right: it required plddt_target > 70 and said
    # nothing about the binder, so a symmetric rule was applied to one side. The floored
    # binder was then found by reading the numbers rather than by any check (§8.17 again,
    # and the second time in the same control).
    #
    # So NOT MET is split in two. A run whose chains are not modelled says NOTHING about
    # separation, and reporting that as NOT MET invites exactly the reading that cost this
    # workflow a day -- that 040 cannot score an interface. "This control cannot tell" and
    # "the scorer failed" are different findings with different remedies, so they get
    # different words and different exit codes.
    MODELLED_FLOOR = 70.0
    CHAINS = ("plddt_target", "plddt_binder_in_complex")
    floored = []
    for name in list(POSITIVE) + list(NEGATIVE):
        for column in CHAINS:
            value = number(rows[name], column)
            if value is None or value < MODELLED_FLOOR:
                floored.append((name, column, value))

    if floored:
        print("\n" + "=" * 78)
        print("Gate B Phase 2: NOT INTERPRETABLE -- this control cannot tell")
        # EVERY refusal this judge makes prints as a "  - " bullet, this one included.
        # falsify_controls_phase2.sh quotes the first such line as the reason a mutation
        # was caught. When NOT INTERPRETABLE was first added with an output shape of its
        # own, the falsifier recorded that mutation as refused with an EMPTY reason: a
        # pass that looks like a pass and tells nobody what was caught. The same thing had
        # already happened once, for the early-exit versus bulleted shapes. Two formats
        # were one too many and a third was the point to stop adding greps to the reader
        # and standardise the writer instead.
        for name, column, value in floored:
            shown = "NOT MEASURED" if value is None else f"{value:.2f}"
            print(f"  - {name} {column} is {shown}, below the {MODELLED_FLOOR:.1f} "
                  "floor, so that chain was not modelled and no interface number "
                  "computed against it means anything (8e4772a7).")
        print("  Separation is not reported, because a chain that was not modelled caps "
              "the interface")
        print("  metrics regardless of the design (8e4772a7), so any margin read off them "
              "would be a")
        print("  measurement of the floor. This is NOT a verdict on 040's interface "
              "scoring: it is a")
        print("  statement that this control cannot produce one. See")
        print("  tests/reference/tier2-phase2-controls.md.")
        print("=" * 78)
        return 2

    # ---------------------------------------------------------------- separation
    print("\nseparation -- the worst positive against the best negative")
    verdict_lines = []
    failures = []
    for column, direction, gates in METRICS:
        pos = [number(rows[n], column) for n in POSITIVE]
        neg = [number(rows[n], column) for n in NEGATIVE]
        if any(v is None for v in pos + neg):
            note = (f"  {column:<26} NOT MEASURED on at least one control, so it "
                    "separates nothing")
            print(note)
            if gates:
                failures.append(f"{column} is missing on a control and it is a gating "
                                "metric; a missing metric is not a number (§8.7)")
            continue
        if direction == "higher":
            worst_pos, best_neg = min(pos), max(neg)
            margin = worst_pos - best_neg
        else:
            worst_pos, best_neg = max(pos), min(neg)
            margin = best_neg - worst_pos
        sign = "+" if margin > 0 else ""
        tag = "GATES" if gates else "diagnostic"
        print(f"  {column:<26} positive {worst_pos:>10.4f}   negative {best_neg:>10.4f}"
              f"   margin {sign}{margin:.4f}   [{tag}]")
        verdict_lines.append((column, margin, gates))
        if gates and margin <= 0:
            failures.append(
                f"{column} does not separate: the positive is {worst_pos:.4f} and the "
                f"negative is {best_neg:.4f}, margin {margin:+.4f}. 040 cannot tell a "
                "real binder from the same backbone with its residues shuffled."
            )

    # ---------------------------------------------------------------- the refusal
    print("\nthe guard, shown refusing")
    refuse_log = os.path.join(out, "n3-refuse.log")
    if not os.path.exists(refuse_log):
        failures.append("n3-refuse.log is absent, so node 03's single_sequence refusal "
                        "was never exercised (§8.13)")
        print("  n3-refuse.log ABSENT")
    else:
        text = open(refuse_log, errors="replace").read()
        said = "single_sequence" in text
        nonzero = "exit=0 " not in text and "exit=0\n" not in text
        print(f"  refused single_sequence on the complex path: names it={said}, "
              f"non-zero exit={nonzero}")
        if not (said and nonzero):
            failures.append(
                "node 03 did not refuse msa_mode=single_sequence on the complex path. "
                "colabfold clamps max_seq to 1 in that mode, so the alignment is "
                "discarded while config.json, the re-saved a3m and coverage.png all "
                "record that it was used (§11.1)."
            )

    # ---------------------------------------------------------------- verdict
    print("\n" + "=" * 78)
    if failures:
        print("Gate B Phase 2: NOT MET")
        for line in failures:
            print(f"  - {line}")
        print("=" * 78)
        return 1

    gating = [(c, m) for c, m, g in verdict_lines if g]
    print("Gate B Phase 2: MET -- " + ", ".join(
        f"{c} separates by {m:+.4f}" for c, m in gating))
    print()
    print("What this does NOT say:")
    print("  - Nothing about thresholds. The target alignment in an offline run is")
    print("    synthetic, and §11.1 measured a small synthetic a3m scoring WORSE than no")
    print("    alignment at all, so these absolute values are not calibration inputs.")
    print("    pae_interaction is printed as a diagnostic for that reason and does not")
    print("    gate here, and does not gate in the workflow either (§7).")
    print("  - Nothing about designed binders. C8 is a native pair with deep MSAs on both")
    print("    sides; §4 records 029's 0.100 -> 0.860 as an UPPER bound for the designed")
    print("    case, and nothing in this repo has measured the designed case.")
    print("  - Nothing until this judge has been shown refusing:")
    print("    bash tests/falsify_controls_phase2.sh <this directory>")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
