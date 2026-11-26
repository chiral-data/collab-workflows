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
# The pairs this judge knows, keyed by the token run_controls_phase2.sh writes to pair.txt.
# Named here rather than inferred from the scores, because a judge that worked out which
# control was the positive from which one scored better would pass any run (8.4).
PAIRS = {
    "c8": (
        {"C8": "1BRS barnase on barstar, the deposited sequence (NATIVE binder)"},
        {"C9": "the same backbone, same residues, order shuffled (SCRAMBLED)"},
    ),
    "c13": (
        {"C13": "7XYQ DBL1_03 on human PD-L1, the deposited sequence (DESIGNED binder)"},
        {"C14": "the same backbone, same residues, order shuffled (SCRAMBLED)"},
    ),
}

# What a MET verdict explicitly does not establish, per pair. Kept beside PAIRS so that
# adding a pair without deciding its limits is a KeyError at the verdict rather than the
# previous pair's limits printed under a different experiment.
CAVEATS = {
    "c8": [
        "  - Nothing about thresholds. The target alignment in an offline run is",
        "    synthetic, and \u00a711.1 measured a small synthetic a3m scoring WORSE than no",
        "    alignment at all, so these absolute values are not calibration inputs.",
        "    pae_interaction is printed as a diagnostic for that reason and does not",
        "    gate here, and does not gate in the workflow either (\u00a77).",
        "  - Nothing about designed binders. C8 is a native pair with deep MSAs on both",
        "    sides; \u00a74 records 029's 0.100 -> 0.860 as an UPPER bound for the designed",
        "    case, and nothing in this repo has measured the designed case.",
    ],
    "c13": [
        "  - Nothing about thresholds, and the run says the OPPOSITE of calibration: C13",
        "    is a crystallised, deposited, published binder and 040's shipped gating",
        "    thresholds REJECT it (sc_rmsd_A 5.141 > 2.0; plddt_binder 74.4 < 80.0), while",
        "    the two it refuses to gate on classify both controls correctly. That is a",
        "    measured false negative, it is Phase 3's problem, and separation here does",
        "    not license any cutoff. n=1 pair, and C13 clears 10.0 A by 0.145.",
        "  - Little about designed binders in general. This pair DOES measure the designed",
        "    case for the first time (plddt_binder_in_complex 88.92 with no alignment,",
        "    against C12's 95.05), which retires the C8-era caveat. But it is ONE deposited",
        "    binder at 124 residues where C12 was 60, and \u00a74's 0.860 upper bound was never",
        "    a prediction that a designed binder would reach it. C13's 0.84 landing there",
        "    was pre-registered as explicitly UNKNOWN and is a result, not a confirmation.",
        "  - Nothing about targets without homologs. C13's target is human PD-L1 with a",
        "    4374-sequence alignment and modelled at 90.30. A target whose MSA is shallow",
        "    puts the interface back on the floor (8e4772a7), and that case is unmeasured.",
    ],
}

# The default is the pair this judge was written for, so an older run directory with no
# pair.txt in it still reads correctly rather than erroring or, worse, being judged against
# the wrong control names.
DEFAULT_PAIR = "c8"


def pair_of(out):
    """Which pair this run is, read from the run itself.

    run_controls_phase2.sh writes pair.txt as its first act after choosing. Reading it back
    rather than taking it as an argument means the judge cannot be pointed at a run and
    told it is something else -- the run says what it is, and a run that predates pair.txt
    is the default pair, which is what it was.
    """
    path = os.path.join(out, "pair.txt")
    if not os.path.exists(path):
        return DEFAULT_PAIR, f"no pair.txt; assuming {DEFAULT_PAIR!r} (a pre-pair.txt run)"
    token = (open(path).readline() or "").strip().lower()
    if token not in PAIRS:
        return None, (f"pair.txt says {token!r}, which this judge does not know. "
                      f"Known: {sorted(PAIRS)}. Refusing to guess which control is the "
                      "positive.")
    return token, f"pair {token!r}, from pair.txt"

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

    token, note = pair_of(out)
    if token is None:
        print(f"\nNOT MET: {note}")
        return 1
    POSITIVE, NEGATIVE = PAIRS[token]
    print(f"\n{note}")
    for name, description in list(POSITIVE.items()) + list(NEGATIVE.items()):
        role = "POSITIVE" if name in POSITIVE else "NEGATIVE"
        print(f"  {name:<5} {role:<9} {description}")

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

    # WHICH CHAINS THE FLOOR APPLIES TO, AND WHY IT IS NOT SYMMETRIC OVER THE ROLES.
    #
    # 8e4772a7 is a finding about when an interface number is UNINTERPRETABLE, and the
    # question it actually asks is whether a chain's low confidence is an artifact of the
    # FLAGS or a property of the SEQUENCE. That distinction is not symmetric over roles:
    #
    #   positive's target   should model, it has a real MSA                   -> floor it
    #   positive's binder   should model, a designed sequence folds unaided
    #                       (C12 measured 95.05, C13 measured 88.92)          -> floor it
    #   negative's target   should model, same chain and same MSA as the
    #                       positive's, so a drop here is the flags' fault     -> floor it
    #   negative's binder   should NOT model. scrambling destroys the fold,
    #                       and that IS the intervention                  -> DO NOT floor
    #
    # C8/C9 are the contrast that makes this concrete. Their floored chain was native
    # barnase, a protein that folds perfectly well in reality and was floored only because
    # §4 withholds its alignment: an artifact of the flags, and the floor fires on it
    # correctly. A scramble does not fold IN THE WORLD, not merely under these flags, so
    # requiring it to model is requiring the negative to stop being a negative. Phase 1
    # settled this by precedent without anyone noticing: C3 (scrambled, sc_rmsd 18.83 A)
    # and C4 (115 alanines) are both accepted negatives, and C3 is this same intervention.
    #
    # THIS IS A NARROWING MADE AFTER SEEING A RUN, which is the suspect direction (§8.17),
    # so the two things that keep it honest are written here in the code rather than left
    # in the doc where the next reader of this function will not look:
    #
    #   1. The pre-registration committed in 1a5eb76a, EIGHT HOURS BEFORE the run, already
    #      said that a well-modelled negative is expected and that what matters is iptm and
    #      not the negative's pLDDT. The prose specification and this executable one
    #      contradicted each other in one repo state. The run arbitrated between them; it
    #      did not invent the question, and nothing here was written to fit a number.
    #   2. Run 2's floored chain is C8's plddt_binder_in_complex at 38.55, and C8 is the
    #      POSITIVE, so this narrowing provably does NOT rescue the run that motivated the
    #      floor in the first place. A gate change that cannot save the failure behind it
    #      is not a rescue of that failure. falsify_controls_phase2.sh checks both the
    #      chains that must still floor and that run 2 still comes back NOT INTERPRETABLE.
    #
    # What the narrowing does NOT do is make the negative a good negative. A floored binder
    # on the negative side confounds "does not bind" with "does not fold", so it is carried
    # forward as a stated LIMIT on the verdict rather than deleted from the output.
    FLOORED_CHAINS = {
        "positive": ("plddt_target", "plddt_binder_in_complex"),
        "negative": ("plddt_target",),
    }
    CONFOUNDING_CHAINS = {"negative": ("plddt_binder_in_complex",)}

    roles = ([(name, "positive") for name in POSITIVE]
             + [(name, "negative") for name in NEGATIVE])

    floored = []
    confounded = []
    for name, role in roles:
        for column in FLOORED_CHAINS[role]:
            value = number(rows[name], column)
            if value is None or value < MODELLED_FLOOR:
                floored.append((name, column, value))
        for column in CONFOUNDING_CHAINS.get(role, ()):
            value = number(rows[name], column)
            if value is None or value < MODELLED_FLOOR:
                confounded.append((name, column, value))

    if floored:
        print("\n" + "=" * 78)
        print("Gate B Phase 2: NOT INTERPRETABLE -- this control cannot tell")
        # EVERY refusal this judge makes prints as a "  - " bullet, this one included,
        # AND NOTHING ELSE MAY USE THAT SHAPE. The second half is not decoration: the
        # "limits on whatever verdict follows" block below was added as "  - " bullets and
        # immediately made the falsifier attribute every mutation of run 3 to the honest
        # run's own floored negative. One direction of this rule had been written down for
        # two sessions; the other had not, and that is the half that broke.
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

    # ------------------------------------------------- limits carried into the verdict
    # Printed for BOTH verdicts. The negative's binder being on the floor no longer stops
    # the run being read, but it changes what the reading is worth, and a limit that only
    # appears when the news is bad is not a limit, it is an excuse.
    if confounded:
        print("\nlimits on whatever verdict follows")
        # NOT a "  - " bullet. That shape is RESERVED for "this is why the run was
        # refused", and falsify_controls_phase2.sh reads the first one as the reason a
        # mutation was caught. Printed as a bullet, these lines made every mutation of
        # run 3 record as refused because of C14's floored binder: the honest run's own
        # property, identical in every mutated copy, and nothing to do with the mutation.
        # The refusals were real and the reasons on the record were fiction, which is §8.4
        # in the reporting layer rather than the checking layer.
        for name, column, value in confounded:
            shown = "NOT MEASURED" if value is None else f"{value:.2f}"
            print(f"  limit: {name} {column} is {shown}, below the "
                  f"{MODELLED_FLOOR:.1f} floor.")
        print("    This does NOT make the run uninterpretable: a scramble is allowed not "
              "to fold,")
        print("    and that is the intervention rather than an artifact of the flags. It "
              "DOES")
        print("    confound the negative, whose low interface numbers now mean \"does not "
              "bind\"")
        print("    AND \"does not fold\" together. So this is an EASY negative, and part "
              "of any")
        print("    margin below is the chimera being badly modelled on both sides "
              "(8e4772a7:")
        print("    a scrambled partner cost C13's target 10.3 points of its own pLDDT). "
              "The hard")
        print("    negative, a design that models well and does not bind, does not exist "
              "in this")
        print("    repo and this run does not create one.")

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
    # Per-pair, because these sentences are claims about the controls and the C8 set are
    # FALSE of C13. "Nothing in this repo has measured the designed case" was accurate
    # while C8/C9 were the only pair; printed underneath a C13 verdict it would be the
    # judge contradicting the run it is judging. A caveat that is wrong is worse than no
    # caveat, for the same reason a check that passes for the wrong reason is (§8.4).
    for line in CAVEATS[token]:
        print(line)
    print("  - Nothing until this judge has been shown refusing:")
    print("    bash tests/falsify_controls_phase2.sh <this directory>")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
