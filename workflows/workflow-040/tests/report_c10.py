#!/usr/bin/env python3
"""Judge C10: what ColabFold does with a complex input whose ':' is missing.

Four claims. Three gate; the fourth is the MEASUREMENT this control exists to take, and it
is reported rather than gated because its value is not known in advance and a judge that
decided it beforehand would be asserting the thing it was sent to find out.

  1  gates    arm A, the honest colon input, came back as TWO chains of 87 and 108
  2  gates    arm B, the colon deleted, did NOT come back as two chains
  3  gates    neither arm reached the network
  4  RECORDS  did arm B fail loudly, or succeed quietly -- and did it still report an iptm

Claim 1 is first and gates for the reason in §8.16: if the honest arm does not produce a
two-chain complex, arm B producing something else says nothing about the colon, and every
number below would be about a broken control.

Claim 2 gates because a control in which the mutation changed nothing is not a control. If
deleting the colon still yielded two chains, this whole exercise measured nothing.

Claim 4 is the point, and the first version of this docstring got its own premise wrong.
It said node 03's "stated justification" was that a colon-less complex comes back as one
fused chain while still producing an iptm. Node 03 says no such thing -- its message is
only that the colon count is wrong and the file must have changed since node 02 verified
it. The iptm half was my assertion, and the measurement refuted it: the fusion is real and
silent, and no iptm is reported, because ColabFold omits it for a single-chain prediction.

Which is why claim 4 records rather than gates. Whichever way it falls, the consequence for
node 03's guard is stated here rather than left to a reader:

  silent, with an iptm   the guard is LOAD-BEARING: the failure is invisible downstream
  silent, no iptm        the guard still matters, and node 04 would have nulled the metric
  loud error             the guard is defence in depth, and the predictor refuses too

Usage:  python3 -I tests/report_c10.py <RUN_DIR>
"""

import glob
import json
import os
import sys

EXPECTED_CHAIN_LENGTHS = [87, 108]   # barstar target, barnase binder -- from the probe


def chains_in_pdb(path):
    """Residues per chain, by counting CA records. Ordered as they appear in the file."""
    order, counts = [], {}
    with open(path, errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            chain = line[21]
            if chain not in counts:
                order.append(chain)
                counts[chain] = 0
            counts[chain] += 1
    return [(c, counts[c]) for c in order]


def arm(root, name):
    """Everything readable about one arm, without deciding anything about it."""
    out = os.path.join(root, name, "out")
    pdbs = sorted(glob.glob(os.path.join(out, "*_unrelaxed_rank_001*.pdb")))
    scores = sorted(glob.glob(os.path.join(out, "*_scores_rank_001*.json")))
    log_path = os.path.join(root, f"{name}.log")
    log = open(log_path, errors="replace").read() if os.path.exists(log_path) else ""

    record = {
        "pdb": pdbs[0] if pdbs else None,
        "chains": chains_in_pdb(pdbs[0]) if pdbs else None,
        "scores_file": scores[0] if scores else None,
        "iptm": None,
        "ptm": None,
        "log_bytes": len(log),
        "exited_non_zero": "exited non-zero" in log,
        # Deliberately broad: the question is whether it complained at all, and a judge
        # that grepped for one exact string would report "quiet" for a different message.
        "log_mentions_error": any(
            marker in log.lower()
            for marker in ("traceback", "error", "exception", "invalid", "refus")
        ),
    }
    if scores:
        try:
            data = json.load(open(scores[0]))
            record["iptm"] = data.get("iptm")
            record["ptm"] = data.get("ptm")
        except (ValueError, OSError) as error:
            record["scores_unreadable"] = str(error)
    return record


def main(argv):
    if len(argv) != 2:
        print(__doc__)
        return 2
    root = argv[1]

    print("=" * 78)
    print("workflow-040 C10: what ColabFold returns for a colon-less complex")
    print(f"  run: {root}")
    print("=" * 78)

    staged_path = os.path.join(root, "inputs.json")
    if not os.path.exists(staged_path):
        print(f"\nNOT MET: no {staged_path}; the staging step did not run.")
        return 1
    staged = json.load(open(staged_path))

    print("\ninputs, as staged")
    for name in ("A", "B"):
        entry = staged.get(name, {})
        print(f"  {name}  colons={entry.get('colons')}  "
              f"residues={entry.get('residues_excluding_colon')}  "
              f"declared chain lengths={entry.get('declared_chain_lengths')}")

    arms = {name: arm(root, name) for name in ("A", "B")}

    print("\nwhat came back")
    for name in ("A", "B"):
        got = arms[name]
        print(f"  {name}  chains={got['chains']}  iptm={got['iptm']}  "
              f"ptm={got['ptm']}  non-zero exit={got['exited_non_zero']}  "
              f"complained={got['log_mentions_error']}")
        if got["pdb"] is None:
            print(f"      no rank_001 pdb under {name}/out")

    failures = []

    # ---------------------------------------------------------------- claim 1
    a_chains = arms["A"]["chains"]
    a_ok = a_chains is not None and [n for _, n in a_chains] == EXPECTED_CHAIN_LENGTHS
    if not a_ok:
        failures.append(
            f"the honest colon input did not come back as two chains of "
            f"{EXPECTED_CHAIN_LENGTHS}; it gave {a_chains}. Nothing about arm B means "
            "anything until this holds (§8.16)."
        )

    # ---------------------------------------------------------------- claim 2
    b_chains = arms["B"]["chains"]
    b_is_two = b_chains is not None and [n for _, n in b_chains] == EXPECTED_CHAIN_LENGTHS
    if b_is_two:
        failures.append(
            "deleting the colon changed nothing: arm B came back as the same two chains. "
            "A control whose mutation had no effect has measured nothing."
        )

    # ---------------------------------------------------------------- claim 3
    for name in ("A", "B"):
        log_path = os.path.join(root, f"{name}.log")
        text = open(log_path, errors="replace").read() if os.path.exists(log_path) else ""
        if "api.colabfold.com" in text:
            failures.append(
                f"arm {name}'s log mentions api.colabfold.com, so this run cannot be "
                "called offline."
            )

    # ---------------------------------------------------------------- claim 4, recorded
    print("\n" + "-" * 78)
    print("THE MEASUREMENT: what a colon-less complex actually does")
    b = arms["B"]
    if arms["A"]["pdb"] is None and b["pdb"] is None:
        # NEITHER arm produced anything, so this is not evidence about the predictor at
        # all -- it is a run that did not happen. Saying "the predictor refused" here
        # would be §8.3's exact conflation: an empty output directory and a tool that
        # declined are indistinguishable until something tells them apart, and the thing
        # that tells them apart is whether the HONEST arm worked. Caught by running the
        # judge against a staged-but-unpredicted tree before spending any GPU time.
        print("  NOT AVAILABLE: neither arm produced a structure, so nothing here is")
        print("  evidence about the predictor. A run in which the honest arm also")
        print("  produced nothing is a run that did not happen, not a refusal.")
        verdict_note = "no measurement: neither arm produced a structure"
    elif b["pdb"] is None:
        print("  Arm B produced no structure while arm A did. The predictor refused or")
        print("  crashed on the colon-less input specifically.")
        print("  CONSEQUENCE: node 03's guard is defence in depth -- the predictor")
        print("  rejects this input too, so the failure would not have been silent.")
        verdict_note = "arm B produced no structure while arm A did; the predictor refused it"
    else:
        lengths = [n for _, n in (b["chains"] or [])]
        total = sum(lengths)
        print(f"  Arm B produced a structure: {len(lengths)} chain(s), lengths {lengths}, "
              f"{total} residues total.")
        if len(lengths) == 1:
            print("  The two sequences were FUSED into one chain, which is what node 03's")
            print("  refusal message says happens.")
        if b["iptm"] is not None:
            print(f"  And it still reported iptm = {b['iptm']}.")
            print("  CONSEQUENCE: node 03's guard is LOAD-BEARING. The input declares no")
            print("  interface, the predictor invents a single chain, and an interface")
            print("  number comes out anyway -- computable on garbage and ordinary-looking")
            print("  (§8.6). Nothing downstream could tell this from a real complex.")
            verdict_note = (f"arm B fused to {len(lengths)} chain(s) and still reported "
                            f"iptm={b['iptm']}")
        else:
            print("  It reported no iptm.")
            print("  CONSEQUENCE: the guard still matters -- the structure is wrong -- but")
            print("  node 04 would have nulled the interface metrics rather than")
            print("  reporting a plausible one.")
            verdict_note = f"arm B fused to {len(lengths)} chain(s) and reported no iptm"
        if not b["log_mentions_error"]:
            print("  The predictor did not complain. The failure is SILENT.")
    print("-" * 78)

    print("\n" + "=" * 78)
    if failures:
        print("C10: NOT MET")
        for line in failures:
            print(f"  - {line}")
        print("=" * 78)
        return 1
    print("C10: MET")
    print(f"  The control is sound and the measurement is recorded: {verdict_note}.")
    print("  Arm A gave the expected two chains, arm B did not, neither reached the")
    print("  network, and the consequence for node 03's guard is stated above.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
