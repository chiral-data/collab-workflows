#!/usr/bin/env python3
"""Write the prediction inputs: one design per file, always.

Contract
--------
Reads   design_inventory.json
Writes  mono_<design>.fasta      the binder sequence alone, exactly one record
        prediction_inputs.json   what was written, and the facts node 04 needs

Validation happens here because it is free here. A bad alphabet, a length outside bounds
or a missing chain costs nothing to catch on CPU and costs GPU minutes to catch inside
the predictor -- which is why workflow-037 made spec validation its own CPU node
(37c23bb5).

One design per file, never a batch
----------------------------------
Not a style preference. Workflow-026 takes a multi-record FASTA without complaining and
then produces a confident report about a design that does not exist: its 02_predict
copies all ColabFold output flat into one directory, and its 03_analyze_confidence globs
*_scores_rank_*.json across every design and takes job_id from whichever file sorts
first. Three designs at five models makes the reported "top 5" ranks 1,1,1,2,2 drawn from
three different designs, filed under one design's name (818954e9).

So: one record per file, one invocation per file, one output directory per design, and the
mapping back by filename prefix -- never by directory listing order.

Record ids must be filesystem-safe, because ColabFold names its outputs after them
(54bf21b6). The real sample confirms the shape:
sequences_unrelaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb.

The record id is the design name and NOTHING else -- no biosafety tag, no run id, no
length. 037 can afford a tag in its exported FASTA because nothing downstream of it
reads the file back by name; here the id becomes a filename prefix that node 03 writes
and node 04 parses, so anything appended to it has to be stripped by two more nodes that
would then each have their own idea of the format. The biosafety notice travels in
prediction_inputs.json and in node 05's export instead (7709809b).

prediction_inputs.json, per design
----------------------------------
    design              keyed on, everywhere
    path                the file written, basename only
    mode                "mono" in Phase 1
    record_id           what ColabFold will name its outputs after
    chain_order         ["binder"] in Phase 1; ["target", "binder"] in Phase 2
    chain_lengths       parallel to chain_order
    expected_chains     how many chains the prediction must come back with
    expected_residues   sum(chain_lengths); node 04 asserts len(plddt) equals this

chain_order and chain_lengths exist for one reason: the scores file carries NO chain
boundaries. Measured first-hand, not inferred -- see tests/reference/: the real monomer
file carries plddt as a flat list[35] and pae as 35x35, the real complex file carries
plddt as a flat list[63] spanning both chains and pae as 63x63, and neither names a chain
anywhere. So node 04 must slice the per-residue array by position, and the position
depends on the order chosen HERE. This node decides it, so this node records it.

In Phase 1 the binder is the whole array and the slice is trivial. Record it anyway. When
Phase 2 adds TARGET:BINDER the binder becomes the TAIL -- measured: a 34+29 complex
written target-first came back with 34 residues in chain A and 29 in chain B -- and
workflow-029 slices plddt[:binder_length] because its binder is chain A
(04_fold_and_filter/script.py:207). A node 04 that inherits that assumption reports the
target's confidence as the binder's, which on a well-modelled target makes every design
look good.

Colon rules, both directions
----------------------------
A mono input must contain NO colon: a colon makes the predictor treat one design as two
chains. A Phase 2 complex input must contain EXACTLY ONE, with both halves non-empty --
zero silently predicts a monomer and reports its confidence as an interface, which is
026's bug, and two predicts three chains. 026 validates sequences against
ACDEFGHIKLMNPQRSTVWYX and so rejects the only syntax that encodes a complex (d3d707ed);
040 must accept it in the right file and refuse it in the wrong one.

The check is on the bytes written, not on the variable intended to be written. A guard
that reads the value it just assigned is asserting the assignment, not the file (§8.4).

Skips, with a warning rather than silently
------------------------------------------
A design whose sequence is empty or all-unknown is skipped rather than written as junk,
and unknown residues are counted: ColabFold folds X as a sequence gap, which depresses
confidence for reasons that have nothing to do with the design (54bf21b6). Node 01 marks
those designs; this node honours the mark rather than re-deciding it.

Re-running
----------
Must be safe or refused. Workflow-029's chain swap was not idempotent, so re-running over
a non-empty output directory corrupted backbones it had already fixed -- and it fired for
real, during a GPU smoke test, on a directory left by an earlier run (4d1a523e).

Rewriting this node's own outputs is safe, because they are derived from the inventory and
nothing else: same inventory, same bytes. What is NOT safe is a mono_*.fasta left behind
for a design that is not in this run. Node 03 globs mono_*.fasta, so a stale file is
predicted, and node 04 then scores a design the inventory has never heard of. That is
refused, and the refusal leaves every existing file exactly as it was.

Fails when
----------
  * a binder sequence carries a residue outside the accepted alphabet
  * a binder is too short for a superposition to mean anything
  * score_interface is true (Phase 2 is not built -- say so, do not ignore it)
  * the output directory holds a mono_*.fasta for a design not in this run
  * the file it just wrote does not pass the colon and record-count checks
"""
import json
import math
import os
import sys

# 026's alphabet, in 026's order, including X. X is accepted and counted, not rejected: a
# design carrying modified residues is a real design, and refusing it here would discard it
# where recording it costs nothing (54bf21b6). The order is kept because the string is
# quoted back in the refusal and recorded in prediction_inputs.json, and a reader comparing
# it against 026's source should see the same string rather than a sorted paraphrase.
ALPHABET_STRING = "ACDEFGHIKLMNPQRSTVWYX"
ALPHABET = set(ALPHABET_STRING)

# A Kabsch superposition needs three non-collinear points to be defined at all, and four
# before the fit stops being exactly determined by its own inputs. An RMSD over fewer than
# four alpha carbons is a number without a meaning, so it is refused rather than reported.
# This is a geometric floor, not a design-quality threshold -- 037's binder_length default
# is 50-120 and nothing here second-guesses it.
MIN_BINDER_RESIDUES = 4

# Not a limit, a warning line. The only timings this workflow has measured are 35 residues
# in 61 s and a 63-residue complex in 80 s on an RTX 3060 Ti (tests/reference/), which says
# nothing about where memory runs out. AF2 at this size is practically bounded by VRAM
# rather than by anything in this node, so the node says so once and proceeds; §11.5 is the
# open question, and predict_timeout_s is the thing that actually stops a runaway.
LONG_SEQUENCE_WARN_RESIDUES = 1000

MONO_PREFIX = "mono_"
FASTA_SUFFIX = ".fasta"


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def numeric_param(env_name, display_name, default, cast):
    """Read a numeric parameter, failing the same way every other bad parameter does."""
    raw = os.environ.get(env_name)
    if raw is None or raw.strip() == "":
        return default
    try:
        value = cast(raw)
    except ValueError:
        kind = "a whole number" if cast is int else "a number"
        fail(f"{display_name} must be {kind}, got '{raw}'.")
    if not math.isfinite(value):
        fail(
            f"{display_name} must be a finite number, got '{raw}', "
            f"which parses as {value}."
        )
    return value


def bool_param(env_name, display_name, default):
    raw = os.environ.get(env_name)
    if raw is None or raw.strip() == "":
        return default
    lowered = raw.strip().lower()
    if lowered in ("true", "1", "yes", "on"):
        return True
    if lowered in ("false", "0", "no", "off"):
        return False
    fail(f"{display_name} must be true or false, got '{raw}'.")


def load_inventory(path):
    if not os.path.exists(path):
        fail(
            f"No '{os.path.basename(path)}' in ./inputs/. This node reads node 01's "
            "inventory, which is where the designed chain was resolved; it does not parse "
            "structures itself, so there is no fallback and there should not be one -- a "
            "second parser is a second chance to pick a different chain (§8.1)."
        )
    try:
        with open(path) as handle:
            inventory = json.load(handle)
    except Exception as exc:
        fail(f"Could not read '{path}': {exc}")
    if not isinstance(inventory.get("designs"), list):
        fail(f"'{path}' has no 'designs' list; node 01 did not write it, or it is a stub.")
    return inventory


def check_rerun(out_dir, expected_names):
    """Refuse a stale mono_*.fasta before writing anything at all.

    Order matters: the refusal has to happen before the first write, or a run that is
    refused has still half-replaced the directory it refused to touch (4d1a523e).
    """
    if not os.path.isdir(out_dir):
        return
    existing = sorted(
        name for name in os.listdir(out_dir)
        if name.startswith(MONO_PREFIX) and name.endswith(FASTA_SUFFIX)
    )
    expected = {f"{MONO_PREFIX}{name}{FASTA_SUFFIX}" for name in expected_names}
    stale = [name for name in existing if name not in expected]
    if stale:
        fail(
            f"./outputs/ already holds {len(stale)} prediction input(s) for design(s) "
            f"that are not in this run: {', '.join(stale)}. Node 03 globs "
            f"{MONO_PREFIX}*{FASTA_SUFFIX}, so these would be predicted, and node 04 "
            "would then score designs the inventory has never heard of. Nothing has been "
            "written and nothing has been removed -- clear ./outputs/ or run in a fresh "
            "directory. (029's node 04 re-ran over a populated directory during a GPU "
            "smoke test and swapped already-swapped chains back: 4d1a523e.)"
        )
    if existing:
        print(
            f"NOTE: ./outputs/ already holds {len(existing)} input file(s), all of them "
            "for designs in this run. Rewriting them: they are derived from the inventory "
            "and nothing else, so the same inventory produces the same bytes.",
            flush=True,
        )


def write_mono(out_dir, entry):
    """Write one single-record FASTA and return what was recorded about it.

    Returns a dict, or None if the design was skipped upstream.
    """
    design = entry["design"]
    sequence = entry["binder_sequence"]
    path = os.path.join(out_dir, f"{MONO_PREFIX}{design}{FASTA_SUFFIX}")

    with open(path, "w") as handle:
        handle.write(f">{design}\n{sequence}\n")

    return {
        "design": design,
        "path": os.path.basename(path),
        "mode": "mono",
        "record_id": design,
        # Recorded even though Phase 1's slice is the whole array. The value is what node
        # 04 slices by, and Phase 2 changes it; a node 04 that reads it from here cannot
        # inherit 029's head-slice by accident.
        "chain_order": ["binder"],
        "chain_lengths": [len(sequence)],
        "expected_chains": 1,
        "expected_residues": len(sequence),
        "binder_chain_in_design": entry["binder_chain"],
        "designed_chain_source": entry["designed_chain_source"],
        "n_unknown": entry["n_unknown"],
    }


def verify_written(path, expect_colons):
    """Check the FILE, not the string that was meant to be written.

    A guard that re-reads the variable it just assigned asserts the assignment. These
    checks exist because a colon miscount costs GPU minutes and produces a plausible
    number, so they read the bytes back off disk (§8.4).
    """
    with open(path) as handle:
        text = handle.read()

    records = [line for line in text.splitlines() if line.startswith(">")]
    if len(records) != 1:
        fail(
            f"'{os.path.basename(path)}' holds {len(records)} record(s), not 1. One design "
            "per file is the whole defence against 026's conflation bug, where a batch "
            "file made three designs' ranks report as one design's top 5 (818954e9)."
        )

    body = "\n".join(
        line for line in text.splitlines() if line and not line.startswith(">")
    )
    colons = body.count(":")
    if colons != expect_colons:
        if expect_colons == 0:
            fail(
                f"'{os.path.basename(path)}' contains {colons} colon(s) in its sequence. A "
                "colon makes the predictor treat one design as two chains, so this "
                "monomer input would come back as a complex and its confidence would be "
                "read as an interface -- which is 026's bug arriving from the other "
                "direction (d3d707ed)."
            )
        fail(
            f"'{os.path.basename(path)}' contains {colons} colon(s) in its sequence, not "
            f"{expect_colons}. Zero silently predicts a monomer and reports its confidence "
            "as an interface; two predicts three chains."
        )

    if ":" in records[0]:
        fail(
            f"The record id in '{os.path.basename(path)}' contains a colon. 026 reads a "
            "colon in a record id or description as a chain-break marker, which is why "
            "037's UNSCREENED tag carries none (7709809b)."
        )
    return len(records), colons


def main():
    # Phase 1 only. Node 01 refuses this too and refuses it first, which is the point:
    # this node is reachable on its own, so it carries its own copy rather than trusting
    # that something upstream ran.
    if bool_param("PARAM_SCORE_INTERFACE", "score_interface", False):
        fail(
            "score_interface = true is Phase 2 and is not built. This node writes no "
            f"cplx_ input, so the flag would be silently dropped and the run would report "
            "self-consistency numbers under the heading the user asked for (§8.7, §12)."
        )
    numeric_param("PARAM_NUM_MODELS", "num_models", 5, int)

    os.makedirs("./outputs", exist_ok=True)
    inventory = load_inventory("./inputs/design_inventory.json")
    designs = inventory["designs"]
    if not designs:
        fail("node 01's inventory lists no designs. There is nothing to predict.")

    eligible = [e for e in designs if not e.get("skipped")]
    skipped = [e for e in designs if e.get("skipped")]

    # Validate every design before writing any file, so a batch with one bad sequence does
    # not leave a half-written output directory behind.
    for entry in eligible:
        sequence = entry["binder_sequence"]
        bad = sorted({c for c in sequence if c not in ALPHABET})
        if bad:
            fail(
                f"Design '{entry['design']}' has residue code(s) {bad} outside the "
                f"accepted alphabet. 026 validates every sequence against "
                f"{ALPHABET_STRING} and refuses anything else, so this would be "
                "refused inside the predictor after the GPU had already been allocated "
                "(d3d707ed). Node 01 scores an unrecognised residue as X; a code that is "
                "neither a standard residue nor X means something between the two nodes "
                "disagrees about the alphabet."
            )
        if len(sequence) < MIN_BINDER_RESIDUES:
            fail(
                f"Design '{entry['design']}' has {len(sequence)} residue(s), fewer than "
                f"{MIN_BINDER_RESIDUES}. A Kabsch superposition needs three non-collinear "
                "points to be defined and four before the fit stops being exactly "
                "determined by its own inputs, so the self-consistency RMSD this workflow "
                "exists to report would be a number without a meaning."
            )
        if len(entry["binder_ca"]) != len(sequence):
            fail(
                f"Design '{entry['design']}' has {len(sequence)} residue(s) but "
                f"{len(entry['binder_ca'])} coordinate slot(s). Node 01 keeps those arrays "
                "index-aligned on purpose, because node 04's RMSD maps the two structures "
                "by position; a mismatch here means they were not."
            )

    check_rerun("./outputs", [e["design"] for e in eligible])

    records = []
    for entry in eligible:
        record = write_mono("./outputs", entry)
        n_records, n_colons = verify_written(
            os.path.join("./outputs", record["path"]), expect_colons=0
        )
        record["verified_records_in_file"] = n_records
        record["verified_colons_in_sequence"] = n_colons
        records.append(record)

    payload = {
        "n_inputs": len(records),
        "n_skipped": len(skipped),
        "phase": 1,
        "mode": "mono",
        "score_interface": False,
        "alphabet": ALPHABET_STRING,
        "min_binder_residues": MIN_BINDER_RESIDUES,
        # Carried forward from node 01 rather than recomputed, so one wording of the chain
        # provenance reaches the report.
        "designed_chain_sources": inventory.get("designed_chain_sources"),
        "n_target_numbering_gaps": inventory.get("n_target_numbering_gaps"),
        "securedna_screened": False,
        "biosafety_notice": (
            "These sequences are unscreened de novo designs. Passing validation makes "
            "synthesis more likely, not less."
        ),
        "inputs": records,
        "skipped": [
            {"design": e["design"], "reason": e.get("skip_reason")} for e in skipped
        ],
    }
    with open("./outputs/prediction_inputs.json", "w") as handle:
        json.dump(payload, handle, indent=2)

    for record in records:
        print(
            f"  {record['path']}: 1 record, {record['expected_residues']} residue(s), "
            f"{record['verified_colons_in_sequence']} colon(s)"
            + (f", {record['n_unknown']} X" if record["n_unknown"] else ""),
            flush=True,
        )
    for entry in skipped:
        print(f"  SKIPPED {entry['design']}: {entry.get('skip_reason')}", flush=True)

    long_ones = [r for r in records
                 if r["expected_residues"] > LONG_SEQUENCE_WARN_RESIDUES]
    if long_ones:
        print(
            f"WARNING: {len(long_ones)} design(s) exceed {LONG_SEQUENCE_WARN_RESIDUES} "
            "residues. This workflow has measured 35 residues in 61 s and a 63-residue "
            "complex in 80 s on an RTX 3060 Ti, and nothing beyond that, so it cannot say "
            "where this runs out of VRAM -- only that it is the dimension that decides "
            "(§11.5). predict_timeout_s is what stops a runaway, not this line.",
            flush=True,
        )
    if skipped:
        print(
            f"WARNING: {len(skipped)} design(s) were skipped by node 01 and have no "
            "prediction input. They stay in prediction_inputs.json under 'skipped' with "
            "the reason, so node 05 reports them as rows with a stated cause rather than "
            "omitting them (§8.5).",
            flush=True,
        )

    print(
        f"Wrote {len(records)} prediction input(s) and ./outputs/prediction_inputs.json",
        flush=True,
    )


if __name__ == "__main__":
    main()
