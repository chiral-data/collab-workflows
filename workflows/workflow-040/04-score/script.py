#!/usr/bin/env python3
"""Measure each refolded design against the backbone it was designed as.

This is the node the whole workflow exists for, and the node where a wrong answer looks
most like a right one. Everything here is a number, and every number is plausible.

Contract
--------
Reads   design_inventory.json      node 01: binder residue ids and Ca coords, IN ORDER
        design_metrics.json        workflow-037 node 04, passed through by node 01
        *_scores_rank_*.json       node 03, flat and prefixed per design
        *_rank_*.pdb               node 03, all ranks; select rank_001
        prediction_manifest.json   node 03: chain_order, chain_lengths, per-design status
Writes  design_scores.json         one row per design, failures included
        design_scores.csv          the same table

Resolve each design's files by ITS OWN prefix, named in prediction_manifest.json. Never
glob the directory and never pair files by sort order: everything arrives flat, because
Silva's ports match file names rather than paths, and globbing across designs is the
mechanism of 026's conflation bug -- it takes job_id from whichever scores file sorts
first and reports a confident result about a design that does not exist (818954e9).

Metrics
-------
sc_rmsd_A       Ca RMSD of the refolded binder against the designed binder, after Kabsch
                superposition.
plddt_binder    mean pLDDT over the binder residues.
plddt_target    the same for the target. Phase 2. DIAGNOSTIC ONLY -- no threshold, no
                contribution to pass or fail.
iptm            Phase 2.
pae_interaction Phase 2: mean PAE over BOTH off-diagonal blocks, averaged.

Superimpose, and map by position
--------------------------------
Two independent traps here, and the test suite has a case for each.

First: superimpose. A node that computes RMSD without Kabsch returns ~0 for a prediction
that happens to be in the same frame and something enormous for one that is not -- so it
rejects every good design while reporting a batch of plausible-looking failures. That is
exactly how workflow-029's chain-order defect presented: every design failed all four
filters with iPTM near 0.2, binder pLDDT near 40 and Ca RMSD near 16 A, and nothing
errored (d5ba7487). tests/fixtures.json scoring/rmsd-rotated-copy pins it, and pins it
with an answer that is 0 by construction, so the case does not depend on any
superposition code being right. Measured: the same coordinates unsuperposed come out at
259.98 A.

Second: correspondence is by POSITION in the chain, established by node 01 -- never by
residue number. Three renumberings stack: the deposited entry's, RFD4's (it renumbers the
preserved target chain from 1, so barnase in 1BRS, which starts at label 3, arrives
already shifted from what 037 node 02 conditioned on -- ef74646e), and ColabFold's own
renumber-from-1-per-chain. Assert the mapping -- equal lengths, 1:1 -- and refuse a
mismatch. Truncating to the shorter and scoring produces a mis-registered comparison
whose number reads as a mediocre design rather than as a bug.

Slicing pLDDT: read the order, do not assume it
-----------------------------------------------
The scores file carries NO chain boundaries. Measured first-hand, in tests/reference/: the
real monomer file is a flat plddt list[35] with a 35x35 pae, the real complex file a flat
list[63] spanning both chains with a 63x63 pae, and neither names a chain anywhere.

So read chain_order and chain_lengths out of prediction_manifest.json, and ASSERT
len(plddt) == sum(chain_lengths) before slicing. That assertion is the only thing standing
between a position assumption and a plausible wrong number. Workflow-029 slices
plddt[:binder_length] because its binder is chain A; Phase 2's input is TARGET:BINDER --
measured, a 34+29 complex written target-first came back 34 in chain A and 29 in chain B --
so the binder is the TAIL.

Never default a missing metric to a number
------------------------------------------
MEASURED, not hypothetical: in the real monomer run captured in tests/reference/, iptm is
ABSENT -- not zero, not null, absent -- while ptm is present at 0.57. Phase 1 is
monomer-only, so this is the normal case, not an edge case.

Workflow-029 reads scores.get('iptm', 0.0) and sets pae_interaction = float('inf') when
the matrix is absent. The first reports 0.0 -- a measured-looking bad interface -- for
every design in every run; the second fails the threshold, so a missing measurement is
rendered as a verdict. Both are indistinguishable from a real result downstream.

Workflow-026 does it correctly, in one line (03_analyze_confidence/script.py:86):

    iptm = float(data["iptm"]) if "iptm" in data else None

Copy 026. A metric that could not be computed is null, and the reason is recorded in
null_reasons beside it, so "absent because this is a monomer" and "absent because the
file was truncated" do not read alike.

That applies to the rows that give up early too, and at first it did not. score_one()
returns before measuring anything in four cases -- no manifest entry, a prediction node 03
recorded as failed, a prediction recorded 'ok' that named no rank_001 files, and a PDB
that would not parse -- and only the second one filled null_reasons. The other three left
design_scores.json carrying a blank plddt_binder with nothing beside it, which is the
state this whole section exists to forbid, arrived at by a different door. The suite's
no-score-file-is-failure case is what caught it. All four now go through bail(), which
writes the row's own cause into every metric in METRICS_NULLED_ON_BAIL, so a branch added
later cannot forget to.

Merge by name
-------------
Keyed on the design name, never by position. Workflow-037's metrics rows may arrive in a
different order from the predictions, and merging by position puts one design's generation
metrics against another's scores with every row still plausible. 026's conflation bug is
the same error one level up: it takes job_id from whichever scores file sorts first
(818954e9).

037's rows are keyed by the CIF basename (design_0_generated.cif) and this workflow keys
on the stem (design_0), so node 01 recorded source_file for exactly this join. Matching on
the stem by string surgery here would be a second place that knows how 037 names things.

Order and thresholds are decided HERE, once
-------------------------------------------
Node 04 sorts and node 05 inherits, so there is one source of order. Sort NUMERICALLY:
lexicographic puts design_10 between design_1 and design_2 (#51). Likewise pass/fail is
decided here and carried, not recomputed in node 05 from the thresholds -- two nodes
applying the same rule eventually disagree about the same design, which is the defect
179003b9 exists to remove.

Record each threshold beside the value it was compared against, so a reader can re-filter
without re-running. Every threshold is a default to calibrate, not truth: 037's
compactness limit was first drafted at 1.5 and would have fired on three real monomer
designs at 1.44, 1.46 and 1.49, and running it over the archive before shipping is what
caught that (fd75d77b).

Refuse non-finite thresholds by name
------------------------------------
037's numeric_param() accepted nan and inf. Neither is caught by a range check, because
every comparison against nan is false: nan passed "must be positive" and passed "outside
the expected range" without a warning, then made every distance test false, so the node
reported no contacts on a structure that has plenty (55edd8fc). All five numeric
parameters here are compared against measurements, and all of them must fail the same way
-- 037 created that inconsistency itself by adding two clean parameters and leaving two
raw int() casts ten lines above them (975ab379).

Fails when
----------
  * len(plddt) does not match the recorded chain lengths
  * the residue mapping is not 1:1 and equal-length
  * the prediction's chain count does not match what node 02 declared
  * a threshold parameter is nan, inf, non-numeric or non-positive
  * score_interface is true (Phase 2 is not built)
  * designed_chain_source is single-chain and interface scoring was requested
  * every design failed, so there is nothing to report
"""
import csv
import json
import math
import os
import re
import sys

COLUMNS = [
    "design",
    "status",
    "reason",
    "binder_chain",
    "designed_chain_source",
    "binder_length",
    # Carried, not looked up. Node 05 exports the passing sequences and node 04 is its
    # only dependency, so the alternative is node 05 declaring design_inventory.json as an
    # input it would otherwise have no business reading -- and silva would not route it
    # there anyway, since node 05 depends on this node alone. Passing the fact forward is
    # the rule (CLAUDE.md).
    "binder_sequence",
    "sc_rmsd_A",
    "sc_rmsd_n_residues",
    "plddt_binder",
    "plddt_target",
    "ptm",
    "iptm",
    "pae_interaction",
    "max_pae",
    "n_unknown",
    "passes",
    "failed_thresholds",
    "max_sc_rmsd_a",
    "min_plddt_binder",
    # Phase 2. Null on a Phase 1 run, where they are not merely unmeasured but
    # meaningless: there is no interface in a monomer prediction.
    #
    # plddt_binder_in_complex is NOT the same number as plddt_binder and does not replace
    # it. plddt_binder is the binder predicted ALONE, which is the self-consistency
    # question; this is the binder predicted in context. 029 had only one of the two and
    # so could not tell a binder that folds alone but not against its target from one
    # that folds in neither case.
    "plddt_binder_in_complex",
    "interface_scored",
    "interface_gated",
    "interface_reason",
    "min_iptm",
    "max_pae_interaction",
    "cplx_msa_mode",
    "cplx_binder_aligned",
    # Merged from 037, by name. Kept to the few that say whether a bad score has a
    # geometric explanation already on record; the full table stays in design_metrics.json.
    "gen_compactness_ratio",
    "gen_overextended",
    "gen_chain_breaks",
    "gen_steric_clash",
    "gen_n_interface_contacts",
]

# Below this, a target has not been modelled well enough for any interface number computed
# against it to mean anything -- which is 040's single most likely cause of bad Phase 2
# numbers, given §4. It warns and never gates: 029's absence of this number produced a
# wrong diagnosis, and the fix is to report it, not to filter on it (2359ffa4).
LOW_TARGET_PLDDT_WARN = 70.0

# Every metric that score_one() can report, so a row that gives up early can say why each
# one is null in one place rather than per branch. Listed here, beside COLUMNS, because a
# metric added to COLUMNS and forgotten here comes out null and unexplained.
METRICS_NULLED_ON_BAIL = (
    "sc_rmsd_A", "plddt_binder", "plddt_target", "iptm", "pae_interaction",
    "ptm", "max_pae", "plddt_binder_in_complex",
)

# Phase 2's own metrics. The comment above COLUMNS says they are "null on a Phase 1 run,
# where they are not merely unmeasured but meaningless" -- and that sentence has to be in
# the ROW, not only in the source. plddt_binder_in_complex was added to COLUMNS and
# forgotten here, and came out exactly as the comment beside METRICS_NULLED_ON_BAIL
# predicted: null and unexplained, in the JSON, while the CSV omitted the column
# altogether. Two representations of one fact disagreeing is how a reader concludes the
# metric was measured and came back empty.
INTERFACE_METRICS = ("plddt_binder_in_complex", "iptm", "pae_interaction", "plddt_target")


def explain_unset_interface(row, null_reasons, score_interface):
    """Every interface metric still None must say why. Nulls are explained or they lie."""
    for metric in INTERFACE_METRICS:
        if row.get(metric) is not None or null_reasons.get(metric):
            continue
        null_reasons[metric] = (
            "score_interface is false, so this run predicted the binder alone. There is "
            "no interface in a monomer prediction, so this is meaningless here rather "
            "than merely unmeasured (§4)."
            if not score_interface else
            "no complex prediction was scored for this design, and nothing above said "
            "why. This reason is the fallback and its presence is itself a defect: it "
            "means a branch returned without explaining its own null (§8.7)."
        )


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def numeric_param(env_name, display_name, default, cast):
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


def positive_param(env_name, display_name, default, cast):
    value = numeric_param(env_name, display_name, default, cast)
    if value <= 0:
        fail(f"{display_name} must be positive, got {value}.")
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


def kabsch_rmsd(mobile, reference):
    """Ca RMSD after optimal superposition, and the number of positions it is over.

    The superposition is the measurement. Without it this function returns the distance
    between two structures in whatever frames they happen to be in, which for a real
    prediction is hundreds of Angstroms and reads as a batch of failed designs rather
    than as a missing step (d5ba7487).
    """
    import numpy as np

    P = np.asarray(mobile, dtype=float)
    Q = np.asarray(reference, dtype=float)
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    H = Pc.T @ Qc
    U, _, Vt = np.linalg.svd(H)
    # Without this reflection correction the "best fit" can be a mirror image, which fits
    # a mirrored structure perfectly and reports it as a recovered fold.
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    rotation = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    aligned = (rotation @ Pc.T).T
    return float(np.sqrt(((aligned - Qc) ** 2).sum() / len(P))), len(P)


def read_pdb_chains(path):
    """[(chain_id, [ca_xyz, ...]), ...] in file order, from a predicted PDB.

    A different parser from node 01's mmCIF reader on purpose: that file came from the
    generator and this one comes from the predictor, and sharing a parser between the
    thing being measured and the thing it is measured against is how a check passes for
    the wrong reason (§8.4).
    """
    chains, order = {}, []
    try:
        handle = open(path)
    except OSError as exc:
        return None, f"could not open '{os.path.basename(path)}': {exc}"
    with handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[12:16].strip() != "CA":
                continue
            chain = line[21]
            if chain not in chains:
                chains[chain] = []
                order.append(chain)
            try:
                chains[chain].append(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                )
            except ValueError as exc:
                return None, f"malformed coordinate in '{os.path.basename(path)}': {exc}"
    if not order:
        return None, f"'{os.path.basename(path)}' has no CA atoms"
    return [(cid, chains[cid]) for cid in order], None


def slice_by_order(values, chain_order, chain_lengths, what):
    """{role: values_for_that_role}, sliced by POSITION in the recorded order.

    Refuses rather than slicing when the lengths do not add up. That refusal is the only
    thing standing between a position assumption and a plausible wrong number (§8.1).
    """
    total = sum(chain_lengths)
    if len(values) != total:
        fail(
            f"{what} has {len(values)} entr(ies) but node 02 recorded chain order "
            f"{chain_order} with lengths {chain_lengths}, summing to {total}. The scores "
            "file carries no chain boundaries at all, so this array is sliced by position "
            "and nothing else; a length that does not match means the slice would land "
            "somewhere other than where it is recorded to land, and every per-chain "
            "number after it would be plausible and wrong (§8.1)."
        )
    out, at = {}, 0
    for role, length in zip(chain_order, chain_lengths):
        out[role] = values[at:at + length]
        at += length
    return out


def interface_pae(matrix, chain_order, chain_lengths):
    """Mean PAE over BOTH off-diagonal blocks, averaged. Phase 2.

    PAE is not symmetric: entry (i, j) is the predicted error in residue j's position when
    the frame is aligned on residue i, so target->binder and binder->target are different
    measurements. Averaging one block alone reports whichever direction the code happened
    to index first. The fixture sets the two blocks to 6.0 and 14.0 so that the right
    answer is exactly 10.0 and either single-block reading is visibly wrong.
    """
    first = chain_lengths[0]
    total = sum(chain_lengths)
    if len(matrix) != total or any(len(row) != total for row in matrix):
        return None, (
            f"the pae matrix is {len(matrix)}x"
            f"{len(matrix[0]) if matrix else 0}, not {total}x{total}"
        )
    upper = [matrix[i][j] for i in range(first) for j in range(first, total)]
    lower = [matrix[i][j] for i in range(first, total) for j in range(first)]
    if not upper or not lower:
        return None, "one of the off-diagonal blocks is empty"
    means = (sum(upper) / len(upper), sum(lower) / len(lower))
    return (means[0] + means[1]) / 2.0, None


def natural_key(name):
    """Sort runs of digits as numbers. design_10 does not belong between 1 and 2 (#51)."""
    return [int(p) if p.isdigit() else p for p in re.split(r"(\d+)", name)]


def load_json(path, what, required=True):
    if not os.path.exists(path):
        if required:
            fail(f"No '{os.path.basename(path)}' in ./inputs/ ({what}).")
        return None
    try:
        with open(path) as handle:
            return json.load(handle)
    except Exception as exc:
        fail(f"Could not read '{path}' ({what}): {exc}")


def score_one(entry, prediction, max_rmsd, min_plddt):
    """One row. Never raises for a per-design problem: it records it (§8.5)."""
    design = entry["design"]
    row = {column: None for column in COLUMNS}
    row.update({
        "design": design,
        "binder_chain": entry["binder_chain"],
        "designed_chain_source": entry["designed_chain_source"],
        "binder_length": entry["n_binder_residues"],
        "binder_sequence": entry["binder_sequence"],
        "n_unknown": entry["n_unknown"],
        "max_sc_rmsd_a": max_rmsd,
        "min_plddt_binder": min_plddt,
        "passes": False,
        "failed_thresholds": "",
        "status": "failed",
    })
    null_reasons = {}

    def bail(reason):
        """Record the row's cause AND the per-metric reason, for every metric at once.

        Every metric on this row is null, and a null with no recorded reason cannot be
        told apart from a null nobody noticed (§8.7) -- which is the whole argument for
        null_reasons existing. Three of these early returns used to set row["reason"] and
        leave null_reasons empty, so design_scores.json carried a blank plddt_binder with
        nothing beside it; the suite's no-score-file-is-failure case caught it. One
        function now, so a branch added later cannot forget.
        """
        row["reason"] = reason
        for metric in METRICS_NULLED_ON_BAIL:
            null_reasons[metric] = reason
        return row, null_reasons

    if prediction is None:
        return bail(
            "no prediction: node 03's manifest has no entry for this design. It is a row "
            "here with a stated cause rather than an omission (§8.5)."
        )

    if prediction.get("status") != "ok":
        return bail(prediction.get("reason") or (
            f"node 03 recorded status '{prediction.get('status')}', so nothing was "
            "measured"
        ))

    chain_order = prediction["chain_order"]
    chain_lengths = prediction["chain_lengths"]

    scores_name = prediction.get("scores_rank_001")
    pdb_name = prediction.get("pdb_rank_001")
    if not scores_name or not pdb_name:
        return bail(
            "node 03 recorded status 'ok' but named no rank_001 "
            + ("scores file" if not scores_name else "PDB")
            + ". Exit codes see nothing (§8.3): there is nothing here to measure."
        )

    scores = load_json(os.path.join("./inputs", scores_name), "prediction scores")

    chains, error = read_pdb_chains(os.path.join("./inputs", pdb_name))
    if chains is None:
        return bail(f"could not read the predicted structure: {error}")

    if len(chains) != len(chain_order):
        # Fatal for the run, not for this design. A complex that came back one chain means
        # the ':' was dropped in the inputs, and the next design has the same inputs: iptm
        # and pae_interaction computed from it would be one chain's confidence reported as
        # an interface, which is 026's exact failure.
        fail(
            f"'{design}' was predicted as {len(chains)} chain(s) "
            f"{[cid for cid, _ in chains]}, but node 02 recorded chain order "
            f"{chain_order} -- {len(chain_order)} chain(s). For a complex this means the "
            "':' was dropped somewhere and the result is a monomer masquerading as one, "
            "so every interface number computed from it would describe a single chain. "
            "Refusing the whole run rather than this design, because the cause is "
            "upstream of it."
        )

    predicted = {}
    for role, (chain_id, coords) in zip(chain_order, chains):
        predicted[role] = coords
    for role, length in zip(chain_order, chain_lengths):
        if len(predicted[role]) != length:
            fail(
                f"'{design}': the predicted {role} chain has {len(predicted[role])} "
                f"residue(s) but node 02 recorded {length}, so the mapping between the "
                "declaration and the structure is broken. Correspondence here is by "
                "POSITION and nothing else, because three renumberings stack in this "
                "workflow and no correspondence by residue number survives them "
                "(ef74646e). An unequal mapping cannot be scored: truncating to the "
                "shorter one would report an RMSD over a different set of residues than "
                "the one it names (§8.1)."
            )

    # And the other end of the same mapping, which is a DIFFERENT check: above compares
    # the structure against node 02's declaration, this compares it against the design it
    # is supposed to be a refolding of. A prediction can agree with the declaration and
    # still not match the design -- that is what happens when node 02 and node 01 disagree
    # about the binder -- and nothing else in the pipeline looks at both.
    #
    # Found by the suite: with the prediction two residues longer than the design, every
    # check above passed and zip() silently truncated to the shorter array. The RMSD came
    # out 0.0 and the design was reported as passing. This is the exact shape of the
    # defect the case exists for, and it was present until the case was run.
    designed_length = len(entry["binder_ca"])
    if len(predicted["binder"]) != designed_length:
        fail(
            f"'{design}': the predicted binder has {len(predicted['binder'])} residue(s) "
            f"and the designed binder has {designed_length}, so there is no 1:1 mapping "
            "between them. Refusing rather than scoring the overlap: zip() would truncate "
            "to the shorter array and report an RMSD over a prefix, which comes out small "
            "and reads as a design that recovered its fold (§8.1)."
        )

    # Designed Ca, in chain order, dropping only positions the generator left without a
    # Ca -- and dropping the SAME positions from the prediction, so the 1:1 correspondence
    # survives. Node 01 holds null in place rather than compacting the array for exactly
    # this reason.
    pairs = [
        (designed, predicted_xyz)
        for designed, predicted_xyz in zip(entry["binder_ca"], predicted["binder"])
        if designed is not None
    ]
    n_missing = len(entry["binder_ca"]) - len(pairs)
    if len(pairs) < 4:
        row["reason"] = (
            f"only {len(pairs)} residue(s) of '{design}' have a Ca in the designed "
            "structure, which is too few for a superposition to be defined by anything "
            "other than its own inputs."
        )
        null_reasons["sc_rmsd_A"] = row["reason"]
    else:
        rmsd, n = kabsch_rmsd(
            [p for _, p in pairs], [d for d, _ in pairs]
        )
        # Four decimal places, which is one more than the data carries: a PDB ATOM record
        # has three, so 0.001 A is the finest distinction the input can express and a
        # reported RMSD must not be coarser than that. It matters at the bottom of the
        # range, which is where the keystone case lives -- a rotated copy comes back at
        # 0.00051 A from PDB rounding alone, and rounding to three decimals would print
        # that as 0.001, sitting exactly on the fixture's tolerance. Node 05 formats for
        # a reader; this file is the data.
        row["sc_rmsd_A"] = round(rmsd, 4)
        row["sc_rmsd_n_residues"] = n
        if n_missing:
            print(
                f"NOTE: '{design}': {n_missing} designed residue(s) have no Ca, so the "
                f"RMSD is over {n} of {entry['n_binder_residues']} positions. The same "
                "positions were dropped from both structures, so the correspondence is "
                "still 1:1.",
                flush=True,
            )

    # pLDDT. Absent is null with a reason, never a number (§8.7).
    if "plddt" not in scores:
        null_reasons["plddt_binder"] = (
            f"'{scores_name}' carries no 'plddt' key. 029 would report 0.0 here, which "
            "reads as a measured, terrible fold (§8.7)."
        )
    else:
        blocks = slice_by_order(
            scores["plddt"], chain_order, chain_lengths, f"'{scores_name}' plddt"
        )
        binder = blocks["binder"]
        row["plddt_binder"] = round(sum(binder) / len(binder), 2)
        if "target" in blocks and blocks["target"]:
            target = blocks["target"]
            row["plddt_target"] = round(sum(target) / len(target), 2)
        else:
            null_reasons["plddt_target"] = (
                "this is a monomer prediction: there is no target chain in it, so there "
                "is no target confidence to report. Reporting 0.0 would read as a target "
                "that modelled terribly."
            )

    row["ptm"] = float(scores["ptm"]) if "ptm" in scores else None
    if "ptm" not in scores:
        null_reasons["ptm"] = f"'{scores_name}' carries no 'ptm' key"
    row["max_pae"] = float(scores["max_pae"]) if "max_pae" in scores else None
    if "max_pae" not in scores:
        null_reasons["max_pae"] = f"'{scores_name}' carries no 'max_pae' key"

    # 026's line, exactly (03_analyze_confidence/script.py:86). MEASURED: the real monomer
    # scores file has no iptm key at all, so this is Phase 1's normal case.
    row["iptm"] = float(scores["iptm"]) if "iptm" in scores else None
    if "iptm" not in scores:
        null_reasons["iptm"] = (
            "absent from the scores file, which is what a monomer prediction produces: "
            "iptm is an interface quantity and alphafold2_ptm does not emit it. Measured "
            "in tests/reference/. 029's scores.get('iptm', 0.0) would report a "
            "measured-looking 0.0 for every design in every monomer run (§8.7)."
        )

    if "target" not in chain_order:
        null_reasons["pae_interaction"] = (
            "there is no interface in a monomer prediction, so there are no off-diagonal "
            "blocks to average. 029 sets inf here when the matrix is absent, which fails "
            "the threshold and so renders a missing measurement as a verdict (§8.7)."
        )
    elif "pae" not in scores:
        null_reasons["pae_interaction"] = (
            f"'{scores_name}' carries no 'pae' key, so the interface error was not "
            "measured. Not inf: inf fails the threshold and would be read as a verdict."
        )
    else:
        value, why = interface_pae(scores["pae"], chain_order, chain_lengths)
        if value is None:
            null_reasons["pae_interaction"] = why
        else:
            row["pae_interaction"] = round(value, 3)

    # The verdict, decided here and carried. Node 05 renders it rather than recomputing
    # it: two nodes applying the same rule eventually disagree about the same design
    # (179003b9).
    #
    # plddt_target is deliberately absent from this list. It is diagnostic, by decision:
    # 029's wrong diagnosis came from not having the number, and the fix is to report it,
    # not to filter on it (2359ffa4).
    failed = []
    if row["sc_rmsd_A"] is None:
        failed.append("sc_rmsd_A (not measured)")
    elif row["sc_rmsd_A"] > max_rmsd:
        failed.append(f"sc_rmsd_A {row['sc_rmsd_A']} > {max_rmsd}")
    if row["plddt_binder"] is None:
        failed.append("plddt_binder (not measured)")
    elif row["plddt_binder"] < min_plddt:
        failed.append(f"plddt_binder {row['plddt_binder']} < {min_plddt}")

    row["status"] = "scored"
    row["passes"] = not failed
    row["failed_thresholds"] = "; ".join(failed)
    row["reason"] = None if not failed else "; ".join(failed)
    return row, null_reasons


def add_interface(entry, row, null_reasons, prediction, limits, gate):
    """Overlay the interface metrics from the COMPLEX prediction onto a scored row.

    Phase 2 ADDS interface scoring; it does not replace self-consistency. sc_rmsd_A and
    plddt_binder keep coming from the monomer prediction, because whether the binder folds
    as designed is the question this workflow exists to answer and an interface metric
    says nothing about it. This function only ever writes plddt_target, iptm and
    pae_interaction -- the three things a monomer cannot produce.

    GATING IS OPT-IN AND OFF BY DEFAULT, which is a decision and not an omission:

      * §7 ships pae_interaction as a provisional band rather than a gate. Its published
        cutoff comes from protocols using AF2 INITIAL GUESS, which hands the predictor the
        designed complex as a starting point; vanilla ColabFold does not, so the number
        does not transfer. It is also floored by how well the TARGET was modelled -- a
        target at ~42 pLDDT moved the floor from 7.6 to 15.4 A and put 10.0 out of reach
        for a whole batch (8e4772a7, §4).
      * And the Tier 2 controls measured what an uncalibrated threshold on a confidence
        metric is worth: 115 alanines score 95.81 plddt_binder, above the real design, and
        clear the 80.0 cutoff. A metric can be confidently excellent about nothing
        (§8.6). iptm and pae_interaction have had no equivalent control run against them
        at all.

    So the numbers are reported, with their thresholds beside them and a recorded flag
    saying they did not gate. A reader must not be able to mistake `passes` for a verdict
    about the interface, which is why the flag is in the row rather than in the README.
    """
    row["interface_scored"] = False
    row["interface_gated"] = bool(gate)
    row["min_iptm"] = limits["min_iptm"]
    row["max_pae_interaction"] = limits["max_pae_interaction"]

    if prediction is None:
        reason = (
            "score_interface was asked for, but node 03's manifest has no COMPLEX "
            "prediction for this design. The self-consistency numbers above are real; "
            "the interface was not measured, and is null rather than a number (§8.7)."
        )
        for metric in ("iptm", "pae_interaction", "plddt_target"):
            null_reasons[metric] = reason
            row[metric] = None
        row["interface_reason"] = reason
        return row, null_reasons

    if prediction.get("status") != "ok":
        reason = prediction.get("reason") or (
            f"the complex prediction recorded status {prediction.get('status')!r}"
        )
        for metric in ("iptm", "pae_interaction", "plddt_target"):
            null_reasons[metric] = reason
            row[metric] = None
        row["interface_reason"] = reason
        return row, null_reasons

    scores_name = prediction.get("scores_rank_001")
    if not scores_name:
        reason = (
            "the complex prediction is 'ok' but names no rank_001 scores file, so there "
            "is nothing to read (§8.3)."
        )
        for metric in ("iptm", "pae_interaction", "plddt_target"):
            null_reasons[metric] = reason
            row[metric] = None
        row["interface_reason"] = reason
        return row, null_reasons

    scores = load_json(os.path.join("./inputs", scores_name), "complex scores")
    chain_order = prediction["chain_order"]
    chain_lengths = prediction["chain_lengths"]

    # The assertion §8.1 asks for, on the complex array this time. If the recorded lengths
    # do not describe the array, every slice below is of the wrong residues and each one
    # still produces a plausible number.
    if "plddt" not in scores:
        # Both columns come out of this one array, so both are absent when it is, and both
        # say so. Only plddt_target did at first, which left the other reported as simply
        # not there (§8.7).
        for metric in ("plddt_target", "plddt_binder_in_complex"):
            null_reasons[metric] = f"'{scores_name}' carries no 'plddt' key"
    else:
        blocks = slice_by_order(
            scores["plddt"], chain_order, chain_lengths, f"'{scores_name}' plddt"
        )
        target = blocks.get("target") or []
        if target:
            row["plddt_target"] = round(sum(target) / len(target), 2)
            null_reasons.pop("plddt_target", None)
        else:
            null_reasons["plddt_target"] = (
                f"'{scores_name}' has no target block under chain order {chain_order}"
            )
        # The binder's confidence IN THE COMPLEX, which is a different number from the
        # monomer's and is recorded separately rather than overwriting it. 029 had only
        # one and could not tell a binder that folds alone but not in context from one
        # that does neither.
        binder = blocks.get("binder") or []
        if binder:
            row["plddt_binder_in_complex"] = round(sum(binder) / len(binder), 2)
        else:
            # Said, rather than left absent. The target block above records a reason when
            # it is empty and this one did not, so a complex whose binder block could not
            # be sliced produced a row with the column simply missing and nothing saying
            # why -- §8.3's case and §8.7's. An absent metric is an absence with a reason
            # attached, never a gap a reader is left to fill in.
            null_reasons["plddt_binder_in_complex"] = (
                f"'{scores_name}' has no binder block under chain order {chain_order} "
                f"with lengths {chain_lengths}"
            )

    row["iptm"] = float(scores["iptm"]) if "iptm" in scores else None
    if "iptm" in scores:
        null_reasons.pop("iptm", None)
    else:
        null_reasons["iptm"] = (
            f"'{scores_name}' carries no 'iptm' key, although this is a complex "
            f"prediction with model_type {prediction.get('model_type')!r}. Measured in "
            "tests/reference/: a real multimer scores file does carry it, so its absence "
            "here means the prediction was not the multimer it was recorded as."
        )

    if "pae" not in scores:
        null_reasons["pae_interaction"] = (
            f"'{scores_name}' carries no 'pae' key, so the interface error was not "
            "measured. Not inf: inf fails a threshold and renders a missing measurement "
            "as a verdict (§8.7)."
        )
    else:
        value, why = interface_pae(scores["pae"], chain_order, chain_lengths)
        if value is None:
            null_reasons["pae_interaction"] = why
        else:
            row["pae_interaction"] = round(value, 3)
            null_reasons.pop("pae_interaction", None)

    row["interface_scored"] = True
    row["interface_reason"] = None
    row["cplx_msa_mode"] = prediction.get("msa_mode")
    row["cplx_binder_aligned"] = prediction.get("binder_aligned", False)

    # The low-target warning, which never gates. 029's wrong diagnosis came from not
    # HAVING this number, and the fix is to report it rather than filter on it (2359ffa4).
    if row.get("plddt_target") is not None and row["plddt_target"] < LOW_TARGET_PLDDT_WARN:
        print(
            f"WARNING: '{row['design']}': the target modelled at "
            f"{row['plddt_target']} pLDDT, below {LOW_TARGET_PLDDT_WARN}. Interface PAE "
            "can be no better than the worse-modelled of the two chains (8e4772a7), so "
            f"this design's pae_interaction of {row.get('pae_interaction')} is floored by "
            "the target rather than by the design. This is diagnostic and does not gate.",
            flush=True,
        )

    if not gate:
        return row, null_reasons

    # Opt-in gating. Added to the existing failure list rather than replacing it, so a
    # design that fails self-consistency still says so.
    failed = [part for part in (row.get("failed_thresholds") or "").split("; ") if part]
    if row["iptm"] is None:
        failed.append("iptm (not measured)")
    elif row["iptm"] < limits["min_iptm"]:
        failed.append(f"iptm {row['iptm']} < {limits['min_iptm']}")
    if row["pae_interaction"] is None:
        failed.append("pae_interaction (not measured)")
    elif row["pae_interaction"] > limits["max_pae_interaction"]:
        failed.append(
            f"pae_interaction {row['pae_interaction']} > "
            f"{limits['max_pae_interaction']}"
        )
    row["passes"] = not failed
    row["failed_thresholds"] = "; ".join(failed)
    row["reason"] = None if not failed else "; ".join(failed)
    return row, null_reasons


def main():
    score_interface = bool_param("PARAM_SCORE_INTERFACE", "score_interface", False)

    max_rmsd = positive_param("PARAM_MAX_SC_RMSD_A", "max_sc_rmsd_a", 2.0, float)
    min_plddt = numeric_param("PARAM_MIN_PLDDT_BINDER", "min_plddt_binder", 80.0, float)
    if not 0.0 <= min_plddt <= 100.0:
        fail(f"min_plddt_binder is a pLDDT, so it must be 0-100, got {min_plddt}.")

    # Phase 2 limits. Read and range-checked on every run, including a Phase 1 one: a
    # parameter that is only validated when it is used is a parameter that accepts nan
    # quietly until the day it matters (§8.8).
    min_iptm = numeric_param("PARAM_MIN_IPTM", "min_iptm", 0.5, float)
    if not 0.0 <= min_iptm <= 1.0:
        fail(f"min_iptm is an iPTM, so it must be 0-1, got {min_iptm}.")
    max_pae_interaction = positive_param(
        "PARAM_MAX_PAE_INTERACTION", "max_pae_interaction", 10.0, float
    )
    gate_interface = bool_param("PARAM_GATE_ON_INTERFACE", "gate_on_interface", False)
    if gate_interface and not score_interface:
        fail(
            "gate_on_interface = true with score_interface = false. There would be no "
            "interface metric to gate on, so every design would fail for want of a "
            "measurement that was never asked for. Set score_interface = true, or leave "
            "both off."
        )

    os.makedirs("./outputs", exist_ok=True)

    inventory = load_json("./inputs/design_inventory.json", "node 01's inventory")
    manifest = load_json("./inputs/prediction_manifest.json", "node 03's manifest")
    metrics = load_json(
        "./inputs/design_metrics.json", "workflow-037's generation metrics", required=False
    )

    designs = inventory.get("designs") or []
    if not designs:
        fail("node 01's inventory lists no designs.")

    # Keyed by the design name, never by position (818954e9) -- AND by mode, because
    # Phase 2 puts TWO predictions in the manifest for every design. A dict keyed on the
    # design alone would silently keep whichever came last, which is 026's conflation bug
    # arriving by a new route: the complex's confidence would be read as the monomer's
    # self-consistency, and every number would still look reasonable.
    predictions, collisions = {}, []
    for prediction in manifest.get("predictions") or []:
        key = (prediction["design"], prediction.get("mode", "mono"))
        if key in predictions:
            collisions.append(key)
        predictions[key] = prediction
    if collisions:
        fail(
            f"node 03's manifest holds more than one prediction for {collisions}. One "
            "design and one mode name exactly one prediction; two means the manifest was "
            "written twice or a design name repeats, and keying by it would silently take "
            "whichever came last (818954e9)."
        )

    limits = {"min_iptm": min_iptm, "max_pae_interaction": max_pae_interaction}
    if score_interface:
        n_cplx = len([k for k in predictions if k[1] == "cplx"])
        if not n_cplx:
            fail(
                "score_interface = true, but node 03's manifest holds no complex "
                "prediction at all. Scoring the monomers and reporting them under a "
                "Phase 2 heading is the failure this check exists for: every interface "
                "column would be null and the run would otherwise look successful "
                "(§8.3, §8.7)."
            )
        print(
            f"Phase 2: {n_cplx} complex prediction(s) to score alongside the "
            f"self-consistency ones. Interface thresholds min_iptm={min_iptm} "
            f"max_pae_interaction={max_pae_interaction} are "
            + ("GATING" if gate_interface else "REPORTED AND DO NOT GATE")
            + ".",
            flush=True,
        )
        if not gate_interface:
            print(
                "  They do not gate by default on purpose. §7 ships pae_interaction as a "
                "provisional band: its published cutoff assumes AF2 initial guess, which "
                "this workflow does not use, and it is floored by how well the target "
                "modelled (8e4772a7). And the Tier 2 controls measured what an "
                "uncalibrated confidence threshold is worth -- 115 alanines clear "
                "plddt_binder 80.0 at 95.81 (§8.6). Set PARAM_GATE_ON_INTERFACE=true to "
                "gate anyway; every row records which it was.",
                flush=True,
            )
    if manifest.get("dry_run"):
        fail(
            "prediction_manifest.json records dry_run = true, so node 03 resolved its "
            "argv and never ran the predictor. There is nothing to score, and scoring "
            "the absence would produce a report full of nulls that looked like a result."
        )

    # 037 keys its rows by the CIF basename; node 01 recorded source_file so this join
    # does not have to know how 037 names things.
    gen_rows = {}
    if metrics:
        by_file = {
            os.path.basename(r.get("design", "")): r for r in metrics.get("designs") or []
        }
        for entry in designs:
            row = by_file.get(entry["source_file"])
            if row is not None:
                gen_rows[entry["design"]] = row
        missing = [e["design"] for e in designs if e["design"] not in gen_rows]
        if missing:
            print(
                f"WARNING: {len(missing)} design(s) have no row in design_metrics.json: "
                f"{', '.join(missing)}. Their generation metrics are blank in the report. "
                "They are matched by node 01's recorded source_file, never by position: "
                "merging by position puts one design's geometry against another's scores "
                "with every row still plausible (818954e9).",
                flush=True,
            )

    rows, reasons = [], {}
    for entry in sorted(designs, key=lambda e: natural_key(e["design"])):
        if entry.get("skipped"):
            row = {column: None for column in COLUMNS}
            row.update({
                "design": entry["design"],
                "status": "skipped",
                "reason": entry.get("skip_reason"),
                "binder_chain": entry["binder_chain"],
                "designed_chain_source": entry["designed_chain_source"],
                "binder_length": entry["n_binder_residues"],
                "binder_sequence": entry["binder_sequence"],
                "n_unknown": entry["n_unknown"],
                "passes": False,
                "failed_thresholds": "",
                "max_sc_rmsd_a": max_rmsd,
                "min_plddt_binder": min_plddt,
            })
            rows.append(row)
            reasons[entry["design"]] = {
                "sc_rmsd_A": entry.get("skip_reason"),
                "plddt_binder": entry.get("skip_reason"),
            }
            continue

        row, null_reasons = score_one(
            entry, predictions.get((entry["design"], "mono")), max_rmsd, min_plddt
        )
        if score_interface and row["status"] == "scored":
            row, null_reasons = add_interface(
                entry, row, null_reasons,
                predictions.get((entry["design"], "cplx")), limits, gate_interface,
            )
        explain_unset_interface(row, null_reasons, score_interface)
        gen = gen_rows.get(entry["design"])
        if gen:
            row["gen_compactness_ratio"] = gen.get("compactness_ratio")
            row["gen_overextended"] = gen.get("overextended")
            row["gen_chain_breaks"] = gen.get("chain_breaks")
            row["gen_steric_clash"] = gen.get("steric_clash")
            row["gen_n_interface_contacts"] = gen.get("n_interface_contacts")
        rows.append(row)
        reasons[entry["design"]] = null_reasons

    scored = [r for r in rows if r["status"] == "scored"]
    passing = [r for r in scored if r["passes"]]

    payload = {
        "n_designs": len(rows),
        "n_scored": len(scored),
        "n_passing": len(passing),
        "n_failed": len([r for r in rows if r["status"] == "failed"]),
        "n_skipped": len([r for r in rows if r["status"] == "skipped"]),
        "phase": 2 if score_interface else 1,
        "score_interface": score_interface,
        "interface_gated": bool(gate_interface),
        "n_interface_scored": len([r for r in rows if r.get("interface_scored")]),
        # Beside the values they were compared against, so a reader can re-filter without
        # re-running. Defaults to calibrate, not truth (§7).
        "thresholds": {
            "max_sc_rmsd_a": max_rmsd,
            "min_plddt_binder": min_plddt,
            **({"min_iptm": min_iptm,
                "max_pae_interaction": max_pae_interaction} if score_interface else {}),
        },
        "thresholds_are_calibrated": False,
        # What actually decided `passes`, as opposed to what was reported. On a Phase 2
        # run with gating off, iptm and pae_interaction are in the row and not in this
        # list, and that difference is the whole point of recording both.
        "gating_metrics": ["sc_rmsd_A", "plddt_binder"] + (
            ["iptm", "pae_interaction"] if (score_interface and gate_interface) else []
        ),
        "diagnostic_metrics": ["plddt_target", "ptm", "max_pae"] + (
            ["plddt_binder_in_complex"] if score_interface else []
        ) + ([] if (score_interface and gate_interface) else
             ["iptm", "pae_interaction"]),
        "why_interface_does_not_gate": None if (
            not score_interface or gate_interface
        ) else (
            "§7 ships pae_interaction as a provisional band rather than a gate: its "
            "published cutoff comes from protocols using AF2 initial guess, which this "
            "workflow does not use, and it is floored by how well the TARGET modelled -- "
            "a target at ~42 pLDDT moved the floor from 7.6 to 15.4 A and put 10.0 out of "
            "reach for a whole batch (8e4772a7). The Tier 2 controls also measured what "
            "an uncalibrated confidence threshold is worth: 115 alanines score 95.81 "
            "plddt_binder, above the real design, and clear the 80.0 cutoff (§8.6). iptm "
            "and pae_interaction have had no equivalent control run against them. They "
            "are reported with their thresholds beside them; set "
            "PARAM_GATE_ON_INTERFACE=true to gate anyway."
        ),
        "order": "numeric by design name, decided here; node 05 inherits it",
        "designed_chain_sources": inventory.get("designed_chain_sources"),
        "designed_chain_convention": (metrics or {}).get("designed_chain_convention"),
        "n_target_numbering_gaps": inventory.get("n_target_numbering_gaps"),
        "predictor": {
            "model_type": manifest.get("model_type"),
            "msa_mode": manifest.get("msa_mode"),
            "num_models": manifest.get("num_models"),
            "num_recycle": manifest.get("num_recycle"),
            "host_url_passed": manifest.get("host_url_passed"),
            "network_expected": manifest.get("network_expected"),
            "synthetic": manifest.get("synthetic", False),
        },
        "seed": inventory.get("seed"),
        "hf_ckpt_revision": inventory.get("hf_ckpt_revision"),
        "rfd4_proteina_version": inventory.get("rfd4_proteina_version"),
        "securedna_screened": False,
        "null_reasons": reasons,
        "designs": rows,
    }
    with open("./outputs/design_scores.json", "w") as handle:
        json.dump(payload, handle, indent=2)
    with open("./outputs/design_scores.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            f"  {row['design']}: {row['status']}"
            + (f" sc_rmsd={row['sc_rmsd_A']} A" if row["sc_rmsd_A"] is not None else "")
            + (f" plddt_binder={row['plddt_binder']}"
               if row["plddt_binder"] is not None else "")
            + (f" plddt_target={row['plddt_target']}"
               if row["plddt_target"] is not None else "")
            + (f" iptm={row['iptm']}" if row["iptm"] is not None else "")
            + ("  PASS" if row["passes"] else "")
            + (f"  [{row['reason'][:120]}]" if row.get("reason") else ""),
            flush=True,
        )

    low_target = [r for r in scored
                  if r["plddt_target"] is not None
                  and r["plddt_target"] < LOW_TARGET_PLDDT_WARN]
    if low_target:
        named = ", ".join(
            "{0} ({1})".format(r["design"], r["plddt_target"]) for r in low_target
        )
        print(
            f"WARNING: {len(low_target)} design(s) sit against a target modelled below "
            f"{LOW_TARGET_PLDDT_WARN} pLDDT: {named}. "
            "No interface number computed against a target that was not modelled means "
            "anything, and given §4 this is the most likely cause of bad Phase 2 figures. "
            "This number never gates: 029's wrong diagnosis came from not having it, and "
            "the fix was to report it rather than to filter on it (2359ffa4).",
            flush=True,
        )
    sources = inventory.get("designed_chain_sources") or []
    if "assumed-first-chain" in sources:
        print(
            "WARNING: at least one design's binder chain was ASSUMED, not resolved, so "
            "every number above for it describes whichever chain came first in the file. "
            "The value travels into design_scores.json and the report (§8.2).",
            flush=True,
        )
    if not scored:
        fail(
            f"No design was scored: {payload['n_failed']} failed and "
            f"{payload['n_skipped']} were skipped. Every row and its reason is in "
            "design_scores.json. Reporting success with nothing scored would make an "
            "empty report look like a clean one."
        )

    print(
        f"Wrote design_scores.json and design_scores.csv: {len(passing)} of "
        f"{len(scored)} scored design(s) pass "
        f"(sc_rmsd_A <= {max_rmsd}, plddt_binder >= {min_plddt}). These thresholds are "
        "literature conventions, not calibrated for this workflow.",
        flush=True,
    )


if __name__ == "__main__":
    main()
