#!/usr/bin/env python3
"""Measure each refolded design against the backbone it was designed as.

STUB -- not implemented. Exits 1.

This is the node the whole workflow exists for, and the node where a wrong answer looks
most like a right one. Everything here is a number, and every number is plausible.

Contract
--------
Reads   design_inventory.json      node 01: binder residue ids and Ca coords, IN ORDER
        design_metrics.json        workflow-037 node 04, passed through by node 01
        *_scores_rank_*.json       node 03, flat and prefixed per design
        *_rank_*.pdb               node 03, all ranks; select rank_001
        prediction_manifest.json   node 03: chain_order, chain_lengths, per-design status

Resolve each design's files by ITS OWN prefix, named in prediction_manifest.json. Never
glob the directory and never pair files by sort order: everything arrives flat, because
Silva's ports match file names rather than paths, and globbing across designs is the
mechanism of 026's conflation bug -- it takes job_id from whichever scores file sorts
first and reports a confident result about a design that does not exist (818954e9).
Writes  design_scores.json        one row per design, failures included
        design_scores.csv         the same table

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
superposition code being right.

Second: correspondence is by POSITION in the chain, established by node 01 -- never by
residue number. Three renumberings stack: the deposited entry's, RFD4's (it renumbers the
preserved target chain from 1, so barnase in 1BRS, which starts at label 3, arrives
already shifted from what 037 node 02 conditioned on -- ef74646e), and ColabFold's own
renumber-from-1-per-chain. Assert the mapping -- equal lengths, 1:1 -- and refuse a
mismatch. Truncating to the shorter and scoring produces a mis-registered comparison
whose number reads as a mediocre design rather than as a bug.

Slicing pLDDT: read the order, do not assume it
-----------------------------------------------
The scores file carries NO chain boundaries. Measured, from the real ColabFold 1.6.1 run
in workflow-026/sample_outputs/confidence_summary.json: plddt_per_residue is length 164,
exactly chain_lengths[0], pae is 164x164, and nothing names a chain.

So read chain_order and chain_lengths out of prediction_manifest.json, and ASSERT
len(plddt) == sum(chain_lengths) before slicing. That assertion is the only thing standing
between a position assumption and a plausible wrong number. Workflow-029 slices
plddt[:binder_length] because its binder is chain A; Phase 2's input is TARGET:BINDER, so
the binder is the tail.

Never default a missing metric to a number
------------------------------------------
MEASURED, not hypothetical: in that same real monomer run, iptm is ABSENT from all five
models, along with interface_pae_mean, while ptm is present at 0.88-0.90. Phase 1 is
monomer-only, so this is the normal case, not an edge case.

Workflow-029 reads scores.get('iptm', 0.0) and sets pae_interaction = float('inf') when
the matrix is absent. The first reports 0.0 -- a measured-looking bad interface -- for
every design in every run; the second fails the threshold, so a missing measurement is
rendered as a verdict. Both are indistinguishable from a real result downstream.

Workflow-026 does it correctly, in one line
(03_analyze_confidence/script.py:86):

    iptm = float(data["iptm"]) if "iptm" in data else None

Copy 026. A metric that could not be computed is null, with a recorded reason.

Merge by name
-------------
Keyed on the design name, never by position. Workflow-037's metrics rows may arrive in a
different order from the predictions, and merging by position puts one design's generation
metrics against another's scores with every row still plausible. 026's conflation bug is
the same error one level up: it takes job_id from whichever scores file sorts first
(818954e9).

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
  * a threshold parameter is nan, inf, non-numeric or non-positive
  * score_interface is true (Phase 2 is not built)
  * designed_chain_source is single-chain and interface scoring was requested
"""
import sys

print(
    "ERROR: 04-score is not implemented yet.\n"
    "       Scaffold only -- see workflows/workflow-040/validation_pipeline.md §6\n"
    "       and workflows/workflow-040/tests/fixtures.json (section 'scoring').",
    file=sys.stderr,
)
sys.exit(1)
