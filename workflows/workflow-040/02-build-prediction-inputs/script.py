#!/usr/bin/env python3
"""Write the prediction inputs: one design per file, always.

STUB -- not implemented. Exits 1.

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

prediction_inputs.json, per design
----------------------------------
    design              keyed on, everywhere
    path                the file written
    mode                "mono" in Phase 1
    chain_order         ["binder"] in Phase 1; ["target", "binder"] in Phase 2
    chain_lengths       parallel to chain_order
    expected_chains     how many chains the prediction must come back with

chain_order and chain_lengths exist for one reason: the scores file carries NO chain
boundaries. Measured, from a real ColabFold 1.6.1 run in
workflow-026/sample_outputs/confidence_summary.json -- plddt_per_residue has length 164,
exactly chain_lengths[0], and pae is 164x164, with nothing anywhere naming a chain. So
node 04 must slice the per-residue array by position, and the position depends on the
order chosen HERE. This node decides it, so this node records it.

In Phase 1 the binder is the whole array and the slice is trivial. Record it anyway. When
Phase 2 adds TARGET:BINDER the binder becomes the TAIL, and workflow-029 slices
plddt[:binder_length] because its binder is chain A (04_fold_and_filter/script.py:207).
A node 04 that inherits that assumption would report the target's confidence as the
binder's, which on a well-modelled target makes every design look good.

Colon rules, both directions
----------------------------
A mono input must contain NO colon: a colon makes the predictor treat one design as two
chains. A Phase 2 complex input must contain EXACTLY ONE, with both halves non-empty --
zero silently predicts a monomer and reports its confidence as an interface, which is
026's bug, and two predicts three chains. 026 validates sequences against
ACDEFGHIKLMNPQRSTVWYX and so rejects the only syntax that encodes a complex (d3d707ed);
040 must accept it in the right file and refuse it in the wrong one.

Skips, with a warning rather than silently
------------------------------------------
A design whose sequence is empty or all-unknown is skipped rather than written as junk,
and unknown residues are counted: ColabFold folds X as a sequence gap, which depresses
confidence for reasons that have nothing to do with the design (54bf21b6).

Re-running
----------
Must be safe or refused. Workflow-029's chain swap was not idempotent, so re-running over
a non-empty output directory corrupted backbones it had already fixed -- and it fired for
real, during a GPU smoke test, on a directory left by an earlier run (4d1a523e). A
refusal must leave what is already there untouched.

Fails when
----------
  * a binder sequence carries a residue outside the accepted alphabet
  * a binder length falls outside the configured bounds
  * score_interface is true (Phase 2 is not built -- say so, do not ignore it)
  * the output directory holds a previous run and the node cannot proceed safely
"""
import sys

print(
    "ERROR: 02-build-prediction-inputs is not implemented yet.\n"
    "       Scaffold only -- see workflows/workflow-040/validation_pipeline.md §6\n"
    "       and workflows/workflow-040/tests/fixtures.json (section 'ingest').",
    file=sys.stderr,
)
sys.exit(1)
