#!/usr/bin/env python3
"""Refold each designed sequence with a predictor that never saw the design.

STUB -- not implemented. Exits 1. GPU.

Contract
--------
Reads   mono_*.fasta, prediction_inputs.json
Writes  <design>_scores_rank_*.json   flat, prefixed per design -- see below
        <design>_*rank_*.pdb          all five ranks; node 04 selects rank_001
        prediction_manifest.json      per design: status, output file names, chain_order,
                                      chain_lengths, argv, elapsed, failure reason

Flat output, prefixed per design
--------------------------------
Silva's ports match file NAMES, not paths: `silva validate` refuses an input like
`predictions/*` because it "can never match". So a per-design subdirectory would be
flattened on the way into node 04's inputs anyway, and the separation between designs has
to live in the filename instead.

That is not a workaround, it is the lesson. Flattening per-design output into one
directory is exactly how workflow-026 conflates designs -- its 02_predict "copies all
colabfold output flat into ./outputs/, so N designs' files land in one directory
undistinguished", and its 03_analyze_confidence then globs across all of them and takes
job_id from whichever sorts first (818954e9). The defence is that ColabFold names its
outputs after the record id and node 02 writes exactly one record per design, so every
file already carries the design it belongs to.

Run colabfold_batch into a scratch directory per design, then move its files out under
that prefix. On a name collision, keep a disambiguating prefix rather than overwriting:
workflow-037 node 03 flattens the same way and had to learn this twice -- first that the
real layout nests one level deeper than it assumed, then that flattening a tree can
collide and the old code would have overwritten one file with another silently
(cb3f6d23).

The manifest carries chain_order and chain_lengths forward from node 02 so node 04 reads
them from its own direct dependency rather than taking an extra edge back to node 02.
Passing the fact forward is the rule; re-deriving it is the defect (CLAUDE.md).

argv
----
    colabfold_batch <input> <outdir>
      --model-type    alphafold2_ptm          (mono; alphafold2_multimer_v3 for Phase 2)
      --msa-mode      single_sequence         (PARAM_MSA_MODE)
      --num-models    5                       (PARAM_NUM_MODELS)
      --num-recycle   3                       (PARAM_NUM_RECYCLE)

--model-type is not optional. The multimer model is what produces iptm at all, and
running a single chain through it is not the monomer prediction anyone means. The real
sample confirms alphafold2_ptm for the monomer path, in its own filename:
sequences_unrelaxed_rank_001_alphafold2_ptm_model_4_seed_000.pdb.

--host-url is NEVER passed. Workflow-026 passes it unconditionally, defaulting to
https://api.colabfold.com (02_predict/script.py:58), which is how a workflow reaches a
third-party service without anyone deciding that it should. The weights are baked into
the image -- workflow-026's Dockerfile runs `python -m colabfold.download` at build time
-- so single_sequence mode has no reason to touch the network. Assert that it does not.

Phase 1 is monomer-only and single_sequence is CORRECT here rather than a compromise: a
de novo binder has no homologs, and designed proteins are idealised, which is why the
self-consistency convention works without an alignment. It is wrong for a complex, for a
reason that is about the target rather than the binder -- see validation_pipeline.md §4.

Timeout: kill the process group, not the child
----------------------------------------------
subprocess.run's own timeout is not enough. Workflow-037 found that generation spawns
workers which go on holding the GPU after the direct child is signalled: a run outlived
the harness that started it, the client was killed, the container was not, and the device
had to be freed by hand (755d4af2). So: Popen in its own process group, under
PARAM_PREDICT_TIMEOUT_S, with kill_tree() escalating SIGTERM then SIGKILL. Verified
against a child that installs SIG_IGN for SIGTERM, which waits out the grace period and
dies of SIGKILL at rc=-9 (55edd8fc).

A stalled prediction holds the device exactly the same way a stalled generation does. The
budget here is unmeasured for this workflow (§11.5): recycles multiply it, and
workflow-026 records 6-12 recycles as the recommendation for multimers against the 3 that
suits a monomer.

Retain the output
-----------------
Stream it and keep it, so a failure is diagnosed rather than inferred. Both halves of
that have cost someone a day already: 037 could only infer an allocator OOM from
warnings the CUDA caching allocator prints while carrying on (755d4af2), and 029
discarded a prediction's ColabFold output unless the design passed, so a failing funnel
could not afterwards distinguish a bad design from a target that was never modelled --
an ambiguity that produced one wrong diagnosis and had to be resolved by computing the
missing number by hand, outside the pipeline (2359ffa4).

Assertions
----------
  * the prediction came back with prediction_inputs.json's expected_chains. A Phase 2
    complex that returns one chain means the ':' was dropped somewhere and the result is
    a monomer masquerading as a complex -- which is 026's exact failure, where
    independent monomers are labelled a multimer and --pair-mode unpaired_paired is set,
    and the flag changes nothing because there is no complex. Fail loudly.
  * a rank_001 scores file and a rank_001 PDB exist per design, found by the design's own
    filename prefix. Never glob across designs.
  * the network was not reached.

Partial failure
---------------
Recorded per design; the batch continues. A design with no prediction is named in the
manifest and must reach the report as a failed row with a stated reason, not vanish from
it (§8.5).

Fails when
----------
  * no GPU (fail at scheduling via [resources], and again at check_gpu() to be legible)
  * every design failed to predict
  * the chain-count assertion fails for any design
"""
import sys

print(
    "ERROR: 03-predict is not implemented yet.\n"
    "       Scaffold only -- see workflows/workflow-040/validation_pipeline.md §6\n"
    "       and workflows/workflow-040/tests/fixtures.json (section 'predict').",
    file=sys.stderr,
)
sys.exit(1)
