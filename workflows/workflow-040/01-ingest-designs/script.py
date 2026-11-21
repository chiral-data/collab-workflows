#!/usr/bin/env python3
"""Read the designed structures and say, per design, which chain is the binder.

STUB -- not implemented. Exits 1. See AI_WORKFLOW_GUIDE.md §4: nodes land one commit at
a time.

This node exists so that nothing downstream ever has to decide which chain is the
design. That decision was already made upstream, twice: workflow-037's node 03 recorded
the target's chain ids and its node 04 resolved the designed chain from them. Both
answers are in design_metrics.json. Read them.

Contract
--------
Reads   *_generated.cif          workflow-037 node 03, one structure per design
        design_metrics.json      workflow-037 node 04, optional
Writes  design_inventory.json    one entry per design
        design_metrics.json      passed through unchanged for nodes 04 and 05 to merge

design_inventory.json, per design:

    design                  the design's name, which is what everything downstream keys on
    binder_chain            chain id of the designed chain
    binder_sequence         one-letter, in chain order
    binder_residue_ids      label_seq_id values, IN CHAIN ORDER
    binder_ca               Ca coordinates, in that same order
    target_chains           [{chain, sequence, residue_ids, gaps}, ...]
    designed_chain_source   carried through VERBATIM from design_metrics.json
    n_unknown               residues with no one-letter code, scored X

The Ca coordinates live here on purpose. Node 04 needs the designed backbone to compute
an RMSD against, and the alternative -- passing the CIFs through and letting node 04
parse them -- would make node 04 resolve the binder chain a second time. That is the
defect class this repo keeps fixing. One node reads the structures; nobody else parses
them.

designed_chain_source
---------------------
Three values, and they are workflow-037's, not ours (04-analyze/script.py:248). Do not
invent a fourth and do not collapse them:

    manifest-target-chains   resolved by set difference against the recorded target
                             chains. The good path.
    single-chain             one chain in the file, so it is the design by construction:
                             a monomer or motif_scaffold run. Valid for self-consistency.
                             Node 02 must NOT write a complex input for it -- there is no
                             target to score an interface against.
    assumed-first-chain      037 declined to answer and fell back to chain order.
                             Proceed, say so loudly, and carry the value forward.

Absent design_metrics.json is the stated-fallback path, not an error: record
assumed-first-chain, warn, and continue.

Why the fallback matters more than it looks
-------------------------------------------
037's first-chain fallback was in force on EVERY real run for weeks. Its node 03 recorded
target_chains by listing every label_asym_id in target.cif, which on 4ZQK is the target
chain, the chain the contig dropped, a sodium and two solvent groups -- so the designed
chain was in that list, the set difference found no candidate, and the fallback fired
every time (ef74646e). And it was right every time, because the design does come out
first. Nothing failed. Nothing looked wrong.

So this node's own test is not that the fallback works. It is that the good path was
taken: tests/fixtures.json, ingest/resolved-chain-good-path.

Fails when
----------
  * no *_generated.cif inputs
  * design_metrics.json names a chain the structure does not contain
  * a design's binder chain resolves to zero residues

Does NOT fail when design_metrics.json is absent, or when a target chain carries internal
gaps -- gaps are measured and reported, because node 03 has to disclose that it folds the
chain whole and fuses them. Trimming around a gap was measured to be much worse
(8e4772a7, §8.12).
"""
import sys

print(
    "ERROR: 01-ingest-designs is not implemented yet.\n"
    "       Scaffold only -- see workflows/workflow-040/validation_pipeline.md §6\n"
    "       and workflows/workflow-040/tests/fixtures.json (section 'ingest').",
    file=sys.stderr,
)
sys.exit(1)
