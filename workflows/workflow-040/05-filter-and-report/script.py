#!/usr/bin/env python3
"""Rank what passed, export it, and show the numbers the ranking rests on.

STUB -- not implemented. Exits 1.

Contract
--------
Reads   design_scores.json        node 04, already sorted and already decided
Writes  passing_sequences.fasta   the designs that passed, tagged
        validation_metrics.csv    the full table, failures included
        report.html               the artifact that leaves the run directory

This node renders. It does not decide: node 04 sorted, node 04 applied the thresholds,
and re-applying them here would be two nodes deciding the same thing and eventually
disagreeing about the same design (179003b9).

040 is allowed to rank, which is the point of it -- workflow-037 is not, and says so
(e8326b55). But say what the order means, and until the Phase 3 controls have run, present
any uncalibrated metric as a value with a provisional band rather than as a gate (§7).

Every gating metric is a column
-------------------------------
Not a log line, not an embedded JSON field. Workflow-037 measured CA-CA spacing and
compactness in node 04 and warned about an unfolded design in the node log, while the
report rendered it as an ordinary row -- so the one artifact built to be emailed away from
the run directory was the one that did not carry the evidence (179003b9).

Failed designs appear as rows with blank metrics and a stated reason. A report that omits
what it could not score reads as a batch where everything worked (§8.5).

designed_chain_source is rendered, not just recorded. 037's report prints which path was
taken rather than a fixed phrase, precisely because its fallback had been silently in
force on every real run (ef74646e).

Isolate the renderers, draw provenance first
--------------------------------------------
037's report called renderTable, renderPlot and renderRunInfo as a bare sequence. Plotly
comes from a CDN, so with no network renderPlot threw, and that exception aborted the rest
of the script -- taking the Run block with it, which is where the seed, the revision and
"SecureDNA screened: no" live. The table still rendered, so the page looked fine while the
disclaimer had silently gone (cb418341).

Rows beyond top_n are inert and say why: 037 made every metrics row clickable while
embedding only top_n structures, so clicking row 11 of 10 silently displayed design 1
(7fba161c).

The screening notice matters MORE here than upstream
----------------------------------------------------
A sequence that has passed validation is more likely to be synthesised, not less. 037's
unscreened warning reached the README, the report, the node output and the manifest --
everywhere except the FASTA, the one artifact built to leave the workflow (7709809b).

Two mechanical constraints on the tag, both load-bearing:

  * No leading ';' comment banner. Biopython 1.85 accepts comment lines with a
    DeprecationWarning saying it will become a ValueError; 1.87 already raises. Shipping
    one would break the next workflow to read the file, and it would be diagnosed as a bug
    in that workflow.
  * No colon in the tag. A colon in a record id or description is read as a chain-break
    marker. Sharper for 040 than for 037, because 040 uses colons deliberately in its
    Phase 2 complex inputs -- so check the tag for one at runtime rather than relying on
    the next person to remember.

Narrow the claim to what was measured
-------------------------------------
sc_rmsd_A under its cutoff means the sequence recovers its own backbone. It does not mean
the design binds. 037 presented its compactness ratio as a folded-or-not verdict until
somebody measured it against real structures and found it flags a textbook alpha-helix at
every length in the default binder range; the threshold survived and the claim was
narrowed to "not compact for its length" (7968abeb). Same discipline: the report says what
each number establishes and what it does not.

Fails when
----------
  * design_scores.json is absent or empty
  * top_n is not a positive whole number
"""
import sys

print(
    "ERROR: 05-filter-and-report is not implemented yet.\n"
    "       Scaffold only -- see workflows/workflow-040/validation_pipeline.md §6\n"
    "       and workflows/workflow-040/tests/fixtures.json (section 'scoring').",
    file=sys.stderr,
)
sys.exit(1)
