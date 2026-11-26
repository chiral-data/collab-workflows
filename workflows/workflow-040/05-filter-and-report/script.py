#!/usr/bin/env python3
"""Rank what passed, export it, and show the numbers the ranking rests on.

Contract
--------
Reads   design_scores.json        node 04, already sorted and already decided
Writes  passing_sequences.fasta   the designs that passed, tagged
        validation_metrics.csv    the full table, failures included
        report.html               the artifact that leaves the run directory

This node renders. It does not decide: node 04 sorted, node 04 applied the thresholds,
and re-applying them here would be two nodes deciding the same thing and eventually
disagreeing about the same design (179003b9). `passes` is read, never recomputed.

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

The fix for cb418341 is not try/catch, it is not needing JavaScript
-------------------------------------------------------------------
037's report called renderTable, renderPlot and renderRunInfo as a bare sequence. Plotly
comes from a CDN, so with no network renderPlot threw, and that exception aborted the rest
of the script -- taking the Run block with it, which is where the seed, the revision and
"SecureDNA screened: no" live. The table still rendered, so the page looked fine while the
disclaimer had silently gone. 037 fixed it by isolating each renderer and drawing
provenance first.

040 goes one step further, because the lesson generalises: the parts of this page that
matter are written as STATIC HTML, server-side, by this node. The provenance block, the
biosafety notice, the thresholds and the whole metrics table are in the file. There is no
arrangement of a missing CDN, a disabled script or a hostile cache that can remove them,
and the offline case is checkable by reading the file rather than by driving a browser.

JavaScript is enhancement only -- sorting the table -- and it is still isolated in 037's
RENDERERS loop with provenance first, so the pattern is carried rather than discarded. The
stylesheet is the only CDN reference left, and losing it costs appearance and nothing else.

DATA is embedded as well, for the suite: a check that greps a report's prose passes on a
page that merely mentions a word. The table and the DATA payload are generated from the
same column list, so they cannot disagree.

Rows beyond top_n
-----------------
top_n marks the shortlist, and the report says so. It does not control an embedded
structure viewer, which is a deliberate narrowing of the spec's revision 1 -- see the note
below. 037 made every metrics row clickable while embedding only top_n structures, so
clicking row 11 of 10 silently displayed design 1 (7fba161c). Nothing here is clickable,
so there is nothing to mis-dispatch.

Why Phase 1's report embeds no structures
-----------------------------------------
Revision 1 of the spec had top_n mean "structures embedded in the report". It does not,
and the reason is worth more than the feature.

The comparison a reader actually wants is the refolded binder laid over the designed one.
That is not a structure, it is a SUPERPOSITION -- and the superposition is node 04's, the
rotation it already computed to produce sc_rmsd_A. A viewer here would either re-derive
it, which is the defect this whole workflow is built to avoid (§8.1), or show two
structures in unrelated frames and invite the reader to eyeball an RMSD, which is worse
than showing nothing: it looks like evidence.

It also costs two passthroughs. The PDBs reach node 04 and stop there, so putting them on
this page means adding them to node 04's declared outputs and this node's inputs, to carry
files neither node reads.

So Phase 1's report carries numbers and provenance, which is what §7 says its job is while
the ranking metric is uncalibrated. Writing superposed coordinates from node 04 and
viewing them here is a well-scoped later change, and it belongs with Phase 3's
calibration, where there is something to look at.

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
    Phase 2 complex inputs -- so the tag is checked for one at runtime rather than relying
    on the next person to remember.

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
  * the biosafety tag contains a colon, or any FASTA line would start with ';'
  * a design that passed carries no sequence
"""
import csv
import html
import json
import math
import os
import sys

# Carried on every exported record's description line. passing_sequences.fasta is the only
# artifact built to leave this workflow, and a sequence that has passed validation is MORE
# likely to be synthesised than one that has not.
#
# A description, not a leading ';' comment banner: Biopython 1.85 accepts comment lines
# with a DeprecationWarning saying it will become a ValueError, 1.87 already raises, and
# workflow-026 calls SeqIO.parse with no try/except. And no colon anywhere, because 026
# reads a colon in a record id or description as a chain-break marker and flips the run to
# multimer pairing (7709809b).
UNSCREENED_TAG = (
    "UNSCREENED - passed in-silico self-consistency only, no binding evidence, "
    "requires biosecurity review before synthesis"
)

# Only a stylesheet, and only for appearance. Nothing on this page needs it to be
# readable, which is the whole point of rendering the content server-side (cb418341).
BOOTSTRAP_CSS = "https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css"

# (key in design_scores.json, column heading). The table and the DATA payload are both
# generated from this list, so a check that reads either one is reading the same thing.
# Every metric that gates is here, as a column, because the report is the artifact that
# leaves the run directory and 037's did not carry the evidence it was built on (179003b9).
COLUMNS = [
    ("design", "Design"),
    ("status", "Status"),
    ("passes", "Passes"),
    ("sc_rmsd_A", "sc_rmsd (A)"),
    ("plddt_binder", "pLDDT binder"),
    ("plddt_target", "pLDDT target"),
    ("plddt_binder_in_complex", "pLDDT binder in cplx"),
    ("iptm", "ipTM"),
    ("pae_interaction", "PAE interface"),
    ("ptm", "pTM"),
    ("binder_chain", "Binder chain"),
    ("designed_chain_source", "Chain resolved by"),
    ("binder_length", "Length"),
    ("sc_rmsd_n_residues", "RMSD over"),
    ("n_unknown", "X residues"),
    ("gen_compactness_ratio", "037 compactness"),
    ("gen_overextended", "037 overextended"),
    ("gen_chain_breaks", "037 chain breaks"),
    ("gen_steric_clash", "037 clash"),
    ("failed_thresholds", "Failed on"),
    ("reason", "Reason"),
]

# What each number does and does not establish, rendered beside it. 037 presented its
# compactness ratio as a folded-or-not verdict until it was measured against real
# structures; the threshold survived and the claim was narrowed (7968abeb).
CLAIMS = [
    ("sc_rmsd_A", "max_sc_rmsd_a",
     "literature convention, uncalibrated for this workflow",
     "Gates. The refolded binder's Ca RMSD against the backbone it was designed as, "
     "after superposition. Under the cutoff this says the sequence encodes its own "
     "structure. It says NOTHING about binding."),
    ("plddt_binder", "min_plddt_binder",
     "literature convention, uncalibrated for this workflow",
     "Gates. The predictor's own confidence in the binder's fold. Confidence in a fold, "
     "not evidence of one, and not evidence of an interface."),
    ("plddt_target", None,
     "no threshold, by decision",
     "Diagnostic only, and deliberately so. No threshold and no contribution to pass or "
     "fail. It is here because 029's absence of this number produced a wrong diagnosis -- "
     "interface filters read as unreachable when the real cause was a target that had "
     "been modelled blind (2359ffa4)."),
    ("plddt_binder_in_complex", None,
     "Phase 2, diagnostic",
     "The binder's confidence predicted IN CONTEXT, which is not the same number as "
     "plddt_binder and does not replace it. plddt_binder is the binder predicted ALONE, "
     "which is the self-consistency question. 029 had only one of the two and so could "
     "not tell a binder that folds alone but not against its target from one that folds "
     "in neither case."),
    ("iptm", None,
     "Phase 2, reported, does not gate by default",
     "Absent from a monomer prediction entirely -- not zero, absent -- so a blank here on "
     "a Phase 1 run is the normal case and not a failure. On a Phase 2 run it is measured "
     "and REPORTED; it decides nothing unless gate_on_interface was set, and the Run "
     "table below says which it was. No control has yet been run against this metric."),
    ("pae_interaction", None,
     "provisional band, not a gate",
     "Phase 2. Mean predicted aligned error over BOTH off-diagonal blocks. A value with a "
     "provisional band rather than a gate: its published cutoff comes from protocols "
     "using AF2 initial guess, which this workflow does not use, and it is floored by how "
     "well the TARGET modelled -- a target at ~42 pLDDT moved the floor from 7.6 to 15.4 A "
     "and put 10.0 out of reach for a whole batch (8e4772a7). Read the pLDDT target column "
     "beside it, always."),
]


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


def show(value):
    if value is None or value == "":
        return ""
    if isinstance(value, bool):
        return "yes" if value else "no"
    return str(value)


def write_fasta(path, rows):
    """The designs that passed, tagged, one record each.

    Checked rather than assumed: the tag is searched for a colon and every written line is
    searched for a leading ';'. Both constraints are about what the NEXT reader of this
    file does with it, and neither is visible in the file itself, so a comment would not
    have been enough (7709809b).
    """
    if ":" in UNSCREENED_TAG:
        fail(
            "the biosafety tag contains a colon. 026 reads a colon in a record id or "
            "description as a chain-break marker, and 040 uses colons deliberately in its "
            "Phase 2 complex inputs, so this is not a hypothetical collision (7709809b)."
        )

    lines = []
    for row in rows:
        sequence = row.get("binder_sequence")
        if not sequence:
            fail(
                f"'{row['design']}' passed but carries no sequence. Node 04 puts "
                "binder_sequence in every row precisely so this node does not have to "
                "reach past its own dependency for it."
            )
        description = (
            f"sc_rmsd_A={show(row.get('sc_rmsd_A'))} "
            f"plddt_binder={show(row.get('plddt_binder'))} "
            f"max_sc_rmsd_A={show(row.get('max_sc_rmsd_a'))} "
            f"min_plddt_binder={show(row.get('min_plddt_binder'))} | {UNSCREENED_TAG}"
        )
        lines.append(f">{row['design']} {description}")
        for at in range(0, len(sequence), 60):
            lines.append(sequence[at:at + 60])

    for line in lines:
        if line.startswith(";"):
            fail(
                f"a FASTA line would begin with ';': {line[:60]!r}. Biopython 1.85 accepts "
                "comment lines with a DeprecationWarning saying it will become a "
                "ValueError, and 1.87 already raises, so this would break the next "
                "workflow to read the file and be diagnosed as a bug in that workflow."
            )
        if line.startswith(">") and ":" in line:
            fail(
                f"a FASTA header would contain a colon: {line[:80]!r}. A colon in a record "
                "id or description is read as a chain-break marker (7709809b)."
            )

    with open(path, "w") as handle:
        for line in lines:
            handle.write(line + "\n")
    return len(rows)


def table_html(rows, top_n):
    head = "".join(f"<th>{html.escape(label)}</th>" for _, label in COLUMNS)
    body = []
    for index, row in enumerate(rows, start=1):
        classes = []
        if row.get("passes"):
            classes.append("pass")
        if row.get("status") != "scored":
            classes.append("unscored")
        if index > top_n:
            classes.append("below-shortlist")
        cells = "".join(
            f"<td>{html.escape(show(row.get(key)))}</td>" for key, _ in COLUMNS
        )
        body.append(
            "<tr class=\"" + " ".join(classes) + f"\"><td>{index}</td>{cells}</tr>"
        )
    return (
        "<table id=\"metrics\"><thead><tr><th>#</th>" + head + "</tr></thead>"
        "<tbody>" + "".join(body) + "</tbody></table>"
    )


def provenance_html(scores):
    predictor = scores.get("predictor") or {}
    thresholds = scores.get("thresholds") or {}
    pairs = [
        ("Designs", scores.get("n_designs")),
        ("Scored", scores.get("n_scored")),
        ("Passing", scores.get("n_passing")),
        ("Failed prediction", scores.get("n_failed")),
        ("Skipped", scores.get("n_skipped")),
        ("Phase", scores.get("phase")),
        ("max_sc_rmsd_A", thresholds.get("max_sc_rmsd_a")),
        ("min_plddt_binder", thresholds.get("min_plddt_binder")),
        ("min_iptm", thresholds.get("min_iptm")),
        ("max_pae_interaction", thresholds.get("max_pae_interaction")),
        ("Thresholds calibrated", scores.get("thresholds_are_calibrated")),
        ("Gating metrics", ", ".join(scores.get("gating_metrics") or [])),
        ("Diagnostic only", ", ".join(scores.get("diagnostic_metrics") or [])),
        ("Interface scored", scores.get("score_interface")),
        ("Interface gated", scores.get("interface_gated")),
        ("Designs with an interface", scores.get("n_interface_scored")),
        ("Chain resolved by", ", ".join(scores.get("designed_chain_sources") or [])),
        ("Target numbering gaps", scores.get("n_target_numbering_gaps")),
        ("Predictor model", predictor.get("model_type")),
        ("MSA mode", predictor.get("msa_mode")),
        ("Models / recycles",
         f"{predictor.get('num_models')} / {predictor.get('num_recycle')}"),
        ("--host-url passed", predictor.get("host_url_passed")),
        ("Network expected", predictor.get("network_expected")),
        ("Synthetic predictions", predictor.get("synthetic")),
        ("RFD4 seed", scores.get("seed")),
        ("RFD4 checkpoint", scores.get("hf_ckpt_revision")),
        ("RFD4 version", scores.get("rfd4_proteina_version")),
        ("SecureDNA screened", scores.get("securedna_screened")),
        ("Ordering", scores.get("order")),
    ]
    return "".join(
        f"<tr><th>{html.escape(label)}</th><td>{html.escape(show(value))}</td></tr>"
        for label, value in pairs
    )


def claims_html(scores):
    thresholds = scores.get("thresholds") or {}
    rows = []
    for metric, threshold_key, provenance, claim in CLAIMS:
        value = thresholds.get(threshold_key) if threshold_key else None
        rows.append(
            f"<tr><th>{html.escape(metric)}</th>"
            f"<td>{html.escape(show(value))}</td>"
            f"<td>{html.escape(provenance)}</td>"
            f"<td>{html.escape(claim)}</td></tr>"
        )
    return "".join(rows)


def nulls_html(scores):
    rows = []
    for design, reasons in sorted((scores.get("null_reasons") or {}).items()):
        for metric, why in sorted((reasons or {}).items()):
            if not why:
                continue
            rows.append(
                f"<tr><td>{html.escape(design)}</td>"
                f"<td>{html.escape(metric)}</td>"
                f"<td>{html.escape(str(why))}</td></tr>"
            )
    if not rows:
        return "<tr><td colspan=\"3\">Every metric was measured for every design.</td></tr>"
    return "".join(rows)


PAGE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>workflow-040 validation report</title>
<link rel="stylesheet" href="{css}">
<style>
  body {{ margin: 0; padding: 1.5rem; font-family: system-ui, sans-serif; }}
  h1, h2 {{ margin: 1.4rem 0 0.6rem; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 0.86rem; }}
  th, td {{ border: 1px solid #c9ced4; padding: 0.3rem 0.45rem; text-align: left;
            vertical-align: top; }}
  thead th {{ background: #eef1f4; position: sticky; top: 0; }}
  tr.pass td {{ background: #eaf7ec; }}
  tr.unscored td {{ background: #fdf3e7; }}
  tr.below-shortlist td {{ color: #666; }}
  .notice {{ border: 2px solid #a4262c; background: #fdf0f0; padding: 0.8rem 1rem;
             margin: 1rem 0; }}
  .scroll {{ max-height: 70vh; overflow: auto; }}
  caption {{ caption-side: bottom; text-align: left; padding-top: 0.4rem;
             color: #444; font-size: 0.82rem; }}
  @media (prefers-color-scheme: dark) {{
    body {{ background: #15181c; color: #e6e9ec; }}
    th, td {{ border-color: #3a4048; }}
    thead th {{ background: #23272d; }}
    tr.pass td {{ background: #1b2a1e; }}
    tr.unscored td {{ background: #2b2317; }}
    tr.below-shortlist td {{ color: #9aa2ab; }}
    .notice {{ background: #2a1718; border-color: #d1686d; }}
    caption {{ color: #b6bdc5; }}
  }}
</style>
</head>
<body>

<h1>workflow-040 &mdash; binder design validation</h1>

<div class="notice">
<strong>These sequences are unscreened.</strong> Nothing in this workflow performs
biosecurity screening, and nothing in it is experimental evidence. A design that passes
every threshold here is a better hypothesis, not a result &mdash; and a sequence that has
passed validation is <em>more</em> likely to be synthesised than one that has not, so this
notice matters more at this stage than at generation.
<br><br>
<strong>SecureDNA screened: {screened}.</strong>
<br><br>
Every record in <code>passing_sequences.fasta</code> carries this tag, verbatim, on its
description line:
<br>
<code>{tag}</code>
</div>

<h2>What this report establishes, and what it does not</h2>
<p>
Self-consistency: each designed sequence was refolded by a predictor that never saw the
design, and the refold compared against the backbone the sequence was designed as. Passing
means <em>the sequence encodes its own structure</em>. It does not mean the design binds
anything.
</p>
<p>
Interface scoring, when it ran, is an ADDITION to that and not a replacement for it. The
Run table says whether it ran and whether it gated. By default it does not gate: ipTM and
PAE interface are measured and reported, and <em>pass</em> still means self-consistency
alone. Every threshold here is a literature convention rather than a figure calibrated
against this workflow own controls, and the one control set that has run measured what
that is worth -- a 115-residue poly-alanine, a sequence encoding nothing, scored 95.8
pLDDT and cleared the binder-confidence cutoff. Self-consistency RMSD was the only metric
that rejected it.
</p>

<h2>Run</h2>
<table><tbody>
{provenance}
</tbody>
<caption>
Rendered server-side, in the file. Losing a stylesheet changes how this page looks and
nothing about what it says &mdash; 037's equivalent block was drawn by JavaScript and
vanished with it when a CDN was unreachable (cb418341).
</caption>
</table>

<h2>Metrics, their thresholds, and the claim each one supports</h2>
<table><thead><tr>
<th>Metric</th><th>Threshold</th><th>Where the threshold comes from</th>
<th>What the number establishes</th>
</tr></thead><tbody>
{claims}
</tbody></table>

<h2>Designs</h2>
<p>
{n_designs} design(s), in the numeric order node 04 decided. The first {top_n} are the
shortlist (<code>top_n</code>); rows below it are greyed and are otherwise identical
&mdash; nothing on this page is clickable, so there is nothing that can dispatch to the
wrong design (7fba161c). Designs that could not be scored are kept as rows with their
reason, because a report that omits them reads as a batch where everything worked.
</p>
<div class="scroll">
{table}
</div>

<h2>Metrics that could not be measured, and why</h2>
<table><thead><tr><th>Design</th><th>Metric</th><th>Why it is blank</th></tr></thead>
<tbody>
{nulls}
</tbody>
<caption>
Blank means not measured. It never means zero: a missing ipTM reported as 0.0 reads as a
measured, terrible interface, and a missing PAE reported as infinity fails its threshold
and so renders an absent measurement as a verdict.
</caption>
</table>

<script>
const DATA = {payload};

function enhanceTable() {{
  // Enhancement only. Everything this page says is already in the HTML above; if this
  // throws, or never runs, the report is still complete. That is the difference between
  // isolating a renderer and not depending on one (cb418341).
  const table = document.getElementById("metrics");
  if (!table) throw new Error("no metrics table");
  const headers = table.querySelectorAll("thead th");
  headers.forEach((header, index) => {{
    header.style.cursor = "pointer";
    header.title = "Click to sort. Node 04 decided the order; this only changes the view.";
    header.addEventListener("click", () => {{
      const body = table.querySelector("tbody");
      const rows = Array.from(body.querySelectorAll("tr"));
      const descending = header.dataset.sorted === "asc";
      rows.sort((a, b) => {{
        const x = a.children[index].textContent.trim();
        const y = b.children[index].textContent.trim();
        const nx = parseFloat(x), ny = parseFloat(y);
        const both = !isNaN(nx) && !isNaN(ny);
        const cmp = both ? nx - ny : x.localeCompare(y);
        return descending ? -cmp : cmp;
      }});
      headers.forEach((h) => {{ delete h.dataset.sorted; }});
      header.dataset.sorted = descending ? "desc" : "asc";
      rows.forEach((row) => body.appendChild(row));
    }});
  }});
}}

var RENDERERS = [
  ['sortable table', enhanceTable]
];
// 037's isolation pattern, kept rather than discarded: each renderer runs on its own, so
// one failing cannot abort the rest of the script. In 040 there is nothing essential left
// for it to take with it, which is the stronger version of the same fix.
for (var i = 0; i < RENDERERS.length; i++) {{
  try {{
    RENDERERS[i][1]();
  }} catch (err) {{
    console.error('workflow-040 report: ' + RENDERERS[i][0] + ' failed', err);
  }}
}}
</script>
</body>
</html>
"""


def main():
    top_n = positive_param("PARAM_TOP_N", "top_n", 10, int)

    os.makedirs("./outputs", exist_ok=True)
    path = "./inputs/design_scores.json"
    if not os.path.exists(path):
        fail(
            "No 'design_scores.json' in ./inputs/. This node renders node 04's decisions; "
            "it does not compute them, and there is nothing here it could fall back to."
        )
    try:
        with open(path) as handle:
            scores = json.load(handle)
    except Exception as exc:
        fail(f"Could not read '{path}': {exc}")

    rows = scores.get("designs") or []
    if not rows:
        fail(
            "design_scores.json carries no designs. An empty report is indistinguishable "
            "from a clean one."
        )

    # Read, never recomputed. Two nodes applying the same threshold eventually disagree
    # about the same design (179003b9).
    passing = [row for row in rows if row.get("passes")]

    n_exported = write_fasta("./outputs/passing_sequences.fasta", passing)

    with open("./outputs/validation_metrics.csv", "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["rank"] + [label for _, label in COLUMNS])
        for index, row in enumerate(rows, start=1):
            writer.writerow([index] + [show(row.get(key)) for key, _ in COLUMNS])

    payload = {
        # The suite reads this rather than the prose: a check that greps a page for a word
        # passes on a page that merely mentions it (037's report_payload).
        "columns": [list(pair) for pair in COLUMNS],
        "top_n": top_n,
        "n_designs": len(rows),
        "n_passing": len(passing),
        "n_exported": n_exported,
        "thresholds": scores.get("thresholds"),
        "thresholds_are_calibrated": scores.get("thresholds_are_calibrated"),
        "gating_metrics": scores.get("gating_metrics"),
        "diagnostic_metrics": scores.get("diagnostic_metrics"),
        "designed_chain_sources": scores.get("designed_chain_sources"),
        "securedna_screened": scores.get("securedna_screened"),
        "biosafety_tag": UNSCREENED_TAG,
        "predictor": scores.get("predictor"),
        "null_reasons": scores.get("null_reasons"),
        "designs": rows,
    }

    page = PAGE.format(
        css=BOOTSTRAP_CSS,
        screened="no" if not scores.get("securedna_screened") else "yes",
        # The tag itself, in the page, not only in the embedded payload. A reader looking
        # at the report should be able to see exactly what left the workflow attached to
        # those sequences without opening the FASTA.
        tag=html.escape(UNSCREENED_TAG),
        provenance=provenance_html(scores),
        claims=claims_html(scores),
        table=table_html(rows, top_n),
        nulls=nulls_html(scores),
        n_designs=len(rows),
        top_n=top_n,
        payload=json.dumps(payload),
    )
    with open("./outputs/report.html", "w") as handle:
        handle.write(page)

    for row in rows[:top_n]:
        print(
            f"  {row['design']}: {row.get('status')}"
            + (f" sc_rmsd={show(row.get('sc_rmsd_A'))} A"
               if row.get("sc_rmsd_A") is not None else "")
            + (f" plddt_binder={show(row.get('plddt_binder'))}"
               if row.get("plddt_binder") is not None else "")
            + ("  PASS" if row.get("passes") else "")
            + (f"  [{str(row.get('reason'))[:100]}]" if row.get("reason") else ""),
            flush=True,
        )
    if len(rows) > top_n:
        print(
            f"  ... and {len(rows) - top_n} more row(s) below the shortlist, all of them "
            "in validation_metrics.csv and report.html",
            flush=True,
        )

    unscored = [r for r in rows if r.get("status") != "scored"]
    if unscored:
        print(
            f"WARNING: {len(unscored)} of {len(rows)} design(s) were not scored and are "
            "rendered as rows with their reason rather than omitted. A report that leaves "
            "them out reads as a batch where everything worked (§8.5).",
            flush=True,
        )
    if not passing:
        print(
            "WARNING: no design passed. That is a result, not an error: the thresholds "
            "here are literature conventions and have not been calibrated against this "
            "workflow's own controls (§7), so read the numbers in the table before "
            "concluding anything about the designs. passing_sequences.fasta is empty.",
            flush=True,
        )

    print(
        f"Wrote passing_sequences.fasta ({n_exported} record(s)), "
        f"validation_metrics.csv ({len(rows)} row(s)) and report.html. "
        "Passing means the sequence recovers its own backbone; it does not mean the "
        "design binds.",
        flush=True,
    )


if __name__ == "__main__":
    main()
