#!/usr/bin/env python3
"""Node 05 — Report.

Build a single report.html with the metrics table and the top_n structures inlined as
JavaScript, so the run's own data travels with the file.

Bootstrap, Plotly and Mol* come from a CDN, so the page is not self-contained in the stricter
sense and needs a network for styling, the scatter and the viewer. The blocks that matter
without one -- the metrics table, the run provenance and the biosafety alert -- are rendered
independently of those libraries and do not depend on them loading.

The report carries the biosafety notice prominently: the RFD4 CLI and Python API do not run
SecureDNA screening (only the HF Spaces apps do), so designs leaving this workflow are
unscreened and require review before synthesis.

Designs are ordered by interface contacts (descending) for binder-style runs, falling back
to design name when nothing has a target to contact. That ordering is a convenience, not a
ranking: with `~metrics` there is no confidence or affinity estimate to rank by.
"""
import glob
import json
import os
import sys

BOOTSTRAP_CSS = "https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css"
PLOTLY_JS = "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"
MOLSTAR_CSS = "https://cdn.jsdelivr.net/npm/molstar@4.9.0/build/viewer/molstar.css"
MOLSTAR_JS = "https://cdn.jsdelivr.net/npm/molstar@4.9.0/build/viewer/molstar.js"

TABLE_COLUMNS = [
    ("design", "Design"),
    ("designed_chain", "Chain"),
    ("designed_length", "Length"),
    ("hydrophobic_fraction", "Hydrophobic"),
    ("net_charge", "Net charge"),
    ("radius_of_gyration_A", "Rg (Å)"),
    ("n_interface_contacts", "Contacts"),
    ("n_target_residues_contacted", "Target res."),
    ("min_interchain_ca_distance_A", "Min CA-CA (Å)"),
    ("steric_clash", "Clash"),
]


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def int_param(env_name, display_name, default):
    """Read a whole-number parameter, failing the way every other bad parameter does.

    Casting inline at the point of use turns a typo into a raw ValueError traceback, which
    reads as a crash in the workflow rather than as the user error it is.
    """
    raw = os.environ.get(env_name)
    if raw is None or raw.strip() == "":
        return default
    try:
        return int(raw)
    except ValueError:
        fail(f"{display_name} must be a whole number, got '{raw}'.")


def build_html(metrics, structures, top_n):
    designs = metrics.get("designs", [])
    payload = json.dumps(
        {
            "metrics": metrics,
            "structures": structures,
            "columns": TABLE_COLUMNS,
        }
    ).replace("</", "<\\/")

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>RFD4-Proteina Designs</title>
<link rel="stylesheet" href="{BOOTSTRAP_CSS}">
<link rel="stylesheet" href="{MOLSTAR_CSS}">
<script src="{PLOTLY_JS}"></script>
<script src="{MOLSTAR_JS}"></script>
<style>
  body {{ padding: 2rem 1rem; }}
  #viewer {{ position: relative; width: 100%; height: 460px; }}
  .table-sm td, .table-sm th {{ font-size: 0.9rem; }}
  .seq {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
          font-size: 0.8rem; word-break: break-all; }}
</style>
</head>
<body>
<div class="container-xl">

  <h1 class="h3">RFD4-Proteina de novo design</h1>
  <p class="text-secondary">
    {len(designs)} design(s) generated. Metrics are geometric only — generation ran with
    <code>~metrics</code>, so there is no confidence, affinity or self-consistency estimate
    anywhere in this report.
  </p>

  <div class="alert alert-danger" role="alert">
    <strong>These designs are not screened.</strong> RFD4's HuggingFace Spaces screen inputs
    and designs with SecureDNA; the CLI and Python API used by this workflow do not. Every
    sequence here requires biosecurity review before it is ordered or synthesised.
  </div>

  <div class="row g-4">
    <div class="col-lg-7">
      <h2 class="h5">Designs</h2>
      <p class="text-muted small mb-2">
        Ordered by interface contacts, descending, with clashing designs moved to the
        bottom &mdash; a convenience, not a ranking. Generation ran with
        <code>~metrics</code>, so there is no confidence, affinity or self-consistency
        estimate to rank by, and contact counts rise with designed length.
      </p>
      <p class="text-muted small mb-2">
        <strong>Clash</strong> marks a design whose closest inter-chain C&alpha;&ndash;C&alpha;
        distance is under 3.8&nbsp;&Aring; &mdash; shorter than the C&alpha;&ndash;C&alpha;
        virtual bond between two covalently joined residues, so the chains are passing
        through each other rather than touching. Such a design has the highest contact
        count of all, which is why it is ranked last here instead of first. Treat a
        flagged design as discarded.
      </p>
      <div class="table-responsive">
        <table class="table table-sm table-hover align-middle" id="metrics-table">
          <thead><tr id="metrics-head"></tr></thead>
          <tbody id="metrics-body"></tbody>
        </table>
      </div>
      <div id="length-plot" style="height:300px"></div>
    </div>
    <div class="col-lg-5">
      <h2 class="h5">Structure</h2>
      <select class="form-select mb-2" id="structure-select"></select>
      <div id="viewer" class="border rounded"></div>
      <p class="seq mt-2" id="sequence"></p>
    </div>
  </div>

  <h2 class="h5 mt-4">Run</h2>
  <ul class="list-unstyled text-secondary small" id="run-info"></ul>

</div>
<script>
const DATA = {payload};
const TOP_N = {top_n};

function fmt(value) {{
  if (value === null || value === undefined || value === "") return "–";
  return value;
}}

function renderTable() {{
  const head = document.getElementById("metrics-head");
  DATA.columns.forEach(([, label]) => {{
    const th = document.createElement("th");
    th.textContent = label;
    head.appendChild(th);
  }});
  const body = document.getElementById("metrics-body");
  const embedded = new Set(DATA.structures.map((s) => s.name));
  DATA.metrics.designs.forEach((row) => {{
    const tr = document.createElement("tr");
    const shown = embedded.has(row.design);
    DATA.columns.forEach(([key]) => {{
      const td = document.createElement("td");
      td.textContent = fmt(row[key]);
      tr.appendChild(td);
    }});
    if (shown) {{
      tr.style.cursor = "pointer";
      tr.addEventListener("click", () => select(row.design));
    }} else {{
      // Only the first top_n structures are embedded. Clicking a row whose structure is
      // absent used to fall back to the first design, silently showing the wrong one.
      tr.classList.add("text-secondary");
      tr.title = "Structure not embedded in this report — raise top_n to view it";
    }}
    body.appendChild(tr);
  }});
}}

function renderPlot() {{
  const rows = DATA.metrics.designs;
  Plotly.newPlot("length-plot", [{{
    type: "scatter",
    mode: "markers",
    x: rows.map((r) => r.designed_length),
    y: rows.map((r) => r.n_interface_contacts),
    text: rows.map((r) => r.design),
    marker: {{ size: 11 }}
  }}], {{
    margin: {{ t: 10, r: 10, b: 45, l: 55 }},
    xaxis: {{ title: "Designed length (residues)" }},
    yaxis: {{ title: "Interface contacts" }}
  }}, {{ displayModeBar: false, responsive: true }});
}}

let viewer = null;

async function select(name) {{
  const entry = DATA.structures.find((s) => s.name === name);
  if (!entry) return;  // not embedded: leave the viewer where it is
  document.getElementById("structure-select").value = entry.name;
  const row = DATA.metrics.designs.find((r) => r.design === entry.name);
  document.getElementById("sequence").textContent = row ? row.designed_sequence : "";
  if (!viewer) {{
    viewer = await molstar.Viewer.create("viewer", {{
      layoutIsExpanded: false,
      layoutShowControls: false,
      layoutShowSequence: false,
      layoutShowLog: false
    }});
  }}
  await viewer.plugin.clear();
  await viewer.loadStructureFromData(entry.cif, "mmcif");
}}

function renderRunInfo() {{
  const m = DATA.metrics;
  const items = [
    ["Designs", m.n_designs],
    ["Seed", m.seed],
    ["Checkpoint revision", m.hf_ckpt_revision],
    ["rfd4-proteina version", m.rfd4_proteina_version],
    ["Contact cutoff", m.contact_cutoff_angstrom + " Å (CA–CA)"],
    ["Designed chain", m.designed_chain_convention],
    ["SecureDNA screened", m.securedna_screened ? "yes" : "no"],
    ["Structures shown", DATA.structures.length + " of " + m.n_designs + " (top_n = " + TOP_N + ")"]
  ];
  const list = document.getElementById("run-info");
  items.forEach(([label, value]) => {{
    const li = document.createElement("li");
    li.textContent = label + ": " + fmt(value);
    list.appendChild(li);
  }});
}}

const select_el = document.getElementById("structure-select");
DATA.structures.forEach((s) => {{
  const option = document.createElement("option");
  option.value = s.name;
  option.textContent = s.name;
  select_el.appendChild(option);
}});
select_el.addEventListener("change", (event) => select(event.target.value));

var RENDERERS = [
  ['table', renderTable],
  ['run info', renderRunInfo],
  ['plot', renderPlot]
];
// Rendered independently and run-info first. These draw with libraries loaded from a CDN,
// so with no network Plotly is undefined and renderPlot throws. Called in a bare sequence
// that exception aborted the rest of the script, and the block it took with it was Run --
// the seed, the checkpoint revision, and the line recording that nothing was screened.
// Losing the screening disclaimer because a chart library did not load is not acceptable,
// so each renderer is isolated and the provenance is drawn before the decoration.
for (var i = 0; i < RENDERERS.length; i++) {{
  try {{
    RENDERERS[i][1]();
  }} catch (err) {{
    console.error('workflow-037 report: ' + RENDERERS[i][0] + ' failed to render', err);
  }}
}}
if (DATA.structures.length) select(DATA.structures[0].name);
</script>
</body>
</html>
"""


def main():
    os.makedirs("./outputs", exist_ok=True)

    top_n = int_param("PARAM_TOP_N", "top_n", 10)
    if top_n < 1:
        fail(f"top_n must be at least 1, got {top_n}")

    metrics_path = "./inputs/design_metrics.json"
    if not os.path.exists(metrics_path):
        fail(f"Metrics not found: {metrics_path}")
    try:
        with open(metrics_path) as handle:
            metrics = json.load(handle)
    except Exception as exc:
        fail(f"Could not read '{metrics_path}': {exc}")

    designs = sorted(glob.glob("./inputs/*_generated.cif"))
    if not designs:
        fail("No generated structures found in ./inputs/ (*_generated.cif)")

    rows = metrics.get("designs", [])
    by_name = {row["design"]: row for row in rows}

    def order_key(path):
        row = by_name.get(os.path.basename(path), {})
        # A clashing design has backbone passing through the target, so it scores the
        # *highest* raw contact count there is — ordering on contacts alone puts the
        # least physical design at the top of the report and the first structure in the
        # viewer. Sink them instead: they are still shown, still counted, and flagged in
        # their own column, but they stop being what the reader sees first.
        return (
            bool(row.get("steric_clash")),
            -(row.get("n_interface_contacts") or 0),
            os.path.basename(path),
        )

    # Convenience ordering, not a ranking — see the module docstring.
    ordered = sorted(designs, key=order_key)
    metrics["designs"] = [
        by_name[os.path.basename(path)]
        for path in ordered
        if os.path.basename(path) in by_name
    ] or rows

    structures = []
    for path in ordered[:top_n]:
        with open(path) as handle:
            structures.append({"name": os.path.basename(path), "cif": handle.read()})

    html = build_html(metrics, structures, top_n)
    with open("./outputs/report.html", "w") as handle:
        handle.write(html)

    size_mb = round(os.path.getsize("./outputs/report.html") / 1e6, 2)
    print(
        f"Reporting {len(structures)} of {len(designs)} design(s); "
        f"wrote ./outputs/report.html ({size_mb} MB)",
        flush=True,
    )
    if size_mb > 50:
        print(
            "WARNING: report is large because structures are inlined; lower top_n.",
            flush=True,
        )


if __name__ == "__main__":
    main()
