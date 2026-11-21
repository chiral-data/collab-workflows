#!/usr/bin/env python3
"""Node 07 — Report.

One self-contained HTML report. All data is embedded inline; Bootstrap, Plotly and 3Dmol come
from CDN and there are no local assets.

Design notes, which are deliberate rather than stylistic:

* **The caveats lead.** A developability report that buries its interpretation limits gets
  discounted by anyone who knows the field, and rightly so. The flag legend sits next to the
  statement that these thresholds describe where clinical-stage antibodies sit rather than
  predicting success, and the calibration offset is stated where the flags are, not in a
  footnote.
* **Flags are never colour alone.** Every flag chip carries an icon and a word as well as a
  colour, which is also what the light-mode contrast of the amber status step requires: the
  status palette's warning step sits below 3:1 on a light surface by design, and the icon+label
  pairing plus the full table view is the documented mitigation.
* **Model error travels with the metric.** Every structure-derived number is shown beside the
  per-residue error of the model it came from, because CDR-H3 is both the loop that matters
  most and the one modelled worst.
* **No dual-axis charts, one measure per plot.** PSH gets its own plot with the published
  threshold bands and the clinical-stage reference mean drawn on it, so a value can be read
  against the distribution it is being judged by.

The viewer shows the shortlist with CDRs highlighted and flagged liability residues called out,
so a red PSH flag can be looked at rather than taken on faith.
"""
import csv
import glob
import json
import os
import sys

BOOTSTRAP_CSS = "https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css"
PLOTLY_JS = "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"
THREEDMOL_JS = "https://3Dmol.org/build/3Dmol-min.js"

# Status palette (fixed, never themed) + categorical slot 1, with the chart chrome inks.
PALETTE = {
    "good": "#0ca30c",
    "warning": "#fab219",
    "critical": "#d03b3b",
    "series_1_light": "#2a78d6",
    "series_1_dark": "#3987e5",
    "surface_light": "#fcfcfb",
    "surface_dark": "#1a1a19",
    "ink_light": "#0b0b0b",
    "ink_dark": "#f0efec",
    "muted": "#898781",
    "grid_light": "#e1e0d9",
    "grid_dark": "#2c2c2a",
}

METRICS = [
    ("flag_total_cdr_length", "Total CDR length"),
    ("flag_psh_cdr_vicinity", "Surface hydrophobicity (PSH)"),
    ("flag_ppc_cdr_vicinity", "Positive charge patches (PPC)"),
    ("flag_pnc_cdr_vicinity", "Negative charge patches (PNC)"),
    ("flag_sfvcsp", "Fv charge symmetry (SFvCSP)"),
]


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def warn(message):
    print(f"WARNING: {message}", flush=True)


def read_csv(path, label, required=True):
    if not os.path.exists(path):
        if required:
            fail(f"{label} not found: {path}. The upstream node must run first.")
        warn(f"{label} not found ({path}); the report will omit that section.")
        return []
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path, label, required=False):
    if not os.path.exists(path):
        if required:
            fail(f"{label} not found: {path}.")
        return {}
    try:
        with open(path) as handle:
            return json.load(handle)
    except Exception as exc:
        warn(f"Could not read {label} ('{path}'): {exc}")
        return {}


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def main():
    os.makedirs("./outputs", exist_ok=True)
    top_n = int(os.environ.get("PARAM_TOP_N", "10"))

    panel = read_csv("./inputs/ranked_panel.csv", "ranked_panel.csv")
    summary = read_json("./inputs/panel_summary.json", "panel_summary.json", required=True)
    developability = read_json("./inputs/developability.json", "developability.json")
    humanness = read_json("./inputs/humanness.json", "humanness.json")
    proposals = read_csv("./inputs/proposals.csv", "proposals.csv", required=False)
    scan_manifest = read_json("./inputs/scan_manifest.json", "scan_manifest.json")

    if not panel:
        fail("ranked_panel.csv contains no rows; node 06 must run first.")

    shortlist = [row for row in panel if row.get("shortlist") == "yes"][:top_n]
    structures = {}
    for row in shortlist:
        path = f"./inputs/models/{row['clone_id']}.pdb"
        if os.path.exists(path):
            with open(path) as handle:
                structures[row["clone_id"]] = handle.read()
    if not structures and glob.glob("./inputs/models/*.pdb"):
        warn("No shortlisted clone had a model file; the 3D viewer will be empty.")

    # Per-clone CDR and flagged-residue positions for the viewer.
    viewer_annotations = {}
    for clone_id in structures:
        entry = (developability.get("clones") or {}).get(clone_id) or {}
        cdr_residues, liability_residues = [], []
        for residue in entry.get("per_residue") or []:
            if str(residue.get("region", "")).startswith("CDR"):
                cdr_residues.append(
                    {"chain": residue["chain"], "position": residue["position"]}
                )
        for hit in entry.get("exposed_liabilities") or []:
            if hit.get("exposed"):
                liability_residues.append(
                    {
                        "chain": "H" if hit.get("chain") == "VH" else "L",
                        "position": hit.get("position"),
                        "motif": hit.get("motif"),
                    }
                )
        viewer_annotations[clone_id] = {
            "cdr": cdr_residues,
            "liabilities": liability_residues,
        }

    flag_counts = {
        column: {
            state: sum(1 for row in panel if row.get(column) == state)
            for state in ("green", "amber", "red")
        }
        for column, _ in METRICS
    }

    psh_rows = sorted(
        (
            {
                "clone_id": row["clone_id"],
                "psh": to_float(row.get("psh_cdr_vicinity")),
                "flag": row.get("flag_psh_cdr_vicinity") or "",
                "error": to_float(row.get("cdrh3_mean_error_angstrom")),
            }
            for row in panel
            if to_float(row.get("psh_cdr_vicinity")) is not None
        ),
        key=lambda item: item["psh"],
    )
    humanness_rows = sorted(
        (
            {
                "clone_id": row["clone_id"],
                "content": to_float(row.get("mean_peptide_content")),
                "identity": to_float(row.get("mean_v_identity")),
            }
            for row in panel
            if to_float(row.get("mean_peptide_content")) is not None
        ),
        key=lambda item: item["content"],
        reverse=True,
    )

    calibration = (developability.get("calibration") or {})
    implementation = (developability.get("implementation") or {})
    cst = calibration.get("cst_reference") or {}
    kept_proposals = [row for row in proposals if row.get("verdict") == "kept"]
    dropped_proposals = [row for row in proposals if row.get("verdict") == "dropped"]
    scan_enabled = bool(scan_manifest.get("enabled"))

    payload = {
        "panel": panel,
        "summary": summary,
        "flagCounts": flag_counts,
        "metrics": [{"column": column, "label": label} for column, label in METRICS],
        "pshRows": psh_rows,
        "humannessRows": humanness_rows,
        "structures": structures,
        "viewerAnnotations": viewer_annotations,
        "proposals": kept_proposals[:200],
        "droppedProposals": dropped_proposals[:200],
        "scanEnabled": scan_enabled,
        "scanManifest": scan_manifest,
        "calibration": calibration,
        "implementation": implementation,
        "cst": cst,
        "humannessMeta": {
            "definition": humanness.get("score_definition", ""),
            "caveat": humanness.get("caveat", ""),
            "settings": humanness.get("settings", {}),
        },
        "palette": PALETTE,
        "thresholds": {
            "total_cdr_length": {"amber": [[54, 60]], "red_above": 60, "red_below": None},
            "psh_cdr_vicinity": {
                "amber": [[83.84, 100.71], [156.20, 173.85]],
                "red_above": 173.85,
                "red_below": 83.84,
            },
            "ppc_cdr_vicinity": {"amber": [[1.25, 3.16]], "red_above": 3.16, "red_below": None},
            "pnc_cdr_vicinity": {"amber": [[1.84, 3.50]], "red_above": 3.50, "red_below": None},
            "sfvcsp": {"amber": [[-20.40, -6.30]], "red_above": None, "red_below": -20.40},
        },
    }

    html = build_html(payload, top_n)
    with open("./outputs/report.html", "w") as handle:
        handle.write(html)

    size_kb = round(os.path.getsize("./outputs/report.html") / 1024, 1)
    print(
        f"Wrote report.html ({size_kb} KB) covering {len(panel)} clone(s), "
        f"{len(structures)} embedded structure(s), {len(kept_proposals)} proposal(s)",
        flush=True,
    )


def build_html(payload, top_n):
    data_json = json.dumps(payload)
    palette = payload["palette"]
    calibration = payload["calibration"]
    cst = payload["cst"]
    summary = payload["summary"]
    offset = calibration.get("panel_psh_offset") or ""
    deviations = calibration.get("known_deviations") or []
    deviation_items = "\n".join(f"<li>{escape(item)}</li>" for item in deviations)
    policy_items = "\n".join(
        f"<li>{escape(item)}</li>" for item in (summary.get("ranking_policy") or [])
    )
    screen_items = "\n".join(
        f"<li>{escape(item)}</li>"
        for item in ((payload.get("scanManifest") or {}).get("screen_rules") or [])
    )

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Antibody Developability Triage</title>
<link href="{BOOTSTRAP_CSS}" rel="stylesheet">
<script src="{PLOTLY_JS}"></script>
<script src="{THREEDMOL_JS}"></script>
<style>
  :root {{
    color-scheme: light;
    --surface-1: {palette['surface_light']};
    --page: #f9f9f7;
    --ink: {palette['ink_light']};
    --ink-2: #52514e;
    --muted: {palette['muted']};
    --grid: {palette['grid_light']};
    --series-1: {palette['series_1_light']};
    --good: {palette['good']};
    --warning: {palette['warning']};
    --critical: {palette['critical']};
    --border: rgba(11,11,11,0.10);
  }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{
      color-scheme: dark;
      --surface-1: {palette['surface_dark']};
      --page: #0d0d0d;
      --ink: {palette['ink_dark']};
      --ink-2: #c3c2b7;
      --grid: {palette['grid_dark']};
      --series-1: {palette['series_1_dark']};
      --border: rgba(255,255,255,0.10);
    }}
  }}
  :root[data-theme="dark"] {{
    color-scheme: dark;
    --surface-1: {palette['surface_dark']};
    --page: #0d0d0d;
    --ink: {palette['ink_dark']};
    --ink-2: #c3c2b7;
    --grid: {palette['grid_dark']};
    --series-1: {palette['series_1_dark']};
    --border: rgba(255,255,255,0.10);
  }}
  body {{
    background: var(--page);
    color: var(--ink);
    font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
  }}
  .card-plain {{
    background: var(--surface-1);
    border: 1px solid var(--border);
    border-radius: 10px;
  }}
  .tile-value {{ font-size: 2rem; font-weight: 600; line-height: 1.1; }}
  .tile-label {{ color: var(--ink-2); font-size: .85rem; }}
  .chip {{
    display: inline-flex; align-items: center; gap: .3rem;
    padding: .1rem .45rem; border-radius: 999px;
    font-size: .78rem; font-weight: 600; white-space: nowrap;
  }}
  .chip-good {{ background: color-mix(in srgb, var(--good) 18%, transparent); color: var(--ink); }}
  .chip-amber {{ background: color-mix(in srgb, var(--warning) 26%, transparent); color: var(--ink); }}
  .chip-red {{ background: color-mix(in srgb, var(--critical) 20%, transparent); color: var(--ink); }}
  .chip-none {{ background: color-mix(in srgb, var(--muted) 18%, transparent); color: var(--ink-2); }}
  .dot {{ width: .55rem; height: .55rem; border-radius: 50%; display: inline-block; }}
  table.table {{ color: var(--ink); font-variant-numeric: tabular-nums; }}
  .table > :not(caption) > * > * {{ background: transparent; }}
  .caveat {{ border-left: 3px solid var(--warning); }}
  .seq {{ font-family: ui-monospace, SFMono-Regular, Menlo, monospace; font-size: .8rem; word-break: break-all; }}
  .viewer {{ width: 100%; height: 460px; position: relative; border-radius: 8px; overflow: hidden; }}
  h2 {{ font-size: 1.25rem; margin-top: 2rem; }}
  .muted {{ color: var(--ink-2); }}
  .small-muted {{ color: var(--muted); font-size: .82rem; }}
</style>
</head>
<body>
<div class="container-xxl py-4" style="max-width: 1320px;">

  <header class="mb-3">
    <h1 class="h3 mb-1">Antibody Developability Triage</h1>
    <p class="muted mb-0">
      TAP-style developability metrics, sequence liabilities, humanness and mutation proposals
      across {summary.get('n_clones', 0)} clone(s). Workflow 039.
    </p>
  </header>

  <div class="card-plain caveat p-3 mb-4">
    <strong>How to read this report</strong>
    <ul class="mb-2 mt-2">
      <li><strong>Flags are risk indicators, not predictions of clinical success.</strong>
          Amber means the value sits in the extremes of the clinical-stage therapeutic
          distribution; red means a value not previously observed among those therapeutics.
          Several approved antibodies violate these guidelines.</li>
      <li><strong>These are TAP-<em>style</em> metrics</strong>, reimplemented from Raybould et
          al. (PNAS 2019) because TAP ships as a web service with no redistributable source.
          They carry a calibration offset and will not match the OPIG web server exactly.</li>
      <li><strong>Humanness is not immunogenicity.</strong> No MHC-II binding or T-cell epitope
          prediction is performed anywhere in this workflow.</li>
      <li><strong>No antigen was considered.</strong> Nothing here speaks to potency,
          specificity or epitope coverage, and the mutation proposals are not affinity
          predictions.</li>
      <li><strong>Structures are predictions.</strong> Every structure-derived metric inherits
          the model's error, and CDR-H3 carries most of it. Model error is shown beside the
          metrics for that reason.</li>
    </ul>
    <div class="small-muted">{escape(calibration.get('interpretation', ''))}</div>
  </div>

  <div class="row g-3" id="tiles"></div>

  <h2>Flag overview</h2>
  <p class="small-muted mb-2">
    Counts per metric across the panel. Thresholds from Table 2 of Raybould et al. (2019),
    derived from 242 clinical-stage therapeutics.
  </p>
  <div class="card-plain p-3 mb-2"><div id="flagChart" style="height: 300px;"></div></div>

  <h2>Surface hydrophobicity against the clinical-stage distribution</h2>
  <p class="small-muted mb-2">
    One point per clone, sorted. Shaded bands are the published amber regions; beyond them is
    red. The dashed line is the clinical-stage mean
    ({cst.get('psh_cdr_vicinity_mean', 'n/a')} &plusmn; {cst.get('psh_cdr_vicinity_sd', 'n/a')},
    {cst.get('n_models', 'n/a')} models) — a panel sitting far from it is a sign of a
    calibration problem in this implementation rather than an unusual panel.
  </p>
  <div class="card-plain p-3 mb-2"><div id="pshChart" style="height: 360px;"></div></div>
  <p class="small-muted">{escape(offset)}</p>

  <h2>Humanness</h2>
  <p class="small-muted mb-2">
    {escape(payload['humannessMeta'].get('definition', ''))}
    Higher is more human. <strong>This is a proxy for immunogenicity risk, not a prediction
    of it.</strong>
  </p>
  <div class="card-plain p-3 mb-2"><div id="humannessChart" style="height: 340px;"></div></div>

  <h2>Ranked panel</h2>
  <p class="small-muted mb-2">
    Nothing is filtered: every clone that entered the panel is here. The ranking is a sort over
    the columns shown, in this order:
  </p>
  <ol class="small-muted">{policy_items}</ol>
  <div class="card-plain p-2 mb-2 table-responsive">
    <table class="table table-sm align-middle mb-0" id="panelTable"></table>
  </div>

  <h2>Shortlist structures</h2>
  <p class="small-muted mb-2">
    CDRs in colour, exposed liability residues shown as spheres. The B-factor column of these
    models holds ImmuneBuilder's per-residue error estimate, not a crystallographic B-factor.
  </p>
  <div class="card-plain p-3 mb-2">
    <div class="d-flex flex-wrap gap-2 align-items-center mb-2">
      <label class="small-muted mb-0" for="cloneSelect">Clone</label>
      <select class="form-select form-select-sm" id="cloneSelect" style="max-width: 260px;"></select>
      <span class="small-muted" id="viewerNote"></span>
    </div>
    <div class="viewer" id="viewer"></div>
  </div>

  <h2>Mutation proposals</h2>
  <div id="proposalSection"></div>

  <h2>Implementation and provenance</h2>
  <div class="card-plain p-3 mb-4">
    <p class="mb-2"><strong>Known deviations from TAP's own implementation</strong></p>
    <ul class="small-muted">{deviation_items}</ul>
    <p class="mb-1"><strong>Re-screening rules applied to every mutation proposal</strong></p>
    <ul class="small-muted">{screen_items or '<li>Mutation scan not run.</li>'}</ul>
    <p class="mb-1"><strong>Metric implementation</strong></p>
    <pre class="small-muted mb-0" style="white-space: pre-wrap;">{escape(json.dumps(payload['implementation'], indent=2))}</pre>
  </div>

  <footer class="small-muted pb-4">
    Threshold source: {escape(calibration.get('threshold_source', 'n/a'))}
    Generated by workflow-039 (collab-workflows).
  </footer>
</div>

<script id="report-data" type="application/json">{data_json}</script>
<script>
const DATA = JSON.parse(document.getElementById('report-data').textContent);
const P = DATA.palette;

function ink() {{
  return getComputedStyle(document.body).getPropertyValue('color').trim();
}}
function cssVar(name) {{
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}}

const FLAG_META = {{
  green: {{ label: 'pass', icon: '\\u2713', cls: 'chip-good', color: P.good }},
  amber: {{ label: 'amber', icon: '\\u0021', cls: 'chip-amber', color: P.warning }},
  red: {{ label: 'red', icon: '\\u2715', cls: 'chip-red', color: P.critical }},
  '': {{ label: 'n/a', icon: '\\u2013', cls: 'chip-none', color: P.muted }}
}};

function chip(state) {{
  const meta = FLAG_META[state] || FLAG_META[''];
  return `<span class="chip ${{meta.cls}}"><span aria-hidden="true">${{meta.icon}}</span>${{meta.label}}</span>`;
}}

/* Stat tiles. A tile, not a chart: these are single headline numbers. */
function renderTiles() {{
  const s = DATA.summary;
  const pc = s.peptide_content;
  const tiles = [
    ['Clones', s.n_clones, 'in the panel, nothing filtered'],
    ['Shortlisted', s.n_shortlisted, `top ${{s.top_n}} after ranking`],
    ['No flags, no exposed liabilities', (s.clean_clones || []).length, 'of the scorable clones'],
    ['Red flags', s.n_red_total, 'values unobserved among clinical-stage therapeutics'],
    ['Amber flags', s.n_amber_total, 'values in the extremes of that distribution'],
    ['Mean human peptide content', pc ? pc.mean.toFixed(3) : 'n/a',
     pc ? `range ${{pc.min.toFixed(3)}} to ${{pc.max.toFixed(3)}}` : 'not computed']
  ];
  document.getElementById('tiles').innerHTML = tiles.map(([label, value, note]) => `
    <div class="col-6 col-lg-2">
      <div class="card-plain p-3 h-100">
        <div class="tile-value">${{value}}</div>
        <div class="tile-label">${{label}}</div>
        <div class="small-muted mt-1">${{note}}</div>
      </div>
    </div>`).join('');
}}

const LAYOUT_BASE = () => ({{
  paper_bgcolor: 'rgba(0,0,0,0)',
  plot_bgcolor: 'rgba(0,0,0,0)',
  font: {{ color: ink(), family: 'system-ui, -apple-system, "Segoe UI", sans-serif', size: 12 }},
  margin: {{ l: 160, r: 24, t: 10, b: 44 }},
  hoverlabel: {{ font: {{ family: 'system-ui, sans-serif' }} }},
  xaxis: {{ gridcolor: cssVar('--grid'), zeroline: false, linecolor: cssVar('--grid') }},
  yaxis: {{ gridcolor: cssVar('--grid'), zeroline: false, linecolor: cssVar('--grid') }}
}});

/* Flag overview: three status series, so a legend is present and each bar is directly
   labelled with its count. Status colour never carries the meaning alone. */
function renderFlagChart() {{
  const labels = DATA.metrics.map(m => m.label);
  const states = [['green', 'pass'], ['amber', 'amber'], ['red', 'red']];
  const traces = states.map(([state, name]) => ({{
    type: 'bar',
    orientation: 'h',
    name: `${{FLAG_META[state].icon}} ${{name}}`,
    y: labels,
    x: DATA.metrics.map(m => DATA.flagCounts[m.column][state]),
    marker: {{ color: FLAG_META[state].color, line: {{ color: cssVar('--surface-1'), width: 2 }} }},
    text: DATA.metrics.map(m => DATA.flagCounts[m.column][state] || ''),
    textposition: 'inside',
    insidetextfont: {{ color: state === 'amber' ? '#0b0b0b' : '#ffffff' }},
    hovertemplate: '%{{y}}<br>' + name + ': %{{x}} clone(s)<extra></extra>'
  }}));
  const layout = Object.assign(LAYOUT_BASE(), {{
    barmode: 'stack',
    showlegend: true,
    legend: {{ orientation: 'h', y: -0.18 }},
    xaxis: Object.assign(LAYOUT_BASE().xaxis, {{ title: 'clones' }})
  }});
  Plotly.newPlot('flagChart', traces, layout, {{ displayModeBar: false, responsive: true }});
}}

/* PSH: one measure, one axis. Threshold bands are recessive background shapes; the
   clinical-stage mean is a dashed reference line. */
function renderPshChart() {{
  const rows = DATA.pshRows;
  if (!rows.length) {{ document.getElementById('pshChart').innerHTML =
    '<p class="small-muted">No PSH values available.</p>'; return; }}
  const t = DATA.thresholds.psh_cdr_vicinity;
  const shapes = t.amber.map(([lo, hi]) => ({{
    type: 'rect', xref: 'x', yref: 'paper', x0: lo, x1: hi, y0: 0, y1: 1,
    fillcolor: P.warning, opacity: 0.12, line: {{ width: 0 }}, layer: 'below'
  }}));
  shapes.push({{
    type: 'rect', xref: 'x', yref: 'paper', x0: t.red_above, x1: Math.max(t.red_above * 1.4,
      Math.max.apply(null, rows.map(r => r.psh)) + 10), y0: 0, y1: 1,
    fillcolor: P.critical, opacity: 0.10, line: {{ width: 0 }}, layer: 'below'
  }});
  shapes.push({{
    type: 'rect', xref: 'x', yref: 'paper', x0: 0, x1: t.red_below, y0: 0, y1: 1,
    fillcolor: P.critical, opacity: 0.10, line: {{ width: 0 }}, layer: 'below'
  }});
  if (DATA.cst && DATA.cst.psh_cdr_vicinity_mean) {{
    shapes.push({{
      type: 'line', xref: 'x', yref: 'paper',
      x0: DATA.cst.psh_cdr_vicinity_mean, x1: DATA.cst.psh_cdr_vicinity_mean,
      y0: 0, y1: 1, line: {{ color: cssVar('--muted'), width: 2, dash: 'dash' }}
    }});
  }}
  const trace = {{
    type: 'scatter', mode: 'markers',
    x: rows.map(r => r.psh),
    y: rows.map(r => r.clone_id),
    marker: {{
      size: 11, color: cssVar('--series-1'),
      line: {{ color: cssVar('--surface-1'), width: 2 }}
    }},
    customdata: rows.map(r => [r.flag || 'unflagged', r.error === null ? 'n/a' : r.error]),
    hovertemplate: '%{{y}}<br>PSH %{{x:.1f}} (%{{customdata[0]}})' +
      '<br>CDR-H3 model error %{{customdata[1]}} \\u00c5<extra></extra>'
  }};
  const layout = Object.assign(LAYOUT_BASE(), {{
    showlegend: false,
    height: Math.max(300, 26 * rows.length + 90),
    xaxis: Object.assign(LAYOUT_BASE().xaxis, {{ title: 'CDR-vicinity PSH (TAP-style)' }}),
    yaxis: Object.assign(LAYOUT_BASE().yaxis, {{ automargin: true, type: 'category' }})
  }});
  Plotly.newPlot('pshChart', [trace], layout, {{ displayModeBar: false, responsive: true }});
}}

/* Humanness: one series, so no legend box — the title names it. */
function renderHumannessChart() {{
  const rows = DATA.humannessRows;
  if (!rows.length) {{ document.getElementById('humannessChart').innerHTML =
    '<p class="small-muted">No humanness scores available.</p>'; return; }}
  const trace = {{
    type: 'bar', orientation: 'h',
    x: rows.map(r => r.content),
    y: rows.map(r => r.clone_id),
    marker: {{ color: cssVar('--series-1'), line: {{ color: cssVar('--surface-1'), width: 2 }} }},
    text: rows.map(r => r.content.toFixed(3)),
    textposition: 'outside',
    customdata: rows.map(r => r.identity === null ? 'n/a' : r.identity.toFixed(3)),
    hovertemplate: '%{{y}}<br>human peptide content %{{x:.3f}}' +
      '<br>mean germline V identity %{{customdata}}<extra></extra>'
  }};
  const layout = Object.assign(LAYOUT_BASE(), {{
    showlegend: false,
    height: Math.max(280, 26 * rows.length + 90),
    xaxis: Object.assign(LAYOUT_BASE().xaxis, {{
      title: 'fraction of 9-mers found in the human reference', range: [0, 1.08]
    }}),
    yaxis: Object.assign(LAYOUT_BASE().yaxis, {{ automargin: true, type: 'category' }})
  }});
  Plotly.newPlot('humannessChart', [trace], layout, {{ displayModeBar: false, responsive: true }});
}}

/* The table view, which is also the relief for the amber status step's light-mode contrast. */
function renderTable() {{
  const columns = [
    ['rank', 'Rank'], ['clone_id', 'Clone'], ['status', 'Status'],
    ['total_cdr_length', 'CDR len'], ['psh_cdr_vicinity', 'PSH'],
    ['ppc_cdr_vicinity', 'PPC'], ['pnc_cdr_vicinity', 'PNC'], ['sfvcsp', 'SFvCSP'],
    ['mean_peptide_content', 'Human'], ['n_exposed_high_liabilities', 'Exposed high liab.'],
    ['cdrh3_mean_error_angstrom', 'CDR-H3 err (\\u00c5)'], ['n_proposals', 'Proposals']
  ];
  const flagColumns = DATA.metrics.map(m => m.column);
  const head = '<thead><tr>' + columns.map(c => `<th>${{c[1]}}</th>`).join('') +
    DATA.metrics.map(m => `<th>${{m.label.replace(/ \\(.*\\)/, '')}}</th>`).join('') +
    '<th>Note</th></tr></thead>';
  const body = '<tbody>' + DATA.panel.map(row => {{
    const shortlisted = row.shortlist === 'yes';
    const cells = columns.map(([key]) => {{
      let value = row[key] === undefined || row[key] === '' ? '\\u2013' : row[key];
      if (key === 'clone_id' && shortlisted) value = `<strong>${{value}}</strong>`;
      return `<td>${{value}}</td>`;
    }}).join('');
    const flags = flagColumns.map(column => `<td>${{chip(row[column] || '')}}</td>`).join('');
    const note = row.thresholds_calibrated === 'no'
      ? '<span class="small-muted">nanobody: thresholds not calibrated</span>'
      : `<span class="small-muted">${{row.reason || ''}}</span>`;
    return `<tr>${{cells}}${{flags}}<td>${{note}}</td></tr>`;
  }}).join('') + '</tbody>';
  document.getElementById('panelTable').innerHTML = head + body;
}}

function renderProposals() {{
  const container = document.getElementById('proposalSection');
  if (!DATA.scanEnabled) {{
    container.innerHTML =
      '<p class="small-muted">The mutation scan was not run (mutate_scan=false).</p>';
    return;
  }}
  const kept = DATA.proposals;
  if (!kept.length) {{
    container.innerHTML = `<div class="card-plain p-3"><p class="mb-1">No proposal survived
      re-screening.</p><p class="small-muted mb-0">That is the screen working rather than a
      failure: ${{DATA.droppedProposals.length}} candidate(s) were dropped for introducing a
      liability or deepening a flagged risk.</p></div>`;
    return;
  }}
  const rows = kept.map(p => `<tr>
      <td>${{p.clone_id}}</td><td>${{p.chain}}</td><td>${{p.position}}</td>
      <td>${{p.region}}</td>
      <td class="seq">${{p.wild_type}}\\u2192${{p.mutant}}</td>
      <td>${{Number(p.probability_gain).toFixed(3)}}</td>
    </tr>`).join('');
  const dropped = DATA.droppedProposals.slice(0, 40).map(p => `<tr>
      <td>${{p.clone_id}}</td><td>${{p.chain}}</td><td>${{p.position}}</td>
      <td class="seq">${{p.wild_type}}\\u2192${{p.mutant}}</td>
      <td class="small-muted">${{p.reason}}</td>
    </tr>`).join('');
  container.innerHTML = `
    <p class="small-muted">
      AbLang2 probabilities, re-screened. <strong>Likelihood-based, not affinity-based</strong>:
      a high-probability mutant looks like a natural antibody, which is not the same as one that
      binds better, and no antigen was considered. Bench hypotheses only.
    </p>
    <div class="card-plain p-2 mb-3 table-responsive">
      <table class="table table-sm mb-0">
        <thead><tr><th>Clone</th><th>Chain</th><th>Position</th><th>Region</th>
          <th>Substitution</th><th>Probability gain</th></tr></thead>
        <tbody>${{rows}}</tbody>
      </table>
    </div>
    <details class="card-plain p-2 mb-2">
      <summary class="small-muted">Dropped by re-screening
        (${{DATA.droppedProposals.length}} shown up to 40)</summary>
      <div class="table-responsive mt-2">
        <table class="table table-sm mb-0">
          <thead><tr><th>Clone</th><th>Chain</th><th>Position</th><th>Substitution</th>
            <th>Reason</th></tr></thead>
          <tbody>${{dropped}}</tbody>
        </table>
      </div>
    </details>`;
}}

let viewer = null;
function renderViewer() {{
  const select = document.getElementById('cloneSelect');
  const ids = Object.keys(DATA.structures);
  if (!ids.length) {{
    document.getElementById('viewer').innerHTML =
      '<p class="small-muted p-3">No structures were available for the shortlist.</p>';
    select.disabled = true;
    return;
  }}
  select.innerHTML = ids.map(id => `<option value="${{id}}">${{id}}</option>`).join('');
  select.addEventListener('change', () => showStructure(select.value));
  showStructure(ids[0]);
}}

function showStructure(cloneId) {{
  const element = document.getElementById('viewer');
  element.innerHTML = '';
  viewer = $3Dmol.createViewer(element, {{ backgroundColor: 'rgba(0,0,0,0)' }});
  viewer.addModel(DATA.structures[cloneId], 'pdb');
  viewer.setStyle({{}}, {{ cartoon: {{ color: cssVar('--muted'), opacity: 0.75 }} }});

  const annotation = DATA.viewerAnnotations[cloneId] || {{ cdr: [], liabilities: [] }};
  const byChain = {{}};
  annotation.cdr.forEach(residue => {{
    const number = parseInt(residue.position, 10);
    if (!isNaN(number)) {{
      byChain[residue.chain] = byChain[residue.chain] || [];
      byChain[residue.chain].push(number);
    }}
  }});
  Object.keys(byChain).forEach((chain, index) => {{
    viewer.setStyle(
      {{ chain: chain, resi: byChain[chain] }},
      {{ cartoon: {{ color: index === 0 ? cssVar('--series-1') : P.good }} }}
    );
  }});

  annotation.liabilities.forEach(hit => {{
    const number = parseInt(hit.position, 10);
    if (isNaN(number)) return;
    viewer.addStyle(
      {{ chain: hit.chain, resi: [number] }},
      {{ stick: {{ color: P.critical, radius: 0.22 }},
        sphere: {{ color: P.critical, radius: 0.55 }} }}
    );
  }});

  viewer.zoomTo();
  viewer.render();

  const note = annotation.liabilities.length
    ? `${{annotation.liabilities.length}} exposed liability residue(s) shown in red`
    : 'no exposed liability residues in this clone';
  document.getElementById('viewerNote').textContent =
    `CDRs highlighted; ${{note}}.`;
}}

renderTiles();
renderFlagChart();
renderPshChart();
renderHumannessChart();
renderTable();
renderProposals();
renderViewer();
</script>
</body>
</html>
"""


def escape(value):
    return (
        str(value)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


if __name__ == "__main__":
    main()
