#!/usr/bin/env python3
"""Node 06 — Report.

Build ./outputs/report.html: one self-contained file, all libraries from CDN, all data
embedded inline as JavaScript. No external data files, so the report can be copied anywhere.

Chart choices (the data's job picks the form):
  * Rg per frame            — change across frames, 2 px line + >=8 px markers, one series per
                              ensemble (categorical slots 1 and 2, validated for both modes)
  * Rg distribution         — magnitude/shape, histogram, single series
  * Pairwise C-alpha RMSD   — continuous magnitude, heatmap on ONE sequential blue hue
                              (never a rainbow; the lightest step means "identical")
  * RMSF along the chain    — one series per chain per ensemble; per residue when the
                              selection is C-alphas, per atom when node 05 fell back to every
                              atom for a system with no protein chain, and labelled as such
  * Physicality             — state, so the fixed status palette with an icon and a label;
                              never colour alone

A legend is present whenever two ensembles are plotted, and every frame is also listed in the
table below, so identity is never carried by colour alone.

The report leads with the caveats rather than burying them: BioEmu2 frequencies are not
calibrated Boltzmann probabilities and the sampled frames are independent draws, not a
trajectory. A reader who takes a cluster size for a population has been misled by the report,
not by the model.
"""
import csv
import json
import os
import sys

PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.2.min.js"
BOOTSTRAP_CDN = "https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css"
THREEDMOL_CDN = "https://3Dmol.org/build/3Dmol-min.js"


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def read_json(path, label, default=None):
    if not os.path.exists(path):
        if default is not None:
            return default
        fail(f"Required input not found: {path}")
    try:
        with open(path) as handle:
            return json.load(handle)
    except Exception as exc:
        fail(f"Could not read {label}: {exc}")


def read_csv(path):
    if not os.path.exists(path):
        fail(f"Required input not found: {path}")
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def read_text(path, label):
    if not os.path.exists(path):
        print(f"WARNING: {label} not found at {path}; the 3D viewer will be omitted.",
              flush=True)
        return ""
    with open(path) as handle:
        return handle.read()


def number(value, digits=3):
    try:
        return round(float(value), digits)
    except (TypeError, ValueError):
        return None


def fmt(value, digits=3, dash="—"):
    rounded = number(value, digits)
    return dash if rounded is None else f"{rounded:g}"


HTML = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>BioEmu2 Ensemble Report</title>
<link href="{bootstrap}" rel="stylesheet">
<script src="{plotly}"></script>
<script src="{threedmol}"></script>
<style>
  :root {{
    color-scheme: light;
    --surface-1: #fcfcfb;
    --page-plane: #f9f9f7;
    --text-primary: #0b0b0b;
    --text-secondary: #52514e;
    --text-muted: #898781;
    --grid: #e1e0d9;
    --axis: #c3c2b7;
    --border: rgba(11,11,11,0.10);
    --series-1: #2a78d6;
    --series-2: #eb6834;
    --seq-100: #cde2fb;
    --seq-700: #0d366b;
    --status-good: #0ca30c;
    --status-warning: #fab219;
    --status-critical: #d03b3b;
  }}
  @media (prefers-color-scheme: dark) {{
    :root:not([data-theme="light"]) {{
      color-scheme: dark;
      --surface-1: #1a1a19;
      --page-plane: #0d0d0d;
      --text-primary: #f0efec;
      --text-secondary: #c3c2b7;
      --text-muted: #898781;
      --grid: #2c2c2a;
      --axis: #383835;
      --border: rgba(255,255,255,0.10);
      --series-1: #3987e5;
      --series-2: #d95926;
      --seq-100: #0d366b;
      --seq-700: #cde2fb;
    }}
  }}
  :root[data-theme="dark"] {{
    color-scheme: dark;
    --surface-1: #1a1a19;
    --page-plane: #0d0d0d;
    --text-primary: #f0efec;
    --text-secondary: #c3c2b7;
    --text-muted: #898781;
    --grid: #2c2c2a;
    --axis: #383835;
    --border: rgba(255,255,255,0.10);
    --series-1: #3987e5;
    --series-2: #d95926;
    --seq-100: #0d366b;
    --seq-700: #cde2fb;
  }}
  body {{
    background: var(--page-plane);
    color: var(--text-primary);
    font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
    padding: 0 16px 48px;
  }}
  .wrap {{ max-width: 1200px; margin: 0 auto; }}
  h1 {{ font-size: 1.6rem; font-weight: 600; margin: 32px 0 4px; }}
  h2 {{ font-size: 1.15rem; font-weight: 600; margin: 32px 0 12px; }}
  .sub {{ color: var(--text-secondary); margin-bottom: 24px; }}
  .card2 {{
    background: var(--surface-1);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 20px;
    margin-bottom: 20px;
  }}
  .tile-row {{ display: flex; flex-wrap: wrap; gap: 12px; }}
  .tile {{
    background: var(--surface-1);
    border: 1px solid var(--border);
    border-radius: 10px;
    padding: 14px 18px;
    flex: 1 1 160px;
  }}
  .tile .label {{
    color: var(--text-secondary); font-size: .78rem;
    text-transform: uppercase; letter-spacing: .04em;
  }}
  .tile .value {{ font-size: 1.55rem; font-weight: 600; margin-top: 2px; }}
  .tile .unit {{ color: var(--text-muted); font-size: .9rem; font-weight: 400; }}
  .caveat {{
    background: var(--surface-1);
    border: 1px solid var(--border);
    border-left: 4px solid var(--status-warning);
    border-radius: 8px; padding: 14px 18px; margin-bottom: 24px;
  }}
  .caveat ul {{ margin: 8px 0 0; padding-left: 20px; color: var(--text-secondary); }}
  table {{ color: var(--text-primary); font-variant-numeric: tabular-nums; }}
  table.table {{ --bs-table-bg: transparent; --bs-table-color: var(--text-primary); }}
  th {{ color: var(--text-secondary); font-weight: 600; font-size: .85rem; }}
  .pill {{
    display: inline-flex; align-items: center; gap: 5px;
    font-size: .8rem; font-weight: 600;
  }}
  .dot {{ width: 9px; height: 9px; border-radius: 50%; display: inline-block; }}
  .legend {{ display: flex; gap: 18px; margin-bottom: 8px; flex-wrap: wrap; }}
  .legend span {{ display: inline-flex; align-items: center; gap: 6px;
                  color: var(--text-secondary); font-size: .85rem; }}
  #viewer {{ width: 100%; height: 460px; position: relative; border-radius: 8px; }}
  .chart {{ width: 100%; height: 340px; }}
  .meta {{ color: var(--text-muted); font-size: .82rem; }}
  code {{ color: var(--text-secondary); }}
</style>
</head>
<body>
<div class="wrap">

<h1>BioEmu2 Ensemble Report</h1>
<div class="sub">{headline}</div>

<div class="caveat">
  <strong>How to read this report</strong>
  <ul>
    <li><strong>Frame counts are not populations.</strong> BioEmu2's sampling frequency is not
        a calibrated Boltzmann probability, so a cluster's size is a count of frames and
        never a free energy.</li>
    <li><strong>Sampled frames are independent draws, not a trajectory.</strong> Frame order
        carries no time information.{propagation_note}</li>
    <li><strong>The alignment conditions everything else.</strong> {msa_sentence}</li>
    <li><strong>Every number here is geometric.</strong> Nothing in this report is a
        confidence estimate, and BioEmu2 is research software, not a validated commercial
        tool.</li>
  </ul>
</div>

<div class="tile-row">
  {tiles}
</div>

<h2>Conformational spread</h2>
<div class="card2">
  <div class="legend">{legend}</div>
  <div id="rg-chart" class="chart"></div>
  <div class="meta">Radius of gyration per frame. A flat line means every frame is the same
    size, which for an equilibrium run is a sign the ensemble has collapsed.</div>
</div>

<div class="row">
  <div class="col-lg-6">
    <div class="card2">
      <div id="rg-hist" class="chart"></div>
      <div class="meta">Distribution of the radius of gyration across the sampled
        ensemble.</div>
    </div>
  </div>
  <div class="col-lg-6">
    <div class="card2">
      <div id="rmsd-heat" class="chart"></div>
      <div class="meta">Pairwise C-alpha RMSD between sampled frames. The palest cells are
        pairs of near-identical structures; blocks of pale cells are conformational
        clusters.</div>
    </div>
  </div>
</div>

<div class="card2">
  <div class="legend">{legend}</div>
  <div id="rmsf-chart" class="chart"></div>
  <div class="meta">{rmsf_caption}</div>
</div>

<h2>Representative structures</h2>
<div class="card2">
  {viewer_block}
</div>

<h2>Per-frame detail</h2>
<div class="card2">
  <div class="table-responsive">
    <table class="table table-sm align-middle">
      <thead><tr>
        <th>Frame</th><th>Ensemble</th><th>R<sub>g</sub> (nm)</th>
        <th>RMSD to medoid (nm)</th><th>Physicality</th><th>Reason</th>
      </tr></thead>
      <tbody>{rows}</tbody>
    </table>
  </div>
</div>

<h2>Run provenance</h2>
<div class="card2">
  <div class="table-responsive">
    <table class="table table-sm">
      <tbody>{provenance}</tbody>
    </table>
  </div>
</div>

</div>

<script>
const DATA = {data_json};

function ink() {{
  const css = getComputedStyle(document.documentElement);
  return {{
    text: css.getPropertyValue('--text-primary').trim(),
    secondary: css.getPropertyValue('--text-secondary').trim(),
    muted: css.getPropertyValue('--text-muted').trim(),
    grid: css.getPropertyValue('--grid').trim(),
    axis: css.getPropertyValue('--axis').trim(),
    s1: css.getPropertyValue('--series-1').trim(),
    s2: css.getPropertyValue('--series-2').trim(),
    seqLo: css.getPropertyValue('--seq-100').trim(),
    seqHi: css.getPropertyValue('--seq-700').trim()
  }};
}}

function layout(c, xTitle, yTitle) {{
  return {{
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: 'rgba(0,0,0,0)',
    font: {{ family: 'system-ui, -apple-system, "Segoe UI", sans-serif',
             color: c.secondary, size: 12 }},
    margin: {{ l: 62, r: 18, t: 12, b: 48 }},
    showlegend: false,
    hovermode: 'closest',
    xaxis: {{ title: {{ text: xTitle, font: {{ color: c.secondary }} }},
              gridcolor: c.grid, linecolor: c.axis, zeroline: false,
              tickfont: {{ color: c.muted }} }},
    yaxis: {{ title: {{ text: yTitle, font: {{ color: c.secondary }} }},
              gridcolor: c.grid, linecolor: c.axis, zeroline: false,
              tickfont: {{ color: c.muted }} }}
  }};
}}

const CONFIG = {{ displayModeBar: false, responsive: true }};

function draw() {{
  const c = ink();
  const colors = [c.s1, c.s2];

  // Rg per frame — one series per ensemble, 2px lines, 8px markers.
  const rgTraces = DATA.ensembles.map((e, i) => ({{
    x: e.frames, y: e.rg, name: e.label, mode: 'lines+markers',
    type: 'scatter',
    line: {{ color: colors[i % 2], width: 2 }},
    marker: {{ color: colors[i % 2], size: 8,
               line: {{ color: 'rgba(0,0,0,0)', width: 0 }} }},
    hovertemplate: e.label + ' frame %{{x}}<br>R<sub>g</sub> %{{y:.3f}} nm<extra></extra>'
  }}));
  Plotly.newPlot('rg-chart', rgTraces, layout(c, 'Frame', 'Radius of gyration (nm)'), CONFIG);

  // Rg distribution — single series, so no legend box.
  Plotly.newPlot('rg-hist', [{{
    x: DATA.ensembles[0].rg, type: 'histogram',
    marker: {{ color: c.s1, line: {{ color: 'rgba(0,0,0,0)', width: 0 }} }},
    hovertemplate: 'R<sub>g</sub> %{{x:.3f}} nm<br>%{{y}} frame(s)<extra></extra>'
  }}], layout(c, 'Radius of gyration (nm)', 'Frames'), CONFIG);

  // Pairwise RMSD — ONE sequential hue, light (identical) to dark (different).
  if (DATA.rmsdMatrix && DATA.rmsdMatrix.length) {{
    Plotly.newPlot('rmsd-heat', [{{
      z: DATA.rmsdMatrix, type: 'heatmap',
      colorscale: [[0, c.seqLo], [1, c.seqHi]],
      colorbar: {{ title: {{ text: 'nm', font: {{ color: c.secondary }} }},
                   tickfont: {{ color: c.muted }}, thickness: 12, outlinewidth: 0 }},
      hovertemplate: 'frame %{{x}} vs %{{y}}<br>%{{z:.3f}} nm<extra></extra>'
    }}], layout(c, 'Frame', 'Frame'), CONFIG);
  }} else {{
    document.getElementById('rmsd-heat').innerHTML =
      '<p class="meta">Too many frames for a pairwise matrix; see the per-frame table.</p>';
  }}

  // Per-residue RMSF — one trace per chain per ensemble. Plotting every chain as one line
  // against resSeq draws a stroke back to the left edge wherever a new chain restarts its
  // numbering, which reads as a flexible region and is an artefact of the x axis.
  const DASHES = ['solid', 'dash', 'dot', 'dashdot'];
  const rmsfTraces = [];
  DATA.ensembles.forEach((e, i) => {{
    if (!e.rmsf || !e.rmsf.length) return;
    const chains = [];
    e.rmsf.forEach(r => {{ if (!chains.includes(r.chain)) chains.push(r.chain); }});
    chains.forEach((chain, ci) => {{
      const residues = e.rmsf.filter(r => r.chain === chain);
      rmsfTraces.push({{
        // Position along the selection, so consecutive residues are consecutive points
        // whatever the author numbering does. The real residue is in the hover text.
        x: residues.map((r, k) => k + 1),
        y: residues.map(r => r.rmsf_nm),
        name: chains.length > 1 ? e.label + ' chain ' + chain : e.label,
        mode: 'lines', type: 'scatter',
        line: {{ color: colors[i % 2], width: 2, dash: DASHES[ci % DASHES.length] }},
        hovertemplate: '%{{text}}<br>RMSF %{{y:.3f}} nm<extra></extra>',
        text: residues.map(r => e.label + ' chain ' + r.chain + ' ' + r.resName + r.resSeq)
      }});
    }});
  }});
  if (rmsfTraces.length) {{
    const l = layout(c, '{rmsf_axis}', 'RMSF (nm)');
    l.hovermode = 'closest';
    l.showlegend = rmsfTraces.length > 1;
    l.legend = {{ orientation: 'h', y: -0.22, font: {{ color: c.secondary }} }};
    if (l.showlegend) {{ l.margin.b = 86; }}
    Plotly.newPlot('rmsf-chart', rmsfTraces, l, CONFIG);
  }} else {{
    document.getElementById('rmsf-chart').innerHTML =
      '<p class="meta">RMSF unavailable for this system.</p>';
  }}
}}

draw();
window.matchMedia('(prefers-color-scheme: dark)').addEventListener('change', draw);

// 3D viewer: representative frames as animation frames.
if (DATA.pdb) {{
  const viewer = $3Dmol.createViewer('viewer', {{ backgroundColor: 'rgba(0,0,0,0)' }});
  viewer.addModelsAsFrames(DATA.pdb, 'pdb');
  viewer.setStyle({{}}, {{ cartoon: {{ color: 'spectrum' }} }});
  viewer.zoomTo();
  viewer.render();
  let playing = false;
  const button = document.getElementById('animate');
  if (button) {{
    button.addEventListener('click', () => {{
      playing = !playing;
      if (playing) {{
        viewer.animate({{ loop: 'forward', reps: 0, interval: 700 }});
        button.textContent = 'Pause';
      }} else {{
        viewer.stopAnimate();
        button.textContent = 'Animate frames';
      }}
    }});
  }}
}}
</script>
</body>
</html>
"""

VIEWER_BLOCK = """<div class="d-flex justify-content-between align-items-center mb-2">
    <div class="meta">{n} representative frame(s): {provenance}. Each model is a separate
      sampled conformation, not a step along a path. All are superposed on {aligned_on}, so
      what moves on screen is conformational change and not rigid-body tumbling.</div>
    <button id="animate" class="btn btn-sm btn-outline-secondary">Animate frames</button>
  </div>
  <div id="viewer"></div>"""


def status_pill(status):
    if status == "passed":
        return '<span class="pill"><span class="dot" style="background:var(--status-good)">' \
               "</span>&#10003; passed</span>"
    if status == "rejected":
        return '<span class="pill"><span class="dot" ' \
               'style="background:var(--status-critical)"></span>&#10007; rejected</span>'
    if status == "unavailable":
        return ('<span class="pill"><span class="dot" '
                'style="background:var(--status-warning)"></span>! not available</span>')
    if status in ("not_checked", ""):
        return '<span class="pill" style="color:var(--text-muted)">not checked</span>'
    return f'<span class="pill">{status}</span>'


def format_counts(counts):
    """Render a {name: count} mapping, or None so the row shows an em dash."""
    if not isinstance(counts, dict) or not counts:
        return None
    return ", ".join(
        f"{name} &times;{count}"
        for name, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    )


def alignment_sentence(sampling):
    """One sentence about the alignment the model was actually conditioned on.

    BioEmu2 is conditioned on an MSA, so its depth bounds what the ensemble can be. The
    report states the depth beside the results for that reason, and states separately whether
    the alignment reached the model: those two were able to disagree, and when they did, the
    run looked normal.
    """
    depth = sampling.get("msa_min_depth")
    source = sampling.get("msa_source")
    loaded = sampling.get("msa_loaded")

    if loaded is False and sampling.get("msa_how"):
        return (
            f"No alignment was used: {sampling['msa_how']}. BioEmu2 is conditioned on an "
            "alignment, so read everything below as the unconditioned case."
        )
    if depth is None:
        return (
            "The alignment depth was not recorded for this run, so a narrow ensemble cannot "
            "be attributed to a shallow alignment or ruled out as one."
        )

    sentence = (
        f"This run was conditioned on an alignment of {depth} row(s)"
        f"{f' from {source}' if source else ''}."
    )
    if depth < 32:
        sentence += (
            " That is shallow. BioEmu2 is conditioned on the alignment, and upstream names "
            "shallow or poor alignments — orphan and some viral proteins — as a cause of "
            "degraded ensembles, so treat a narrow or odd result here as unexplained rather "
            "than as a finding."
        )
    else:
        sentence += (
            " Depth is shown here because the ensemble is conditioned on it: a narrow result "
            "from a shallow alignment is not evidence of a rigid protein."
        )
    return sentence


def tile(label, value, unit=""):
    unit_html = f' <span class="unit">{unit}</span>' if unit else ""
    return (
        f'<div class="tile"><div class="label">{label}</div>'
        f'<div class="value">{value}{unit_html}</div></div>'
    )


def main():
    os.makedirs("./outputs", exist_ok=True)

    # No top_n here: the number of viewer frames is node 05's decision and is read from
    # ensemble_summary.json. A second copy of the parameter in this node could disagree with
    # the file it is describing, and the one in the file is the one that is true.
    summary = read_json("./inputs/ensemble_summary.json", "ensemble summary")
    sampling = read_json("./inputs/sampling_manifest.json", "sampling manifest", default={})
    propagation = read_json(
        "./inputs/propagation_manifest.json", "propagation manifest", default={}
    )
    rows = read_csv("./inputs/ensemble_metrics.csv")
    pdb_text = read_text("./inputs/ensemble_top.pdb", "multi-model PDB")

    ensembles_summary = summary.get("ensembles", {})
    sampled = ensembles_summary.get("sampled", {})
    if not sampled:
        fail("ensemble_summary.json has no 'sampled' ensemble; node 05 did not run correctly.")

    propagation_enabled = bool(summary.get("propagation_enabled"))

    # Per-ensemble chart series, built from the per-frame table so the charts and the table
    # can never disagree.
    ensembles = []
    for label in ("sampled", "propagated"):
        if label == "propagated" and not propagation_enabled:
            continue
        frame_rows = [r for r in rows if r.get("ensemble") == label]
        if not frame_rows:
            continue
        ensembles.append(
            {
                "label": label,
                "frames": [int(r["frame"]) for r in frame_rows],
                "rg": [number(r["radius_of_gyration_nm"], 4) for r in frame_rows],
                "rmsf": (ensembles_summary.get(label, {}) or {}).get("residue_rmsf", []),
            }
        )
    if not ensembles:
        fail("ensemble_metrics.csv contains no rows for the sampled ensemble.")

    settings = sampling.get("settings", {}) or {}
    physicality = sampled.get("physicality", {}) or {}
    rg_stats = sampled.get("radius_of_gyration_nm", {}) or {}
    pairwise = sampled.get("pairwise_rmsd") or {}

    pass_rate = physicality.get("pass_rate")
    pass_display = "not checked" if pass_rate is None else f"{round(pass_rate * 100)}%"

    propagated = ensembles_summary.get("propagated", {}) or {}
    prop_physicality = propagated.get("physicality", {}) or {}
    prop_pass_rate = prop_physicality.get("pass_rate")

    tiles = "".join(
        [
            tile("Frames", sampled.get("n_frames", "—")),
            tile("Mode", settings.get("mode", "—")),
            tile("Mean R<sub>g</sub>", fmt(rg_stats.get("mean"), 3), "nm"),
            tile("Max pairwise RMSD", fmt(pairwise.get("max_nm"), 3), "nm"),
            tile(
                "Clusters",
                sampled.get("n_clusters") if sampled.get("n_clusters") is not None else "—",
            ),
            tile("Physicality passed", pass_display),
        ]
    )

    legend_items = []
    palette_roles = ["--series-1", "--series-2"]
    for index, ensemble in enumerate(ensembles):
        role = palette_roles[index % 2]
        legend_items.append(
            f'<span><span class="dot" style="background:var({role})"></span>'
            f"{ensemble['label']}</span>"
        )
    legend = "".join(legend_items)

    table_rows = []
    for row in rows:
        table_rows.append(
            "<tr>"
            f"<td>{row.get('frame', '')}</td>"
            f"<td>{row.get('ensemble', '')}</td>"
            f"<td>{fmt(row.get('radius_of_gyration_nm'))}</td>"
            f"<td>{fmt(row.get('rmsd_to_medoid_nm'))}</td>"
            f"<td>{status_pill(row.get('physicality_status', ''))}</td>"
            f"<td class=\"meta\">{row.get('physicality_reason') or ''}</td>"
            "</tr>"
        )

    provenance_pairs = [
        ("System source", (summary.get("system") or {}).get("source")),
        (
            "Chains",
            ", ".join((summary.get("system") or {}).get("chain_ids") or []) or None,
        ),
        # Labelled for what it is. Upstream's own documentation for System.from_cif says
        # "polymer sequences include unresolved residues from the mmCIF entity tables", so this
        # is the construct BioEmu2 modelled, not the residue count of the input structure --
        # 1CLL is a 148-residue protein whose crystal resolves 144, and the ensemble has 148.
        ("Polymer residues (entity sequence)",
         (summary.get("system") or {}).get("n_polymer_residues")),
        ("Tokens", (summary.get("system") or {}).get("n_tokens")),
        # What else was in the system. Every one of these is sampled as ligand tokens, and the
        # difference between a cofactor and a cryoprotectant is not visible in a yes/no flag.
        (
            "Non-polymer components sampled",
            ", ".join(
                f"{item.get('ccd_code') or '?'} ({item.get('asym_id')})"
                for item in ((summary.get("system") or {}).get("ligand_components") or [])
            )
            or ("none" if (summary.get("system") or {}).get("has_ligands") is False else None),
        ),
        ("BioEmu2 version", sampling.get("bioemu2_version")),
        ("Checkpoint", settings.get("checkpoint")),
        ("Checkpoint revision", settings.get("hf_ckpt_revision")),
        ("Mode", settings.get("mode")),
        ("Sampler", settings.get("sampler")),
        (
            "Denoising steps",
            settings.get("num_steps") or "200 (the preset default)",
        ),
        ("Seed", settings.get("seed")),
        ("Batch size", settings.get("batch_size")),
        ("Precision", settings.get("precision")),
        ("Matmul precision", settings.get("matmul_precision")),
        ("Physicality filter", settings.get("physicality_filter")),
        ("Physicality not checked because", sampling.get("physicality_not_checked_because")),
        ("Rejection reasons", format_counts(sampling.get("rejection_reasons"))),
        ("MSA source", sampling.get("msa_source")),
        ("MSA minimum depth", sampling.get("msa_min_depth")),
        # What the model was given, which is not the same claim as what node 02 fetched.
        ("Alignment given to the model", sampling.get("msa_how")),
        ("Alignment files", ", ".join(sampling.get("msa_files") or []) or None),
        ("RMSD atom selection", sampled.get("rmsd_atom_selection")),
        (
            "Chains with no C-&alpha; atom",
            ", ".join(sampled.get("excluded_chains") or []) or None,
        ),
        ("Cluster threshold (nm)", sampled.get("rmsd_cluster_threshold_nm")),
        ("Medoid frame", sampled.get("medoid_frame")),
        ("Pairwise RMSD pairs measured", (sampled.get("pairwise_rmsd") or {}).get("n_pairs")),
        ("Propagation", "enabled" if propagation_enabled else "disabled"),
        ("Lag time (ns)", propagation.get("lag_time_ns") if propagation_enabled else None),
        (
            "Propagation sampler",
            propagation.get("sampler") if propagation_enabled else None,
        ),
        (
            "Propagated endpoints passing physicality",
            (
                None
                if not propagation_enabled
                else (
                    "not checked"
                    if prop_pass_rate is None
                    else f"{round(prop_pass_rate * 100)}%"
                )
            ),
        ),
        (
            "Alignment reused for propagation",
            propagation.get("msa_how") if propagation_enabled else None,
        ),
        ("GPU", sampling.get("gpu")),
        ("Sampling time (s)", sampling.get("elapsed_seconds")),
        (
            "Propagation time (s)",
            propagation.get("elapsed_seconds") if propagation_enabled else None,
        ),
    ]
    provenance = "".join(
        f'<tr><th style="width:38%">{label}</th>'
        f"<td><code>{'—' if value is None or value == '' else value}</code></td></tr>"
        for label, value in provenance_pairs
    )

    n_viewer = len(summary.get("viewer_frames") or [])
    n_from_clusters = summary.get("viewer_frames_from_clusters")
    if n_from_clusters is None:
        viewer_provenance = "selection not recorded by node 05"
    elif n_from_clusters == n_viewer:
        viewer_provenance = "all of them cluster centres"
    elif n_from_clusters == 0:
        viewer_provenance = (
            "none of them cluster centres, because no clustering was done — these are the "
            "first frames of the ensemble"
        )
    else:
        viewer_provenance = (
            f"{n_from_clusters} cluster centre(s) and {n_viewer - n_from_clusters} further "
            "frame(s) taken in order to fill the set"
        )
    aligned_frame = summary.get("viewer_frames_superposed_on_frame")
    aligned_on = (
        f"frame {aligned_frame}" if aligned_frame is not None else "the medoid frame"
    )
    viewer_block = (
        VIEWER_BLOCK.format(n=n_viewer, provenance=viewer_provenance, aligned_on=aligned_on)
        if pdb_text
        else '<p class="meta">No structures available for the 3D viewer.</p>'
    )

    mode = settings.get("mode", "equilibrium")
    headline = (
        f"{sampled.get('n_frames', '?')} conformation(s) sampled in {mode} mode"
        f" &middot; {sampled.get('n_atoms', '?')} atoms"
    )
    if propagation_enabled:
        headline += f" &middot; propagated {propagation.get('lag_time_ns')} ns"

    propagation_note = (
        " The propagated series holds one time-lagged endpoint per starting frame — lagged "
        "endpoints, not a continuous path."
        if propagation_enabled
        else " Propagation was disabled for this run."
    )

    msa_sentence = alignment_sentence(sampling)
    data = {
        "ensembles": ensembles,
        "rmsdMatrix": summary.get("rmsd_matrix") or [],
        "pdb": pdb_text,
    }

    # One RMSF point per selected atom, which is one per residue only when the selection is
    # C-alphas. Node 05 falls back to every atom for a system with no protein chain — 1BNA
    # gives 492 points over 24 nucleotides — and calling that "per-residue" would be a
    # four-fold overstatement of the resolution, with no way for a reader to tell.
    selection = (sampled.get("rmsd_atom_selection") or "").lower()
    per_atom = "ca" not in selection
    rmsf_caption = (
        ("Per-atom" if per_atom else "Per-residue")
        + " RMSF about the mean structure, after superposition on the medoid frame, computed "
        + f"over {sampled.get('rmsd_atom_selection') or 'the selected atoms'}. Peaks are the "
        "flexible regions — loops and termini."
    )
    rmsf_axis = ("Atom" if per_atom else "Residue") + " position along the selection"

    html = HTML.format(
        bootstrap=BOOTSTRAP_CDN,
        plotly=PLOTLY_CDN,
        threedmol=THREEDMOL_CDN,
        headline=headline,
        propagation_note=propagation_note,
        msa_sentence=msa_sentence,
        tiles=tiles,
        legend=legend,
        rmsf_caption=rmsf_caption,
        rmsf_axis=rmsf_axis,
        rows="".join(table_rows),
        provenance=provenance,
        viewer_block=viewer_block,
        data_json=json.dumps(data),
    )

    with open("./outputs/report.html", "w") as handle:
        handle.write(html)

    size_kb = round(os.path.getsize("./outputs/report.html") / 1024, 1)
    print(f"Wrote ./outputs/report.html ({size_kb} KB)", flush=True)
    print(
        f"  {len(rows)} frame row(s), {len(ensembles)} ensemble(s), "
        f"{n_viewer} viewer frame(s) ({n_from_clusters} cluster centre(s))",
        flush=True,
    )


if __name__ == "__main__":
    main()
