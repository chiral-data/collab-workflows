"""Generate self-contained HTML report: Tg, thermal expansion, dimensional
stability, and a Mol* view of the representative cell structure."""
import datetime
import json
import pathlib

tg      = json.loads(pathlib.Path("inputs/tg_report.json").read_text())
build   = json.loads(pathlib.Path("inputs/build_report.json").read_text())
cell_pdb_text = pathlib.Path("inputs/cell.pdb").read_text()


def fmt(val, digits=1, unit=""):
    if val is None:
        return "N/A"
    return f"{val:.{digits}f}{(' ' + unit) if unit else ''}"


resin        = tg["resin_type"]
resin_name   = build.get("resin_name", resin)
crystallinity = tg["crystallinity"]
tg_c         = tg.get("tg_c")
tg_reliable  = tg.get("tg_reliable", True)
tg_lit       = tg.get("literature_tg_c")
cte_glassy   = tg.get("cte_glassy_per_c")
cte_rubbery  = tg.get("cte_rubbery_per_c")
vol_change   = tg.get("volume_change_pct")
selected_t   = tg.get("selected_temperature_c")
matched_t    = tg.get("matched_series_temp_c")
series       = sorted(tg.get("density_temp_series", []), key=lambda s: s["temp_c"])
n_chains     = build.get("n_chains", "?")
melt_temp    = build.get("melt_temp_c", "?")
pack_frac    = build.get("pack_density_frac", "?")

# ── Density/specific-volume vs temperature chart (points + fitted lines) ─────
SVG_W, SVG_H, PAD_L, PAD_B, PAD_T, PAD_R = 720, 260, 56, 34, 14, 14

temps = [s["temp_c"] for s in series]
vols  = [s["avg_specific_volume_cm3_g"] for s in series]
t_min, t_max = min(temps), max(temps)
v_min, v_max = min(vols), max(vols)
t_rng = max(t_max - t_min, 1.0)
v_rng = max(v_max - v_min, 1e-6)


def sx(t):
    return PAD_L + (t - t_min) / t_rng * (SVG_W - PAD_L - PAD_R)


def sy(v):
    return SVG_H - PAD_B - (v - v_min) / v_rng * (SVG_H - PAD_T - PAD_B)


n_glassy  = tg.get("fit_split_n_glassy")
n_rubbery = tg.get("fit_split_n_rubbery")

points_svg = "\n".join(
    f'  <circle cx="{sx(t):.1f}" cy="{sy(v):.1f}" r="4.5" fill="#1a5cff" stroke="#fff" stroke-width="1.5"/>'
    for t, v in zip(temps, vols)
)

fit_lines_svg = ""
if n_glassy and n_rubbery and cte_glassy is not None and cte_rubbery is not None:
    glassy_pts  = list(zip(temps[:n_glassy],  vols[:n_glassy]))
    rubbery_pts = list(zip(temps[n_glassy:],  vols[n_glassy:]))

    def line_through(pts, x0, x1, slope):
        x_anchor, y_anchor = pts[0]
        y0 = y_anchor + slope * (x0 - x_anchor)
        y1 = y_anchor + slope * (x1 - x_anchor)
        return y0, y1

    gx0, gx1 = glassy_pts[0][0], glassy_pts[-1][0]
    gy0, gy1 = glassy_pts[0][1], glassy_pts[-1][1]
    rx0, rx1 = rubbery_pts[0][0], rubbery_pts[-1][0]
    ry0, ry1 = rubbery_pts[0][1], rubbery_pts[-1][1]

    # extend both fitted lines to the estimated Tg so they visibly meet at the
    # kink — only when the fit is reliable (a degenerate near-parallel fit
    # would extrapolate the "intersection" far off the chart)
    if tg_c is not None and tg_reliable:
        gy_at_tg = gy0 + (gy1 - gy0) / max(gx1 - gx0, 1e-9) * (tg_c - gx0)
        ry_at_tg = ry0 + (ry1 - ry0) / max(rx1 - rx0, 1e-9) * (tg_c - rx0)
        fit_lines_svg += (
            f'  <line x1="{sx(gx1):.1f}" y1="{sy(gy1):.1f}" x2="{sx(tg_c):.1f}" y2="{sy(gy_at_tg):.1f}" '
            f'stroke="#0f766e" stroke-width="2" stroke-dasharray="4 3"/>\n'
            f'  <line x1="{sx(rx0):.1f}" y1="{sy(ry0):.1f}" x2="{sx(tg_c):.1f}" y2="{sy(ry_at_tg):.1f}" '
            f'stroke="#d97706" stroke-width="2" stroke-dasharray="4 3"/>\n'
            f'  <line x1="{sx(tg_c):.1f}" y1="{PAD_T}" x2="{sx(tg_c):.1f}" y2="{SVG_H-PAD_B}" '
            f'stroke="#0c1829" stroke-width="1" stroke-dasharray="2 3" opacity="0.5"/>\n'
        )
    fit_lines_svg += (
        f'  <line x1="{sx(gx0):.1f}" y1="{sy(gy0):.1f}" x2="{sx(gx1):.1f}" y2="{sy(gy1):.1f}" stroke="#0f766e" stroke-width="2.5"/>\n'
        f'  <line x1="{sx(rx0):.1f}" y1="{sy(ry0):.1f}" x2="{sx(rx1):.1f}" y2="{sy(ry1):.1f}" stroke="#d97706" stroke-width="2.5"/>\n'
    )

axis_svg = (
    f'  <line x1="{PAD_L}" y1="{SVG_H-PAD_B}" x2="{SVG_W-PAD_R}" y2="{SVG_H-PAD_B}" stroke="#d8e0ef" stroke-width="1"/>\n'
    f'  <line x1="{PAD_L}" y1="{PAD_T}" x2="{PAD_L}" y2="{SVG_H-PAD_B}" stroke="#d8e0ef" stroke-width="1"/>\n'
)
tick_labels_svg = "\n".join(
    f'  <text x="{sx(t):.1f}" y="{SVG_H-PAD_B+16}" font-size="10" fill="#6b7fa3" text-anchor="middle" font-family="monospace">{t:.0f}</text>'
    for t in temps
)
tg_label_svg = (
    f'  <text x="{sx(tg_c):.1f}" y="{PAD_T-2}" font-size="10.5" fill="#0c1829" text-anchor="middle" '
    f'font-family="monospace" font-weight="bold">Tg={tg_c:.0f}&#8451;</text>'
    if (tg_c is not None and tg_reliable) else
    ('  <text x="' + str(SVG_W/2) + '" y="' + str(PAD_T + 10) +
     '" font-size="11" fill="#d97706" text-anchor="middle" font-family="monospace" '
     'font-weight="bold">Tg fit unreliable — segments nearly parallel</text>'
     if tg_c is not None else "")
)

chart_svg = f"""<svg width="100%" height="260" viewBox="0 0 {SVG_W} {SVG_H}" xmlns="http://www.w3.org/2000/svg">
{axis_svg}{fit_lines_svg}{points_svg}
{tick_labels_svg}
{tg_label_svg}
</svg>"""

# ── Cross-resin Tg bar chart vs. literature ──────────────────────────────────
LIT_TG = [
    ("PP (ref.)", -10),
    ("PBT",        45),
    ("PA66",       50),
    ("PPS",        88),
    ("PEEK",      143),
]
# If the MD fit is unreliable, don't let a wild extrapolated value (e.g.
# -1000+ C) corrupt the chart's scale — fall back to plotting this resin's
# literature Tg instead, clearly labeled as literature-only.
if tg_c is not None and tg_reliable:
    this_label = f"{resin} (this run, est.)"
    this_val   = tg_c
else:
    this_label = f"{resin} (literature — MD fit unreliable)"
    this_val   = tg_lit if tg_lit is not None else 0.0

all_bars = [(l, v) for l, v in LIT_TG if l.split()[0] != resin]
BAR_MIN  = min(-30.0, this_val - 10)
BAR_MAX  = max(160.0, this_val + 10)
BAR_RNG  = max(BAR_MAX - BAR_MIN, 1.0)


def bar_row(label, val, highlight=False):
    pct   = max(0.0, min(100.0, (val - BAR_MIN) / BAR_RNG * 100))
    color = "#1a5cff" if highlight else "#b0bdd8"
    n_cls = "bar-name highlight" if highlight else "bar-name"
    v_cls = "bar-val highlight"  if highlight else "bar-val"
    return f"""      <div class="bar-row">
        <div class="{n_cls}">{label}</div>
        <div class="bar-track"><div class="bar-fill" style="width:{pct:.2f}%;background:{color}"></div></div>
        <div class="{v_cls}">{val:.0f} &deg;C</div>
      </div>"""


bar_rows = [bar_row(l, v) for l, v in all_bars]
insert_at = sum(1 for _, v in all_bars if v < this_val)
bar_rows.insert(insert_at, bar_row(this_label, this_val, True))
bars_html = "\n".join(bar_rows)

# ── Dimensional-stability table ──────────────────────────────────────────────
series_rows = "\n".join(
    f"<tr><td>{s['stage']}</td><td class='val'>{s['temp_c']:.0f} &deg;C</td>"
    f"<td class='val'>{s['avg_density_kg_m3']:.1f}</td>"
    f"<td class='val'>{s['avg_specific_volume_cm3_g']:.4f}</td></tr>"
    for s in series
)

# ── Mol* viewer ───────────────────────────────────────────────────────────────
pdb_js_escaped = json.dumps(cell_pdb_text)

molstar_script = f"""
(function() {{
  var pdbData = {pdb_js_escaped};
  molstar.Viewer.create('molstar-viewer', {{
    layoutIsExpanded: false,
    layoutShowControls: false,
    layoutShowSequence: false,
    layoutShowLog: false,
    layoutShowRemoteState: false,
    viewportShowExpand: true,
  }}).then(function(viewer) {{
    var p = viewer.plugin;
    p.builders.data.rawData({{ data: pdbData, label: '{resin} cell' }})
      .then(function(d) {{ return p.builders.structure.parseTrajectory(d, 'pdb'); }})
      .then(function(t) {{ return p.builders.structure.hierarchy.applyPreset(t, 'default'); }})
      .catch(function(e) {{
        document.getElementById('molstar-viewer').innerHTML =
          '<p style="color:#d97706;padding:16px">Viewer error: ' + e + '</p>';
      }});
  }}).catch(function(e) {{
    document.getElementById('molstar-viewer').innerHTML =
      '<p style="color:#d97706;padding:16px">Viewer error: ' + e + '</p>';
  }});
}})();
"""

# ── HTML ──────────────────────────────────────────────────────────────────────
html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Heat-Resistant Plastics Report &mdash; {resin} @ {crystallinity} crystallinity</title>
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/molstar@latest/build/viewer/molstar.css">
<script src="https://cdn.jsdelivr.net/npm/molstar@latest/build/viewer/molstar.js"></script>
<style>
  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
  :root {{
    --ink: #0c1829; --ground: #f5f7fb; --surface: #fff; --accent: #1a5cff;
    --accent-dim: #eef2ff; --mid: #6b7fa3; --border: #d8e0ef; --warn: #d97706;
    --ff-d: Georgia, 'Times New Roman', serif;
    --ff-b: system-ui, -apple-system, BlinkMacSystemFont, sans-serif;
    --ff-m: ui-monospace, 'Cascadia Code', 'Fira Code', 'Courier New', monospace;
  }}
  body {{ font-family: var(--ff-b); font-size: 15px; line-height: 1.65; color: var(--ink); background: var(--ground); padding: 40px 20px 80px; }}
  .page {{ max-width: 800px; margin: 0 auto; }}
  .header {{ background: var(--ink); border-radius: 10px 10px 0 0; padding: 32px 36px 28px; position: relative; overflow: hidden; }}
  .header::before {{ content: ''; position: absolute; inset: 0; background-image: radial-gradient(circle, rgba(26,92,255,.18) 1px, transparent 1px); background-size: 28px 28px; }}
  .header-eyebrow {{ font-family: var(--ff-m); font-size: 11px; letter-spacing: .12em; text-transform: uppercase; color: var(--accent); margin-bottom: 10px; }}
  .header h1 {{ font-family: var(--ff-d); font-size: 1.75rem; font-weight: normal; color: #fff; line-height: 1.25; margin-bottom: 16px; }}
  .badges {{ display: flex; gap: 8px; flex-wrap: wrap; }}
  .badge {{ font-family: var(--ff-m); font-size: 11px; letter-spacing: .06em; padding: 3px 10px; border-radius: 3px; background: rgba(26,92,255,.25); color: #a8c0ff; border: 1px solid rgba(26,92,255,.4); }}
  .badge.primary {{ background: var(--accent); color: #fff; border-color: var(--accent); }}
  .header-date {{ font-family: var(--ff-m); font-size: 11px; color: var(--mid); margin-top: 14px; }}
  .kpi-strip {{ display: grid; grid-template-columns: repeat(4, 1fr); background: var(--surface); border: 1px solid var(--border); border-top: none; }}
  .kpi {{ padding: 20px 20px; border-right: 1px solid var(--border); }}
  .kpi:last-child {{ border-right: none; }}
  .kpi-label {{ font-family: var(--ff-m); font-size: 10px; letter-spacing: .08em; text-transform: uppercase; color: var(--mid); margin-bottom: 6px; }}
  .kpi-value {{ font-family: var(--ff-d); font-size: 1.7rem; font-weight: bold; color: var(--ink); font-variant-numeric: tabular-nums; line-height: 1; }}
  .kpi-unit {{ font-family: var(--ff-m); font-size: 10.5px; color: var(--mid); margin-top: 4px; }}
  .section {{ background: var(--surface); border: 1px solid var(--border); border-top: none; padding: 28px 36px; }}
  .section:last-of-type {{ border-radius: 0 0 10px 10px; }}
  .section-label {{ font-family: var(--ff-m); font-size: 10.5px; letter-spacing: .12em; text-transform: uppercase; color: var(--accent); display: flex; align-items: center; gap: 10px; margin-bottom: 16px; }}
  .section-label::after {{ content: ''; flex: 1; height: 1px; background: var(--accent-dim); }}
  .data-table {{ width: 100%; border-collapse: collapse; }}
  .data-table th, .data-table td {{ text-align: left; padding: 8px 12px; border-bottom: 1px solid var(--border); font-size: 13.5px; }}
  .data-table th {{ font-family: var(--ff-m); font-size: 11px; letter-spacing: .05em; color: var(--mid); background: var(--ground); font-weight: 500; }}
  .data-table .val {{ font-family: var(--ff-m); font-size: 13px; font-variant-numeric: tabular-nums; }}
  .data-table tr:last-child td {{ border-bottom: none; }}
  .chart-wrap {{ background: var(--ground); border: 1px solid var(--border); border-radius: 6px; padding: 16px; overflow-x: auto; }}
  .legend {{ display: flex; gap: 18px; font-family: var(--ff-m); font-size: 11px; color: var(--mid); margin-top: 8px; }}
  .legend-swatch {{ display: inline-block; width: 10px; height: 10px; border-radius: 2px; margin-right: 5px; vertical-align: middle; }}
  .bar-chart {{ display: flex; flex-direction: column; gap: 10px; }}
  .bar-row {{ display: grid; grid-template-columns: 170px 1fr 70px; align-items: center; gap: 10px; }}
  .bar-name {{ font-size: 13px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }}
  .bar-name.highlight {{ font-weight: 600; color: var(--accent); }}
  .bar-track {{ background: var(--ground); border-radius: 3px; height: 18px; overflow: hidden; }}
  .bar-fill {{ height: 100%; border-radius: 3px; min-width: 2px; }}
  .bar-val {{ font-family: var(--ff-m); font-size: 12px; color: var(--mid); text-align: right; font-variant-numeric: tabular-nums; }}
  .bar-val.highlight {{ color: var(--accent); font-weight: 600; }}
  .callout-note {{ background: #fffbeb; border-left: 3px solid var(--warn); padding: 12px 16px; border-radius: 0 4px 4px 0; font-size: 13px; color: #78350f; margin-top: 16px; }}
  .callout-note + .callout-note {{ margin-top: 8px; }}
  .callout-note code {{ font-family: var(--ff-m); font-size: 12px; background: rgba(0,0,0,.06); padding: 1px 5px; border-radius: 3px; }}
  #molstar-viewer {{ position: relative; width: 100%; height: 420px; border-radius: 6px; overflow: hidden; border: 1px solid var(--border); }}
  .footer {{ margin-top: 28px; font-family: var(--ff-m); font-size: 11px; color: var(--mid); text-align: center; letter-spacing: .06em; }}
  @media (max-width: 600px) {{
    .kpi-strip {{ grid-template-columns: 1fr 1fr; }}
    .kpi {{ border-right: none; border-bottom: 1px solid var(--border); }}
    .bar-row {{ grid-template-columns: 110px 1fr 55px; }}
    .header, .section {{ padding: 20px; }}
  }}
</style>
</head>
<body>
<div class="page">

  <div class="header">
    <div class="header-eyebrow">workflow-032 &middot; AMECC Theme 3 &middot; Heat-Resistant Plastics</div>
    <h1>{resin_name} &mdash; {crystallinity} crystallinity<br>Protecting EV Battery Cases &amp; Separators from Heat</h1>
    <div class="badges">
      <span class="badge primary">{resin}</span>
      <span class="badge">{crystallinity}</span>
      <span class="badge">{selected_t:.0f}&deg;C selected</span>
      <span class="badge">GAFF2</span>
      <span class="badge">GROMACS 2023.2</span>
    </div>
    <div class="header-date">
      Generated {datetime.date.today()} &middot; {n_chains} chains &middot; melt-quench ladder (melt&rarr;200&rarr;150&rarr;80&rarr;25&thinsp;&deg;C)
    </div>
  </div>

  <div class="kpi-strip">
    <div class="kpi">
      <div class="kpi-label">Tg (estimated)</div>
      <div class="kpi-value" style="{'' if tg_reliable else 'color:var(--warn)'}">{fmt(tg_c, 0) if tg_reliable else 'N/R'}</div>
      <div class="kpi-unit">&deg;C &middot; lit. {fmt(tg_lit, 0)}&thinsp;&deg;C{'' if tg_reliable else ' &middot; unreliable fit'}</div>
    </div>
    <div class="kpi">
      <div class="kpi-label">CTE Glassy</div>
      <div class="kpi-value">{fmt((cte_glassy or 0)*1e4, 2)}</div>
      <div class="kpi-unit">&times;10&#8315;&#8308; /&deg;C</div>
    </div>
    <div class="kpi">
      <div class="kpi-label">CTE Rubbery</div>
      <div class="kpi-value">{fmt((cte_rubbery or 0)*1e4, 2)}</div>
      <div class="kpi-unit">&times;10&#8315;&#8308; /&deg;C</div>
    </div>
    <div class="kpi">
      <div class="kpi-label">&Delta;Volume @ {matched_t:.0f}&deg;C</div>
      <div class="kpi-value">{fmt(vol_change, 2)}</div>
      <div class="kpi-unit">% vs 25&deg;C</div>
    </div>
  </div>

  <div class="section">
    <div class="section-label">Simulation Setup</div>
    <table class="data-table">
      <tr><th>Parameter</th><th>Value</th></tr>
      <tr><td>Resin</td><td class="val">{resin_name} ({resin})</td></tr>
      <tr><td>Crystallinity</td><td class="val">{crystallinity} (packing density {float(pack_frac)*100:.0f}% of target)</td></tr>
      <tr><td>Selected temperature</td><td class="val">{selected_t:.0f} &deg;C (matched to {matched_t:.0f} &deg;C ladder point)</td></tr>
      <tr><td>Cell composition</td><td class="val">{n_chains} &times; oligomer chains</td></tr>
      <tr><td>Force field</td><td class="val">GAFF2 &mdash; antechamber AM1-BCC charges</td></tr>
      <tr><td>Melt-quench ladder</td><td class="val">melt ({melt_temp}&deg;C) &rarr; 200 &rarr; 150 &rarr; 80 &rarr; 25 &deg;C</td></tr>
    </table>
  </div>

  <div class="section">
    <div class="section-label">Glass Transition &mdash; Density/Specific-Volume Kink Fit</div>
    <div class="chart-wrap">{chart_svg}</div>
    <div class="legend">
      <span><span class="legend-swatch" style="background:#1a5cff"></span>MD series point</span>
      <span><span class="legend-swatch" style="background:#0f766e"></span>Glassy fit</span>
      <span><span class="legend-swatch" style="background:#d97706"></span>Rubbery fit</span>
    </div>
    <div class="callout-note">
      &#9888; Tg is the intersection of the two fitted lines over a 5-point ladder &mdash;
      the minimum viable case for a bilinear kink fit. Increase <code>PARAM_STAGE_TIME_PS</code>
      for a better-converged curve and a more reliable Tg.
    </div>
    {'<div class="callout-note">&#9888; <b>This fit is unreliable</b> &mdash; the glassy and rubbery '
     'segment slopes are nearly parallel (a common symptom of short, under-converged default runs), '
     'so their intersection falls far outside the sampled temperature range and is not a physically '
     'meaningful Tg. Increase <code>PARAM_STAGE_TIME_PS</code> to separate the two regimes.</div>'
     if not tg_reliable else ''}
  </div>

  <div class="section">
    <div class="section-label">Tg Comparison &mdash; Heat-Resistant Resin Ranking</div>
    <p style="font-size:13px;color:var(--mid);margin-bottom:18px;">
      Higher Tg &rarr; retains stiffness/shape at higher temperature &rarr; better battery-case/separator protection.
    </p>
    <div class="bar-chart">
{bars_html}
    </div>
  </div>

  <div class="section">
    <div class="section-label">Dimensional Stability &mdash; Density &amp; Specific Volume by Stage</div>
    <table class="data-table">
      <tr><th>Stage</th><th>Temperature</th><th>Density (kg/m&sup3;)</th><th>Specific volume (cm&sup3;/g)</th></tr>
      {series_rows}
    </table>
  </div>

  <div class="section">
    <div class="section-label">Representative Structure</div>
    <div id="molstar-viewer"></div>
  </div>

  <div class="section">
    <div class="callout-note">
      &#9888; Force field: aromatic backbones (PPS, PEEK, PBT) are chemically better suited to OPLS-AA,
      but this pipeline reuses the GAFF2/antechamber pipeline shared with workflow-030/031 for
      auto-parametrization consistency &mdash; see issue #196 for the tradeoff.
    </div>
    <div class="callout-note">
      &#9888; Crystallinity is approximated via initial packing density fraction, not a true
      semi-crystalline lattice &mdash; trend comparison only.
    </div>
  </div>

  <div class="footer">
    workflow-032 &middot; AMECC Theme 3 &middot; {datetime.date.today()}
  </div>

</div>
<script>{molstar_script}</script>
</body>
</html>
"""

pathlib.Path("outputs/report.html").write_text(html)

summary = {
    "resin_type":             resin,
    "crystallinity":          crystallinity,
    "selected_temperature_c": selected_t,
    "tg_c":                   tg_c,
    "tg_reliable":            tg_reliable,
    "literature_tg_c":        tg_lit,
    "cte_glassy_per_c":       cte_glassy,
    "cte_rubbery_per_c":      cte_rubbery,
    "volume_change_pct":      vol_change,
}
pathlib.Path("outputs/summary.json").write_text(json.dumps(summary, indent=2))

print(f"  Report  -> outputs/report.html")
print(f"  Summary -> outputs/summary.json")
print(f"  Tg: {fmt(tg_c, 0)} C (lit. {fmt(tg_lit, 0)} C)")
