"""Generate HTML report from corrected properties."""
import json, pathlib, datetime

props  = json.loads(pathlib.Path("inputs/corrected_properties.json").read_text())
build  = json.loads(pathlib.Path("inputs/build_report.json").read_text())

def fmt(val, digits=3, unit=""):
    if val is None:
        return "N/A"
    return f"{val:.{digits}f}{(' ' + unit) if unit else ''}"

def tr(label, *cols):
    cells = "".join(f"<td class='val'>{c}</td>" for c in cols)
    return f"<tr><td>{label}</td>{cells}</tr>"

resin       = props["resin_type"]
temp_c      = props["temperature_c"]
fiber_wt    = props["fiber_loading_wt_pct"]
gf_applied  = props["gf_correction_applied"]
E_md        = props.get("youngs_modulus_gpa")
E_corr      = props.get("corrected_youngs_modulus_gpa")
rho_md      = props.get("avg_density_kg_m3")
rho_corr    = props.get("corrected_density_kg_m3")
spec_mod    = props.get("corrected_specific_modulus_kNm_kg")
note        = props.get("note", "")
prod_time   = props.get("prod_time_ps", "?")
deform_time = props.get("deform_time_ps", "?")
n_chains    = build.get("n_chains", "?")
fiber_vf    = props.get("fiber_volume_fraction", 0.0)
melt_temp   = build.get("melt_temp_c", "?")

# ── Density convergence ───────────────────────────────────────────────────────
xvg_lines = [l for l in pathlib.Path("inputs/density_prod.xvg").read_text().splitlines()
             if l and not l.startswith(("#", "@"))]
density_series = [float(l.split()[1]) for l in xvg_lines]
n_third   = max(len(density_series) // 3, 1)
rho_final = sum(density_series[-n_third:]) / n_third
rho_init  = density_series[0] if density_series else 0.0

# SVG sparkline — downsample to ~120 pts, normalise to viewBox 720×80
step  = max(len(density_series) // 120, 1)
spark = density_series[::step]
s_min = min(spark)
s_max = max(spark)
s_rng = max(s_max - s_min, 1.0)
SVG_W, SVG_H, PAD = 720, 80, 6

def sx(i):
    return int(i / max(len(spark)-1, 1) * SVG_W)

def sy(v):
    return int(SVG_H - PAD - (v - s_min) / s_rng * (SVG_H - 2*PAD))

pts      = " ".join(f"{sx(i)},{sy(v)}" for i, v in enumerate(spark))
area_pts = f"0,{SVG_H} " + pts + f" {SVG_W},{SVG_H}"

sparkline = f"""
<svg width="100%" height="80" viewBox="0 0 {SVG_W} {SVG_H}"
     preserveAspectRatio="none" xmlns="http://www.w3.org/2000/svg">
  <line x1="0" y1="{int(SVG_H*0.25)}" x2="{SVG_W}" y2="{int(SVG_H*0.25)}"
        stroke="#d8e0ef" stroke-width="1"/>
  <line x1="0" y1="{int(SVG_H*0.5)}"  x2="{SVG_W}" y2="{int(SVG_H*0.5)}"
        stroke="#d8e0ef" stroke-width="1"/>
  <line x1="0" y1="{int(SVG_H*0.75)}" x2="{SVG_W}" y2="{int(SVG_H*0.75)}"
        stroke="#d8e0ef" stroke-width="1"/>
  <polygon points="{area_pts}" fill="rgba(26,92,255,0.08)"/>
  <polyline points="{pts}" fill="none" stroke="#1a5cff"
            stroke-width="2" stroke-linejoin="round" stroke-linecap="round"/>
  <circle cx="{sx(0)}" cy="{sy(spark[0])}" r="4" fill="#1a5cff"/>
  <circle cx="{sx(len(spark)-1)}" cy="{sy(spark[-1])}" r="4" fill="#1a5cff"/>
</svg>"""

# ── Bar chart ─────────────────────────────────────────────────────────────────
REF = [
    ("PP neat (lit.)",       1.4),
    ("PP+30wt%GF (lit.)",    4.5),
    ("PA6+30wt%GF (lit.)",   8.0),
    ("Steel",               26.0),
    ("Al 6061",             26.0),
    ("CF composite",       150.0),
]
this_abs   = abs(spec_mod) if spec_mod else 0.0
this_label = f"{resin}+{fiber_wt:.0f}wt%GF (this run)"
BAR_MAX    = max(200.0, this_abs * 1.05, 155.0)

def bar_row(label, val, highlight=False):
    pct   = max(0, min(100, abs(val) / BAR_MAX * 100))
    color = "#1a5cff" if highlight else "#b0bdd8"
    n_cls = "bar-name highlight" if highlight else "bar-name"
    v_cls = "bar-val highlight"  if highlight else "bar-val"
    v_str = fmt(spec_mod, 1) if highlight else f"{val:.1f}"
    return f"""      <div class="bar-row">
        <div class="{n_cls}">{label}</div>
        <div class="bar-track"><div class="bar-fill" style="width:{pct:.2f}%;background:{color}"></div></div>
        <div class="{v_cls}">{v_str}</div>
      </div>"""

bar_rows = [bar_row(l, v) for l, v in REF]
insert_at = sum(1 for _, v in REF if v < this_abs)
bar_rows.insert(insert_at, bar_row(this_label, this_abs, True))
bars_html = "\n".join(bar_rows)

# ── GF correction section ─────────────────────────────────────────────────────
gf_section = ""
if gf_applied:
    gf_section = f"""
  <div class="section">
    <div class="section-label">Glass-Fiber Correction (Halpin-Tsai)</div>
    <p style="font-size:13px;color:var(--mid);margin-bottom:12px;">
      Halpin-Tsai longitudinal + transverse, Christensen-Waals random-orientation average.
    </p>
    <table class="data-table">
      <tr><th>Parameter</th><th>Value</th></tr>
      <tr><td>Fiber type</td><td class="val">E-glass</td></tr>
      <tr><td>E<sub>fiber</sub></td><td class="val">{props.get('E_fiber_gpa', 72.0)} GPa</td></tr>
      <tr><td>Aspect ratio (l/d)</td><td class="val">{props.get('aspect_ratio', 20)}</td></tr>
      <tr><td>V<sub>f</sub></td><td class="val">{fiber_vf:.3f}</td></tr>
    </table>
    <div class="callout-note">&#9888; {note}</div>
  </div>"""

# ── HTML ──────────────────────────────────────────────────────────────────────
html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Polymer MD Report &mdash; {resin} @ {temp_c:.0f}&thinsp;&deg;C / {fiber_wt:.0f}&thinsp;wt% GF</title>
<style>
  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
  :root {{
    --ink:        #0c1829;
    --ground:     #f5f7fb;
    --surface:    #ffffff;
    --accent:     #1a5cff;
    --accent-dim: #eef2ff;
    --mid:        #6b7fa3;
    --border:     #d8e0ef;
    --warn:       #d97706;
    --ff-d: Georgia, 'Times New Roman', serif;
    --ff-b: system-ui, -apple-system, BlinkMacSystemFont, sans-serif;
    --ff-m: ui-monospace, 'Cascadia Code', 'Fira Code', 'Courier New', monospace;
  }}
  body {{
    font-family: var(--ff-b); font-size: 15px; line-height: 1.65;
    color: var(--ink); background: var(--ground); padding: 40px 20px 80px;
  }}
  .page {{ max-width: 800px; margin: 0 auto; }}
  .header {{
    background: var(--ink); border-radius: 10px 10px 0 0;
    padding: 32px 36px 28px; position: relative; overflow: hidden;
  }}
  .header::before {{
    content: ''; position: absolute; inset: 0;
    background-image: radial-gradient(circle, rgba(26,92,255,.18) 1px, transparent 1px);
    background-size: 28px 28px; pointer-events: none;
  }}
  .header-eyebrow {{
    font-family: var(--ff-m); font-size: 11px; letter-spacing: .12em;
    text-transform: uppercase; color: var(--accent); margin-bottom: 10px;
  }}
  .header h1 {{
    font-family: var(--ff-d); font-size: 1.75rem; font-weight: normal;
    color: #fff; line-height: 1.25; text-wrap: balance; margin-bottom: 16px;
  }}
  .badges {{ display: flex; gap: 8px; flex-wrap: wrap; }}
  .badge {{
    font-family: var(--ff-m); font-size: 11px; letter-spacing: .06em;
    padding: 3px 10px; border-radius: 3px;
    background: rgba(26,92,255,.25); color: #a8c0ff;
    border: 1px solid rgba(26,92,255,.4);
  }}
  .badge.primary {{ background: var(--accent); color: #fff; border-color: var(--accent); }}
  .header-date {{ font-family: var(--ff-m); font-size: 11px; color: var(--mid); margin-top: 14px; }}
  .kpi-strip {{
    display: grid; grid-template-columns: repeat(3, 1fr);
    background: var(--surface); border: 1px solid var(--border); border-top: none;
  }}
  .kpi {{ padding: 24px 28px; border-right: 1px solid var(--border); }}
  .kpi:last-child {{ border-right: none; }}
  .kpi-label {{
    font-family: var(--ff-m); font-size: 10.5px; letter-spacing: .1em;
    text-transform: uppercase; color: var(--mid); margin-bottom: 6px;
  }}
  .kpi-value {{
    font-family: var(--ff-d); font-size: 2.2rem; font-weight: bold;
    color: var(--ink); font-variant-numeric: tabular-nums; line-height: 1;
  }}
  .kpi-unit {{ font-family: var(--ff-m); font-size: 11px; color: var(--mid); margin-top: 4px; }}
  .section {{
    background: var(--surface); border: 1px solid var(--border);
    border-top: none; padding: 28px 36px;
  }}
  .section:last-of-type {{ border-radius: 0 0 10px 10px; }}
  .section-label {{
    font-family: var(--ff-m); font-size: 10.5px; letter-spacing: .12em;
    text-transform: uppercase; color: var(--accent);
    display: flex; align-items: center; gap: 10px; margin-bottom: 16px;
  }}
  .section-label::after {{ content: ''; flex: 1; height: 1px; background: var(--accent-dim); }}
  .data-table {{ width: 100%; border-collapse: collapse; }}
  .data-table th, .data-table td {{
    text-align: left; padding: 8px 12px;
    border-bottom: 1px solid var(--border); font-size: 13.5px;
  }}
  .data-table th {{
    font-family: var(--ff-m); font-size: 11px; letter-spacing: .05em;
    color: var(--mid); background: var(--ground); font-weight: 500;
  }}
  .data-table .val {{ font-family: var(--ff-m); font-size: 13px; font-variant-numeric: tabular-nums; }}
  .data-table tr:last-child td {{ border-bottom: none; }}
  .spark-wrap {{
    background: var(--ground); border: 1px solid var(--border);
    border-radius: 6px; padding: 16px; overflow-x: auto;
  }}
  .spark-meta {{
    display: flex; justify-content: space-between;
    font-family: var(--ff-m); font-size: 11px; color: var(--mid); margin-bottom: 8px;
  }}
  .bar-chart {{ display: flex; flex-direction: column; gap: 10px; }}
  .bar-row {{ display: grid; grid-template-columns: 190px 1fr 90px; align-items: center; gap: 10px; }}
  .bar-name {{ font-size: 13px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; }}
  .bar-name.highlight {{ font-weight: 600; color: var(--accent); }}
  .bar-track {{ background: var(--ground); border-radius: 3px; height: 18px; overflow: hidden; }}
  .bar-fill {{ height: 100%; border-radius: 3px; min-width: 2px; }}
  .bar-val {{ font-family: var(--ff-m); font-size: 12px; color: var(--mid); text-align: right; font-variant-numeric: tabular-nums; }}
  .bar-val.highlight {{ color: var(--accent); font-weight: 600; }}
  .callout-note {{
    background: #fffbeb; border-left: 3px solid var(--warn);
    padding: 12px 16px; border-radius: 0 4px 4px 0;
    font-size: 13px; color: #78350f; margin-top: 20px;
  }}
  .callout-note code {{
    font-family: var(--ff-m); font-size: 12px;
    background: rgba(0,0,0,.06); padding: 1px 5px; border-radius: 3px;
  }}
  .footer {{
    margin-top: 28px; font-family: var(--ff-m); font-size: 11px;
    color: var(--mid); text-align: center; letter-spacing: .06em;
  }}
  @media (max-width: 600px) {{
    .kpi-strip {{ grid-template-columns: 1fr; }}
    .kpi {{ border-right: none; border-bottom: 1px solid var(--border); }}
    .bar-row {{ grid-template-columns: 110px 1fr 60px; }}
    .header, .section {{ padding: 20px; }}
  }}
</style>
</head>
<body>
<div class="page">

  <div class="header">
    <div class="header-eyebrow">workflow-034 &middot; AMECC Theme 1 &middot; Polymer MD</div>
    <h1>{build.get('resin_name', resin)} &mdash; {temp_c:.0f}&thinsp;&deg;C / {fiber_wt:.0f}&thinsp;wt% GF<br>Lightweight Plastics for EV Applications</h1>
    <div class="badges">
      <span class="badge primary">{resin}</span>
      <span class="badge">{temp_c:.0f} &deg;C</span>
      <span class="badge">{fiber_wt:.0f} wt% GF</span>
      <span class="badge">GAFF2</span>
      <span class="badge">GROMACS 2023.2</span>
    </div>
    <div class="header-date">
      Generated {datetime.date.today()} &middot;
      {n_chains} chains &middot; 5-mer oligomer &middot; {prod_time} ps production
    </div>
  </div>

  <div class="kpi-strip">
    <div class="kpi">
      <div class="kpi-label">Young's Modulus</div>
      <div class="kpi-value">{fmt(E_corr, 3)}</div>
      <div class="kpi-unit">GPa</div>
    </div>
    <div class="kpi">
      <div class="kpi-label">Density</div>
      <div class="kpi-value">{fmt(rho_corr, 1)}</div>
      <div class="kpi-unit">kg / m&sup3;</div>
    </div>
    <div class="kpi">
      <div class="kpi-label">Specific Modulus</div>
      <div class="kpi-value">{fmt(spec_mod, 1)}</div>
      <div class="kpi-unit">kN &middot; m / kg</div>
    </div>
  </div>

  <div class="section">
    <div class="section-label">Simulation Setup</div>
    <table class="data-table">
      <tr><th>Parameter</th><th>Value</th></tr>
      <tr><td>Resin</td><td class="val">{build.get('resin_name', resin)} ({resin})</td></tr>
      <tr><td>Temperature</td><td class="val">{temp_c:.0f} &deg;C ({temp_c+273.15:.2f} K)</td></tr>
      <tr><td>Glass-fiber loading</td><td class="val">{fiber_wt:.0f} wt% &mdash; V<sub>f</sub> = {fiber_vf:.3f}</td></tr>
      <tr><td>Cell composition</td><td class="val">{n_chains} &times; 5-mer oligomer</td></tr>
      <tr><td>Force field</td><td class="val">GAFF2 &mdash; antechamber AM1-BCC charges</td></tr>
      <tr><td>Melt equilibration</td><td class="val">200 ps NPT @ {melt_temp+273.15:.0f} K (Berendsen)</td></tr>
      <tr><td>Target equilibration</td><td class="val">500 ps NPT @ {temp_c+273.15:.2f} K (Berendsen)</td></tr>
      <tr><td>Production NPT</td><td class="val">{prod_time} ps</td></tr>
      <tr><td>Deformation NVT</td><td class="val">{deform_time} ps &mdash; uniaxial, 2% strain target</td></tr>
    </table>
  </div>

  <div class="section">
    <div class="section-label">Mechanical &amp; Physical Properties</div>
    <table class="data-table">
      <tr><th>Property</th><th>MD (neat resin)</th><th>After GF correction</th></tr>
      {tr("Young's modulus", fmt(E_md, 3) + " GPa", fmt(E_corr, 3) + " GPa")}
      {tr("Density", fmt(rho_md, 1) + " kg/m&sup3;", fmt(rho_corr, 1) + " kg/m&sup3;")}
      {tr("Specific modulus", "&mdash;", fmt(spec_mod, 1) + " kN&thinsp;&middot;&thinsp;m/kg")}
    </table>
  </div>

{gf_section}

  <div class="section">
    <div class="section-label">Density Convergence &mdash; Production NPT</div>
    <div class="spark-meta">
      <span>Initial: {rho_init:.1f} kg/m&sup3;</span>
      <span>Final (last-third avg): {rho_final:.1f} kg/m&sup3;</span>
    </div>
    <div class="spark-wrap">{sparkline}</div>
    <div class="callout-note">
      &#9888; Default run times start from a 60%-packed amorphous cell; full density convergence
      typically needs 5&ndash;10 ns. Increase <code>PARAM_EQUIL_TIME_PS</code> and
      <code>PARAM_PROD_TIME_PS</code> for production-quality results.
    </div>
  </div>

  <div class="section">
    <div class="section-label">EV Lightweighting &mdash; Specific Modulus</div>
    <p style="font-size:13px;color:var(--mid);margin-bottom:18px;">
      E / &rho; in kN&thinsp;&middot;&thinsp;m/kg &mdash;
      higher = stiffer per unit mass &rarr; lighter part for the same stiffness.
    </p>
    <div class="bar-chart">
{bars_html}
    </div>
    <p style="font-size:12px;color:var(--mid);margin-top:14px;font-family:var(--ff-m);">
      * This-run value reflects short equilibration; converges toward literature with longer runs.
    </p>
  </div>

  <div class="footer">
    workflow-034 &middot; AMECC Theme 1 &middot; {datetime.date.today()}
  </div>

</div>
</body>
</html>
"""

pathlib.Path("outputs/report.html").write_text(html)

summary = {
    "resin_type":                   resin,
    "temperature_c":                temp_c,
    "fiber_loading_wt_pct":         fiber_wt,
    "youngs_modulus_gpa":           E_corr,
    "density_kg_m3":                rho_corr,
    "specific_modulus_kNm_kg":      spec_mod,
    "md_density_converged_kg_m3":   round(rho_final, 1),
    "gf_correction_applied":        gf_applied,
    "simulation_note":              note,
}
pathlib.Path("outputs/summary.json").write_text(json.dumps(summary, indent=2))

print(f"  Report  -> outputs/report.html")
print(f"  Summary -> outputs/summary.json")
print(f"  Specific modulus: {fmt(spec_mod, 1)} kN*m/kg")
