"""Generate HTML barrier performance report."""
import json, math, pathlib, datetime

diff   = json.loads(pathlib.Path("inputs/diffusion_report.json").read_text())
build  = json.loads(pathlib.Path("inputs/build_report.json").read_text())

resin     = diff["resin_type"]
penetrant = diff["penetrant"]
temp_c    = diff["temperature_c"]
D_val     = diff.get("D_cm2_s")
D_sci     = diff.get("D_cm2_s_sci", "N/A")
D_lit     = diff.get("D_literature_cm2_s")
D_lit_sci = diff.get("D_lit_sci", "N/A")
slope     = diff.get("log_log_slope")
regime    = diff.get("diffusive_regime", "")
sim_time  = diff.get("sim_time_ps", 0)

# ── MSD curve ─────────────────────────────────────────────────────────────────
xvg_lines = [l for l in pathlib.Path("inputs/msd.xvg").read_text().splitlines()
             if l and not l.startswith(("#", "@"))]
msd_pts = [(float(l.split()[0]), float(l.split()[1])) for l in xvg_lines]

# Downsample to ~120 points for SVG
step    = max(len(msd_pts) // 120, 1)
msd_ds  = msd_pts[::step]
t_vals  = [p[0] for p in msd_ds]
m_vals  = [p[1] for p in msd_ds]

# Log-log SVG (skip zeros)
ll_pairs = [(math.log10(t), math.log10(m)) for t, m in zip(t_vals, m_vals) if t > 0 and m > 0]
SVG_W, SVG_H = 680, 120

if ll_pairs:
    lx_min, lx_max = min(p[0] for p in ll_pairs), max(p[0] for p in ll_pairs)
    ly_min, ly_max = min(p[1] for p in ll_pairs), max(p[1] for p in ll_pairs)
    lx_rng = max(lx_max - lx_min, 0.1)
    ly_rng = max(ly_max - ly_min, 0.1)
    PAD = 8

    def spx(lx): return int((lx - lx_min) / lx_rng * (SVG_W - 2*PAD) + PAD)
    def spy(ly): return int(SVG_H - PAD - (ly - ly_min) / ly_rng * (SVG_H - 2*PAD))

    pts  = " ".join(f"{spx(lx)},{spy(ly)}" for lx, ly in ll_pairs)
    area = f"0,{SVG_H} " + pts + f" {SVG_W},{SVG_H}"

    # Slope=1 reference line (diffusive regime)
    ref_lx0, ref_lx1 = lx_min, lx_max
    ref_ly0 = (ly_min + ly_max) / 2 - (lx_max - lx_min) / 2
    ref_ly1 = ref_ly0 + (ref_lx1 - ref_lx0)
    ref_line = (f'<line x1="{spx(ref_lx0)}" y1="{spy(ref_ly0)}" '
                f'x2="{spx(ref_lx1)}" y2="{spy(ref_ly1)}" '
                f'stroke="#d97706" stroke-width="1.5" stroke-dasharray="6,4"/>')

    sparkline = f"""
<svg width="100%" height="{SVG_H}" viewBox="0 0 {SVG_W} {SVG_H}"
     preserveAspectRatio="none" xmlns="http://www.w3.org/2000/svg">
  <line x1="{PAD}" y1="{int(SVG_H*0.33)}" x2="{SVG_W-PAD}" y2="{int(SVG_H*0.33)}"
        stroke="#d8e0ef" stroke-width="1"/>
  <line x1="{PAD}" y1="{int(SVG_H*0.66)}" x2="{SVG_W-PAD}" y2="{int(SVG_H*0.66)}"
        stroke="#d8e0ef" stroke-width="1"/>
  {ref_line}
  <polygon points="{area}" fill="rgba(26,92,255,0.07)"/>
  <polyline points="{pts}" fill="none" stroke="#1a5cff"
            stroke-width="2" stroke-linejoin="round"/>
  <circle cx="{spx(ll_pairs[0][0])}" cy="{spy(ll_pairs[0][1])}" r="3.5" fill="#1a5cff"/>
  <circle cx="{spx(ll_pairs[-1][0])}" cy="{spy(ll_pairs[-1][1])}" r="3.5" fill="#1a5cff"/>
</svg>"""
else:
    sparkline = "<p style='color:#6b7fa3;font-size:13px'>No MSD data available.</p>"

# ── Barrier ranking (literature D values) ─────────────────────────────────────
LIT_D = {
    ("EVOH", "O2"):  2.0e-13, ("EVOH", "H2O"): 1.5e-11,
    ("PET",  "O2"):  3.4e-10, ("PET",  "H2O"): 5.0e-12,
    ("PA6",  "O2"):  1.5e-9,  ("PA6",  "H2O"): 5.0e-11,
    ("PP",   "O2"):  2.0e-8,  ("PP",   "H2O"): 3.0e-10,
    ("LDPE", "O2"):  4.5e-7,  ("LDPE", "H2O"): 1.5e-8,
}
barrier_data = sorted(
    [(r, LIT_D.get((r, penetrant), 0)) for r in ["EVOH","PET","PA6","PP","LDPE"]],
    key=lambda x: x[1]
)
D_max = max(v for _, v in barrier_data if v > 0)

def bar_row(r, d_val, highlight=False):
    if d_val <= 0:
        return ""
    log_frac = (math.log10(D_max) - math.log10(d_val)) / max(math.log10(D_max) - math.log10(1e-14), 1)
    pct   = max(2, min(98, log_frac * 96 + 2))
    color = "#1a5cff" if highlight else "#b0bdd8"
    n_cls = "bar-name highlight" if highlight else "bar-name"
    v_cls = "bar-val highlight"  if highlight else "bar-val"
    tag   = " ← best barrier" if r == "EVOH" else ""
    return f"""      <div class="bar-row">
        <div class="{n_cls}">{r}{tag}</div>
        <div class="bar-track"><div class="bar-fill" style="width:{pct:.1f}%;background:{color}"></div></div>
        <div class="{v_cls}">{d_val:.1e}</div>
      </div>"""

bars_html = "\n".join(bar_row(r, d, r==resin) for r, d in barrier_data)

# ── Regime badge ──────────────────────────────────────────────────────────────
if "diffusive" in regime and "sub" not in regime:
    regime_badge = f'<span class="badge good">&#10003; {regime}</span>'
elif "sub" in regime:
    regime_badge = f'<span class="badge warn">&#9888; {regime}</span>'
else:
    regime_badge = f'<span class="badge mid">{regime}</span>'

# ── HTML ──────────────────────────────────────────────────────────────────────
pen_label = "O₂" if penetrant == "O2" else "H₂O"

html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>Barrier Film Report — {resin} / {pen_label} @ {temp_c:.0f} °C</title>
<style>
  *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
  :root {{
    --ink: #0c1829; --ground: #f5f7fb; --surface: #ffffff;
    --accent: #1a5cff; --accent-dim: #eef2ff; --mid: #6b7fa3;
    --border: #d8e0ef; --warn: #d97706; --good: #16a34a;
    --ff-d: Georgia, 'Times New Roman', serif;
    --ff-b: system-ui, -apple-system, sans-serif;
    --ff-m: ui-monospace, 'Cascadia Code', 'Courier New', monospace;
  }}
  body {{ font-family: var(--ff-b); font-size: 15px; line-height: 1.65;
          color: var(--ink); background: var(--ground); padding: 40px 20px 80px; }}
  .page {{ max-width: 800px; margin: 0 auto; }}
  .header {{ background: var(--ink); border-radius: 10px 10px 0 0;
             padding: 32px 36px 28px; position: relative; overflow: hidden; }}
  .header::before {{ content: ''; position: absolute; inset: 0;
    background-image: radial-gradient(circle, rgba(26,92,255,.18) 1px, transparent 1px);
    background-size: 28px 28px; pointer-events: none; }}
  .header-eyebrow {{ font-family: var(--ff-m); font-size: 11px; letter-spacing: .12em;
    text-transform: uppercase; color: var(--accent); margin-bottom: 10px; }}
  .header h1 {{ font-family: var(--ff-d); font-size: 1.75rem; font-weight: normal;
    color: #fff; line-height: 1.25; text-wrap: balance; margin-bottom: 16px; }}
  .badges {{ display: flex; gap: 8px; flex-wrap: wrap; }}
  .badge {{ font-family: var(--ff-m); font-size: 11px; letter-spacing: .06em;
    padding: 3px 10px; border-radius: 3px;
    background: rgba(26,92,255,.25); color: #a8c0ff;
    border: 1px solid rgba(26,92,255,.4); }}
  .badge.primary {{ background: var(--accent); color: #fff; border-color: var(--accent); }}
  .badge.good {{ background: rgba(22,163,74,.2); color: #16a34a; border-color: rgba(22,163,74,.4); }}
  .badge.warn {{ background: rgba(217,119,6,.15); color: var(--warn); border-color: rgba(217,119,6,.4); }}
  .badge.mid  {{ background: rgba(107,127,163,.15); color: var(--mid); border-color: rgba(107,127,163,.4); }}
  .header-date {{ font-family: var(--ff-m); font-size: 11px; color: var(--mid); margin-top: 14px; }}
  .kpi-strip {{ display: grid; grid-template-columns: repeat(3, 1fr);
    background: var(--surface); border: 1px solid var(--border); border-top: none; }}
  .kpi {{ padding: 24px 28px; border-right: 1px solid var(--border); }}
  .kpi:last-child {{ border-right: none; }}
  .kpi-label {{ font-family: var(--ff-m); font-size: 10.5px; letter-spacing: .1em;
    text-transform: uppercase; color: var(--mid); margin-bottom: 6px; }}
  .kpi-value {{ font-family: var(--ff-d); font-size: 1.9rem; font-weight: bold;
    color: var(--ink); font-variant-numeric: tabular-nums; line-height: 1; }}
  .kpi-unit {{ font-family: var(--ff-m); font-size: 11px; color: var(--mid); margin-top: 4px; }}
  .section {{ background: var(--surface); border: 1px solid var(--border);
    border-top: none; padding: 28px 36px; }}
  .section:last-of-type {{ border-radius: 0 0 10px 10px; }}
  .section-label {{ font-family: var(--ff-m); font-size: 10.5px; letter-spacing: .12em;
    text-transform: uppercase; color: var(--accent);
    display: flex; align-items: center; gap: 10px; margin-bottom: 16px; }}
  .section-label::after {{ content: ''; flex: 1; height: 1px; background: var(--accent-dim); }}
  .data-table {{ width: 100%; border-collapse: collapse; }}
  .data-table th, .data-table td {{ text-align: left; padding: 8px 12px;
    border-bottom: 1px solid var(--border); font-size: 13.5px; }}
  .data-table th {{ font-family: var(--ff-m); font-size: 11px; letter-spacing: .05em;
    color: var(--mid); background: var(--ground); font-weight: 500; }}
  .data-table .val {{ font-family: var(--ff-m); font-size: 13px; font-variant-numeric: tabular-nums; }}
  .data-table tr:last-child td {{ border-bottom: none; }}
  .spark-wrap {{ background: var(--ground); border: 1px solid var(--border);
    border-radius: 6px; padding: 16px; overflow-x: auto; }}
  .spark-meta {{ display: flex; justify-content: space-between;
    font-family: var(--ff-m); font-size: 11px; color: var(--mid); margin-bottom: 8px; }}
  .bar-chart {{ display: flex; flex-direction: column; gap: 10px; }}
  .bar-row {{ display: grid; grid-template-columns: 130px 1fr 90px; align-items: center; gap: 10px; }}
  .bar-name {{ font-size: 13px; white-space: nowrap; }}
  .bar-name.highlight {{ font-weight: 600; color: var(--accent); }}
  .bar-track {{ background: var(--ground); border-radius: 3px; height: 18px; overflow: hidden; }}
  .bar-fill {{ height: 100%; border-radius: 3px; min-width: 2px; }}
  .bar-val {{ font-family: var(--ff-m); font-size: 12px; color: var(--mid);
    text-align: right; font-variant-numeric: tabular-nums; }}
  .bar-val.highlight {{ color: var(--accent); font-weight: 600; }}
  .callout-note {{ background: #fffbeb; border-left: 3px solid var(--warn);
    padding: 12px 16px; border-radius: 0 4px 4px 0;
    font-size: 13px; color: #78350f; margin-top: 20px; }}
  .callout-note code {{ font-family: var(--ff-m); font-size: 12px;
    background: rgba(0,0,0,.06); padding: 1px 5px; border-radius: 3px; }}
  .footer {{ margin-top: 28px; font-family: var(--ff-m); font-size: 11px;
    color: var(--mid); text-align: center; letter-spacing: .06em; }}
  @media (max-width: 600px) {{
    .kpi-strip {{ grid-template-columns: 1fr; }}
    .kpi {{ border-right: none; border-bottom: 1px solid var(--border); }}
    .bar-row {{ grid-template-columns: 80px 1fr 70px; }}
    .header, .section {{ padding: 20px; }}
  }}
</style>
</head>
<body>
<div class="page">

  <div class="header">
    <div class="header-eyebrow">workflow-031 &middot; AMECC Theme 2 &middot; Barrier Films</div>
    <h1>{build.get('resin_name', resin)} &mdash; {pen_label} @ {temp_c:.0f}&thinsp;&deg;C<br>Plastics that Protect Food</h1>
    <div class="badges">
      <span class="badge primary">{resin}</span>
      <span class="badge">{pen_label} penetrant</span>
      <span class="badge">{temp_c:.0f} &deg;C</span>
      <span class="badge">GAFF2</span>
      <span class="badge">GROMACS 2023.2</span>
      {regime_badge}
    </div>
    <div class="header-date">
      Generated {datetime.date.today()} &middot; {sim_time:.0f} ps simulation
    </div>
  </div>

  <div class="kpi-strip">
    <div class="kpi">
      <div class="kpi-label">Diffusion Coeff. D</div>
      <div class="kpi-value" style="font-size:1.5rem">{D_sci}</div>
      <div class="kpi-unit">cm&sup2; / s (MD)</div>
    </div>
    <div class="kpi">
      <div class="kpi-label">Literature Ref.</div>
      <div class="kpi-value" style="font-size:1.5rem">{D_lit_sci}</div>
      <div class="kpi-unit">cm&sup2; / s (lit.)</div>
    </div>
    <div class="kpi">
      <div class="kpi-label">MSD slope (log-log)</div>
      <div class="kpi-value">{f"{slope:.2f}" if slope else "N/A"}</div>
      <div class="kpi-unit">1.0 = diffusive regime</div>
    </div>
  </div>

  <div class="section">
    <div class="section-label">Simulation Setup</div>
    <table class="data-table">
      <tr><th>Parameter</th><th>Value</th></tr>
      <tr><td>Resin</td><td class="val">{build.get('resin_name', resin)} ({resin})</td></tr>
      <tr><td>Penetrant</td><td class="val">{pen_label} &mdash; {'TraPPE 2-site' if penetrant=='O2' else 'SPC/E'}</td></tr>
      <tr><td>Temperature</td><td class="val">{temp_c:.0f} &deg;C ({temp_c+273.15:.2f} K)</td></tr>
      <tr><td>Simulation length</td><td class="val">{sim_time:.0f} ps</td></tr>
      <tr><td>MSD regime</td><td class="val">{regime}</td></tr>
    </table>
  </div>

  <div class="section">
    <div class="section-label">MSD vs Time (log-log) &mdash; Diffusive regime: slope = 1</div>
    <div class="spark-meta">
      <span>Time: {msd_pts[0][0]:.1f} &rarr; {msd_pts[-1][0]:.1f} ps</span>
      <span style="color:var(--warn)">&#9472; &#9472; slope = 1 reference</span>
    </div>
    <div class="spark-wrap">{sparkline}</div>
    {'<div class="callout-note">&#9888; MSD slope = ' + str(round(slope,2)) + ' — not yet in the diffusive regime. Increase <code>PARAM_DIFF_TIME_PS</code> (try 20 000 ps for EVOH). D value is an underestimate.</div>' if slope and slope < 0.8 else ''}
  </div>

  <div class="section">
    <div class="section-label">Barrier Ranking &mdash; {pen_label} Diffusion Coefficient (literature, cm&sup2;/s)</div>
    <p style="font-size:13px;color:var(--mid);margin-bottom:18px;">
      Lower D = better barrier. Bars show relative barrier performance on a log scale (longer bar = lower D = better barrier).
    </p>
    <div class="bar-chart">
{bars_html}
    </div>
    <p style="font-size:12px;color:var(--mid);margin-top:14px;font-family:var(--ff-m);">
      Literature values at ~25 °C. Your simulation result: {D_sci} cm&sup2;/s.
    </p>
  </div>

  <div class="footer">
    workflow-031 &middot; AMECC Theme 2 &middot; {datetime.date.today()}
  </div>

</div>
</body>
</html>
"""

pathlib.Path("outputs/report.html").write_text(html)

summary = {
    "resin_type":       resin,
    "penetrant":        penetrant,
    "temperature_c":    temp_c,
    "D_cm2_s":          D_val,
    "D_cm2_s_sci":      D_sci,
    "D_literature":     D_lit,
    "log_log_slope":    slope,
    "diffusive_regime": regime,
    "sim_time_ps":      sim_time,
}
pathlib.Path("outputs/summary.json").write_text(json.dumps(summary, indent=2))

print(f"  Report  -> outputs/report.html")
print(f"  Summary -> outputs/summary.json")
print(f"  D = {D_sci} cm²/s  (lit: {D_lit_sci})")
