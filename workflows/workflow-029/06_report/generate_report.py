#!/usr/bin/env python3
"""
Binder design campaign dashboard.

Reads inputs from ./inputs/ with hardcoded paths (no argparse).
Writes a single self-contained HTML report to ./outputs/report/index.html.
All data embedded inline as JavaScript variables.
Libraries loaded from CDN: Bootstrap 5.3.3, Plotly 2.35.2, Molstar latest.

Zero-results handling: renders a banner without crashing when
prodigy_all_designs.json is empty or top10/ has no PDBs.
"""
import glob
import json
import os

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.offline as pyo
from plotly.subplots import make_subplots

os.makedirs('./outputs/report', exist_ok=True)

top_n = int(os.environ.get('PARAM_TOP_N', '10'))
min_iptm = float(os.environ.get('PARAM_MIN_IPTM', '0.6'))
min_plddt_binder = float(os.environ.get('PARAM_MIN_PLDDT_BINDER', '80.0'))
max_bb_rmsd = float(os.environ.get('PARAM_MAX_BB_RMSD', '1.5'))
max_pae_interaction = float(os.environ.get('PARAM_MAX_PAE_INTERACTION', '10.0'))

# ── Load inputs ───────────────────────────────────────────────────────────────

with open('./inputs/folded/filter_report.json') as f:
    filter_report = json.load(f)

with open('./inputs/results/prodigy_all_designs.json') as f:
    prodigy_results = json.load(f)

# Glob top10 PDB files — do not assume exactly top_n files exist
top10_pdbs = sorted(glob.glob(f'./inputs/results/top{top_n}/*.pdb'))

print(f'Filter report: {len(filter_report)} entries', flush=True)
print(f'PRODIGY results: {len(prodigy_results)} entries', flush=True)
print(f'Top-N PDB files found: {len(top10_pdbs)}', flush=True)

no_results = len(prodigy_results) == 0

# ── Design funnel numbers ─────────────────────────────────────────────────────

n_backbones = len({r['design_id'].rsplit('_seq_', 1)[0] for r in filter_report})
n_sequences = len(filter_report)
n_passed_filter = sum(1 for r in filter_report if r['pass'])
n_passed_prodigy = sum(1 for r in prodigy_results if not r.get('weak_binder', True))

# ── Colours ───────────────────────────────────────────────────────────────────

COLORS = {
    'excellent': '#0f766e',
    'good':      '#0369a1',
    'moderate':  '#fb7c3c',
    'poor':      '#dc2626',
    'pass':      '#0369a1',
    'fail':      '#dc2626',
}


def iptm_quality(v):
    if v is None:
        return 'poor'
    if v >= 0.8:
        return 'excellent'
    if v >= 0.6:
        return 'good'
    if v >= 0.4:
        return 'moderate'
    return 'poor'


def fmt_kd(kd):
    if kd is None:
        return '—'
    coeff, power = f'{kd:.2e}'.split('e')
    return f'{coeff} &times; 10<sup>{int(power)}</sup>'


# ── Funnel chart ─────────────────────────────────────────────────────────────

def make_funnel():
    labels = ['Backbones generated', 'Sequences designed', 'Passed fold filter', 'Below ΔG cutoff']
    values = [n_backbones, n_sequences, n_passed_filter, n_passed_prodigy]
    fig = go.Figure(go.Funnel(
        y=labels, x=values,
        textinfo='value+percent initial',
        marker=dict(color=['#075985', '#0369a1', '#0f766e', '#0d9488']),
    ))
    fig.update_layout(
        title=dict(text='Design Campaign Funnel', x=0.5, font=dict(size=18)),
        height=350,
        template='plotly_white',
        margin=dict(l=20, r=20, t=60, b=20),
    )
    return pyo.plot(fig, output_type='div', include_plotlyjs=False)


# ── iPTM vs ΔG scatter ────────────────────────────────────────────────────────

def make_scatter():
    if not prodigy_results:
        return '<p class="text-muted text-center py-4">No data to display.</p>'

    df = pd.DataFrame(prodigy_results)
    colors = [COLORS['good'] if not r else COLORS['poor'] for r in df.get('weak_binder', [False] * len(df))]

    fig = go.Figure(go.Scatter(
        x=df.get('iptm', []),
        y=df.get('dg', []),
        mode='markers',
        marker=dict(size=10, color=colors, line=dict(width=1, color='white')),
        text=df.get('design_id', []),
        hovertemplate=(
            '<b>%{text}</b><br>'
            'iPTM: %{x:.3f}<br>'
            'ΔG: %{y:.2f} kcal/mol<extra></extra>'
        ),
    ))
    fig.update_layout(
        title=dict(text='iPTM vs ΔG', x=0.5, font=dict(size=18)),
        xaxis=dict(title='iPTM', tickformat='.2f'),
        yaxis=dict(title='ΔG (kcal/mol)'),
        height=420,
        template='plotly_white',
        margin=dict(l=60, r=20, t=60, b=60),
    )
    return pyo.plot(fig, output_type='div', include_plotlyjs=False)


# ── Ranking bar chart ─────────────────────────────────────────────────────────

def make_ranking_chart():
    if not prodigy_results:
        return '<p class="text-muted text-center py-4">No data to display.</p>'

    display = prodigy_results[:20]
    labels = [r['design_id'] for r in display]
    dg_vals = [r['dg'] for r in display]
    colors = [COLORS['good'] if not r.get('weak_binder') else COLORS['moderate'] for r in display]

    fig = go.Figure(go.Bar(
        x=dg_vals, y=labels,
        orientation='h',
        marker_color=colors,
        text=[f'{v:.2f}' for v in dg_vals],
        textposition='inside',
        textfont=dict(color='white', size=11),
    ))
    fig.update_layout(
        title=dict(text=f'Top {len(display)} Designs by ΔG', x=0.5, font=dict(size=18)),
        xaxis=dict(title='ΔG (kcal/mol)'),
        yaxis=dict(categoryorder='total descending'),
        height=max(350, len(display) * 25),
        template='plotly_white',
        margin=dict(l=160, r=20, t=60, b=60),
    )
    return pyo.plot(fig, output_type='div', include_plotlyjs=False)


# ── 3D viewer section (Molstar) ──────────────────────────────────────────────

def make_viewer_section():
    if not top10_pdbs:
        return '<p class="text-muted text-center py-4">No structure files available for 3D viewing.</p>'

    structures = []
    for pdb_path in top10_pdbs[:top_n]:
        design_id = os.path.splitext(os.path.basename(pdb_path))[0].replace('_rank001', '')
        try:
            with open(pdb_path, encoding='utf-8') as fh:
                raw = fh.read()
            escaped = (raw.replace('\\', '\\\\')
                          .replace('`', '\\`')
                          .replace('${', '\\${')
                          .replace('</script>', '<\\/script>'))
            structures.append({'id': design_id, 'data': escaped})
        except Exception as e:
            print(f'WARNING: Could not read {pdb_path}: {e}', flush=True)

    if not structures:
        return '<p class="text-muted text-center py-4">No structure files could be loaded.</p>'

    options_html = ''.join(
        f'<option value="{i}"{" selected" if i == 0 else ""}>{s["id"]}</option>\n'
        for i, s in enumerate(structures)
    )
    struct_entries = ',\n'.join(
        f'    {{label: {json.dumps(s["id"])}, data: `{s["data"]}`}}'
        for s in structures
    )

    return f"""
        <select id="structureSelect" class="form-select form-select-sm mb-3" style="max-width:400px;"
                onchange="loadStructure(this.value)">
            {options_html}
        </select>
        <div id="viewer3d" style="height:500px;width:100%;position:relative;border-radius:10px;"></div>
        <script>
        (function(){{
            var structs = [
{struct_entries}
            ];
            var molPlugin = null;

            window.loadStructure = function(idx) {{
                if (!molPlugin || !structs[idx]) return;
                molPlugin.clear();
                var s = structs[idx];
                molPlugin.builders.data.rawData({{data: s.data, label: s.label}})
                    .then(function(d) {{ return molPlugin.builders.structure.parseTrajectory(d, 'pdb'); }})
                    .then(function(t) {{ return molPlugin.builders.structure.hierarchy.applyPreset(t, 'default'); }})
                    .catch(function(e) {{ console.error('Mol* load error:', e); }});
            }};

            molstar.Viewer.create('viewer3d', {{
                layoutIsExpanded: false,
                layoutShowControls: false,
                layoutShowLeftPanel: false,
                layoutShowSequence: false,
                layoutShowLog: false,
                layoutShowRemoteState: false,
                viewportShowAnimation: false,
                viewportShowExpand: true,
                viewportShowSelectionMode: false
            }}).then(function(viewer) {{
                molPlugin = viewer.plugin;
                loadStructure(0);
            }}).catch(function(e) {{
                document.getElementById('viewer3d').innerHTML =
                    '<p style="color:red;padding:16px;">Mol* failed to initialize: ' + e + '</p>';
            }});
        }})();
        </script>
"""


# ── Ranking table rows ────────────────────────────────────────────────────────

def make_table_rows():
    if not prodigy_results:
        return '<tr><td colspan="8" class="text-center text-muted">No results</td></tr>'
    rows = []
    for i, r in enumerate(prodigy_results, 1):
        weak = r.get('weak_binder', False)
        badge = (f'<span class="badge" style="background:{COLORS["poor"]};color:white;'
                 f'border-radius:12px;padding:3px 8px;font-size:11px;">Weak</span>'
                 if weak else
                 f'<span class="badge" style="background:{COLORS["good"]};color:white;'
                 f'border-radius:12px;padding:3px 8px;font-size:11px;">Strong</span>')
        iptm_val = r.get('iptm')
        plddt_val = r.get('plddt_binder')
        rmsd_val = r.get('bb_rmsd')
        pae_val = r.get('pae_interaction')
        rows.append(f"""
            <tr>
                <td>{i}</td>
                <td><strong>{r['design_id']}</strong></td>
                <td>{r['dg']:.3f}</td>
                <td>{fmt_kd(r.get('kd'))}</td>
                <td>{'—' if iptm_val is None else f'{iptm_val:.3f}'}</td>
                <td>{'—' if plddt_val is None else f'{plddt_val:.1f}'}</td>
                <td>{'—' if rmsd_val is None else f'{rmsd_val:.2f}'}</td>
                <td>{'—' if pae_val is None else f'{pae_val:.1f}'}</td>
                <td>{badge}</td>
            </tr>""")
    return '\n'.join(rows)


# ── Assemble HTML ─────────────────────────────────────────────────────────────

funnel_div = make_funnel()
scatter_div = make_scatter()
ranking_div = make_ranking_chart()
viewer_html = make_viewer_section()
table_rows = make_table_rows()

zero_results_banner = ''
if no_results:
    zero_results_banner = """
        <div class="alert alert-warning" role="alert" style="font-size:1.1rem;">
            <strong>No designs passed all filters.</strong>
            Consider relaxing thresholds: increase <code>max_bb_rmsd</code>, decrease <code>min_iptm</code>,
            or generate more backbones by increasing <code>num_designs</code>.
        </div>"""

best_dg = f'{prodigy_results[0]["dg"]:.2f}' if prodigy_results else '—'
best_iptm = f'{max((r.get("iptm") or 0) for r in prodigy_results):.3f}' if prodigy_results else '—'
n_strong = sum(1 for r in prodigy_results if not r.get('weak_binder', True))

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Binder Design Campaign Report</title>
    <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
    <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/molstar@5.9.0/build/viewer/molstar.css">
    <script src="https://cdn.jsdelivr.net/npm/molstar@5.9.0/build/viewer/molstar.js"></script>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css" rel="stylesheet">
    <style>
        :root {{
            --primary: #075985;
            --secondary: #0284c7;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: linear-gradient(135deg, var(--primary) 0%, var(--secondary) 100%);
            min-height: 100vh;
            margin: 0;
            color: #1f2937;
        }}
        .main-container {{ max-width: 1400px; margin: 0 auto; padding: 20px; }}
        .glass {{
            background: rgba(255,255,255,0.95);
            backdrop-filter: blur(10px);
            border: 1px solid rgba(255,255,255,0.2);
            border-radius: 15px;
            box-shadow: 0 8px 32px rgba(0,0,0,0.1);
            margin-bottom: 28px;
        }}
        .glass-header {{
            background: linear-gradient(135deg, var(--primary) 0%, var(--secondary) 100%);
            color: white;
            border-radius: 15px 15px 0 0;
            padding: 18px 24px;
            font-weight: 600;
            font-size: 1.05rem;
        }}
        .glass-body {{ padding: 24px; }}
        .stat-card {{
            background: rgba(255,255,255,0.95);
            border-radius: 14px;
            padding: 22px;
            text-align: center;
            box-shadow: 0 4px 16px rgba(0,0,0,0.08);
            margin-bottom: 20px;
        }}
        .stat-value {{ font-size: 2rem; font-weight: 700; color: var(--primary); }}
        .stat-label {{ color: #64748b; font-size: 0.85rem; text-transform: uppercase; letter-spacing: 1px; }}
        .table th {{ background: #f8fafc; font-weight: 600; color: var(--primary); border-top: none; }}
        .table tbody tr:hover {{ background: rgba(7,89,133,0.04); }}
        .plot-wrap {{ background: white; border-radius: 10px; padding: 16px; }}
    </style>
</head>
<body>
<div class="main-container">

    <!-- Header -->
    <div class="glass">
        <div class="glass-body">
            <h1 style="color:var(--primary);font-size:2.2rem;font-weight:700;margin:0 0 8px;">
                <i class="fas fa-dna"></i> Binder Design Campaign Report
            </h1>
            <p style="color:#64748b;margin:0;">
                RFdiffusion → ProteinMPNN → ColabFold → PRODIGY &nbsp;|&nbsp;
                <strong>{n_backbones}</strong> backbones &nbsp;·&nbsp;
                <strong>{n_sequences}</strong> sequences &nbsp;·&nbsp;
                <strong>{n_passed_filter}</strong> passed filter &nbsp;·&nbsp;
                <strong>{len(prodigy_results)}</strong> scored
            </p>
        </div>
    </div>

    {zero_results_banner}

    <!-- Summary cards -->
    <div class="row">
        <div class="col-6 col-md-3">
            <div class="stat-card">
                <div style="font-size:2rem;margin-bottom:8px;">🏆</div>
                <div class="stat-value">{best_dg}</div>
                <div class="stat-label">Best ΔG (kcal/mol)</div>
            </div>
        </div>
        <div class="col-6 col-md-3">
            <div class="stat-card">
                <div style="font-size:2rem;margin-bottom:8px;">📈</div>
                <div class="stat-value">{best_iptm}</div>
                <div class="stat-label">Best iPTM</div>
            </div>
        </div>
        <div class="col-6 col-md-3">
            <div class="stat-card">
                <div style="font-size:2rem;margin-bottom:8px;">✅</div>
                <div class="stat-value">{n_passed_filter}</div>
                <div class="stat-label">Passed fold filter</div>
            </div>
        </div>
        <div class="col-6 col-md-3">
            <div class="stat-card">
                <div style="font-size:2rem;margin-bottom:8px;">💊</div>
                <div class="stat-value">{n_strong}</div>
                <div class="stat-label">Strong binders</div>
            </div>
        </div>
    </div>

    <!-- Funnel -->
    <div class="glass">
        <div class="glass-header"><i class="fas fa-filter"></i> Design Campaign Funnel</div>
        <div class="glass-body">
            <div class="plot-wrap">{funnel_div}</div>
        </div>
    </div>

    <!-- Ranking + Scatter -->
    <div class="row">
        <div class="col-lg-6">
            <div class="glass">
                <div class="glass-header"><i class="fas fa-chart-bar"></i> ΔG Ranking</div>
                <div class="glass-body">
                    <div class="plot-wrap">{ranking_div}</div>
                </div>
            </div>
        </div>
        <div class="col-lg-6">
            <div class="glass">
                <div class="glass-header"><i class="fas fa-chart-scatter"></i> iPTM vs ΔG</div>
                <div class="glass-body">
                    <div class="plot-wrap">{scatter_div}</div>
                </div>
            </div>
        </div>
    </div>

    <!-- 3D Viewer -->
    <div class="glass">
        <div class="glass-header"><i class="fas fa-cube"></i> 3D Structure Viewer (Top {top_n})</div>
        <div class="glass-body">
            {viewer_html}
        </div>
    </div>

    <!-- Full ranking table -->
    <div class="glass">
        <div class="glass-header"><i class="fas fa-table"></i> All Scored Designs</div>
        <div class="glass-body">
            <div class="table-responsive">
                <table class="table table-striped table-hover">
                    <thead>
                        <tr>
                            <th>#</th><th>Design ID</th>
                            <th>ΔG (kcal/mol)</th><th>Kd (M)</th>
                            <th>iPTM</th><th>pLDDT</th>
                            <th>RMSD (Å)</th><th>PAE</th><th>Status</th>
                        </tr>
                    </thead>
                    <tbody>{table_rows}</tbody>
                </table>
            </div>
        </div>
    </div>

    <!-- Methodology note -->
    <div class="glass">
        <div class="glass-header"><i class="fas fa-info-circle"></i> Methodology & Caveats</div>
        <div class="glass-body" style="font-size:0.9rem;color:#374151;">
            <p><strong>Pipeline:</strong> RFdiffusion (Watson et al. 2023, Nature) generates poly-Gly binder backbones.
            ProteinMPNN (SolubleMPNN mode) designs sequences for chain A while holding chain B fixed.
            ColabFold (AlphaFold2 multimer v3) folds each sequence as a binder:receptor complex.
            PRODIGY scores the predicted interface.</p>
            <p><strong>Filter thresholds applied:</strong>
            iPTM ≥ {min_iptm} &nbsp;·&nbsp;
            pLDDT (binder) ≥ {min_plddt_binder} &nbsp;·&nbsp;
            Cα RMSD ≤ {max_bb_rmsd} Å &nbsp;·&nbsp;
            Interface PAE ≤ {max_pae_interaction} Å
            </p>
            <p><strong>PRODIGY accuracy:</strong> Vangone &amp; Bonvin 2016 report r = 0.73,
            RMSE = 1.89 kcal/mol on natural crystal structures. Performance on AlphaFold2-predicted
            structures is lower and unbenchmarked. Treat ΔG values as <em>relative rankings</em>
            within this campaign, not absolute affinity predictions.</p>
            <p><strong>Recommended next steps:</strong> Validate top candidates by SPR, ITC, or
            co-crystallisation. Consider affinity maturation or MD refinement for promising hits.</p>
        </div>
    </div>

</div>
</body>
</html>"""

out_path = './outputs/report/index.html'
with open(out_path, 'w', encoding='utf-8') as f:
    f.write(html)

print(f'Report written to {out_path}', flush=True)
