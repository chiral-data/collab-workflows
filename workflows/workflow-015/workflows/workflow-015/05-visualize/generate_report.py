#!/usr/bin/env python3
"""
Node 05: Boltz-2 vs Chai-1 comparison dashboard.

Reads boltz_summary.json and chai_summary.json from prediction nodes,
generates a self-contained HTML report using Plotly.js + Bootstrap 5.

All metrics are displayed on a unified 0-1 scale:
  - Boltz pLDDT is already 0-1
  - Chai pLDDT is 0-100 in B-factor column → divided by 100 here
"""

import argparse
import json
from datetime import datetime
from pathlib import Path



# ── Parsing ───────────────────────────────────────────────────────────────────

def load_summary(path):
    """Load a summary JSON; exit with error if file is missing."""
    p = Path(path)
    if not p.exists():
        print(f"Error: {path} not found — did the prediction node run successfully?", file=__import__('sys').stderr)
        __import__('sys').exit(1)
    with open(p) as f:
        data = json.load(f)
    print(f"  Loaded: {path}")
    return data


def normalize_plddt(confidence_list, tool):
    """Normalize pLDDT to 0-1 scale. Chai stores 0-100 in B-factor."""
    for entry in confidence_list:
        if entry.get('plddt') is not None:
            if tool == 'chai1' and entry['plddt'] > 1.0:
                entry['plddt'] = round(entry['plddt'] / 100.0, 4)
    return confidence_list


def best_model(confidence_list):
    """Return the model entry with the highest pLDDT."""
    return max(confidence_list, key=lambda x: x.get('plddt', 0))


# ── Chart data helpers ────────────────────────────────────────────────────────

def metric_comparison_data(boltz_conf, chai_conf):
    """Build side-by-side metric arrays for Plotly grouped bar chart."""
    metrics = ['plddt', 'ptm', 'iptm']
    labels  = ['pLDDT', 'pTM', 'ipTM']

    b_best = best_model(boltz_conf)
    c_best = best_model(chai_conf)

    boltz_vals = [round(b_best.get(m, 0) or 0, 4) for m in metrics]
    chai_vals  = [round(c_best.get(m, 0) or 0, 4) for m in metrics]

    return labels, boltz_vals, chai_vals


def all_models_table(boltz_conf, chai_conf):
    """Build HTML table rows for all models from both tools."""
    rows = []
    for entry in boltz_conf:
        rows.append(_table_row('Boltz-2', entry, '#0284c7'))
    for entry in chai_conf:
        rows.append(_table_row('Chai-1', entry, '#7c3aed'))
    return '\n'.join(rows)


def _table_row(tool, entry, color):
    plddt = entry.get('plddt', 'N/A')
    ptm   = entry.get('ptm',   'N/A')
    iptm  = entry.get('iptm',  'N/A')
    pae   = entry.get('pae_mean', 'N/A')
    pde   = entry.get('pde_mean', 'N/A')

    def fmt(v): return f'{v:.4f}' if isinstance(v, float) else v

    return f"""
        <tr>
            <td><span style="background:{color};color:white;padding:2px 8px;
                border-radius:10px;font-size:11px;">{tool}</span></td>
            <td>{entry.get('sample', 'N/A')}</td>
            <td>{fmt(plddt)}</td>
            <td>{fmt(ptm)}</td>
            <td>{fmt(iptm)}</td>
            <td>{fmt(pae)}</td>
            <td>{fmt(pde)}</td>
        </tr>"""


# ── HTML generation ───────────────────────────────────────────────────────────

def generate_html(boltz_data, chai_data):
    boltz_conf = normalize_plddt(boltz_data['confidence'], 'boltz2')
    chai_conf  = normalize_plddt(chai_data['confidence'],  'chai1')

    labels, boltz_vals, chai_vals = metric_comparison_data(boltz_conf, chai_conf)
    b_best = best_model(boltz_conf)
    c_best = best_model(chai_conf)
    table_rows = all_models_table(boltz_conf, chai_conf)

    # Winner callout per metric
    def winner(b, c, lower_is_better=False):
        if b is None or c is None:
            return '—'
        return 'Boltz-2 ✓' if (b < c if lower_is_better else b > c) else 'Chai-1 ✓'

    mock_banner = ""

    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Boltz-2 vs Chai-1 Comparison</title>
    <script src="https://cdn.plot.ly/plotly-2.24.1.min.js"></script>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css" rel="stylesheet">
    <style>
        :root {{
            --boltz: #0284c7;
            --chai:  #7c3aed;
            --excellent: #0f766e;
            --good:      #0369a1;
            --moderate:  #fb7c3c;
            --poor:      #dc2626;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: linear-gradient(135deg, #075985 0%, #0284c7 100%);
            min-height: 100vh; margin: 0; color: #1f2937;
        }}
        .main-container {{ max-width: 1400px; margin: 0 auto; padding: 20px; }}
        .header {{
            background: rgba(255,255,255,0.95); border-radius: 15px;
            padding: 30px; margin-bottom: 30px;
            box-shadow: 0 8px 32px rgba(0,0,0,0.1);
        }}
        .header h1 {{ color: #075985; font-size: 2.2rem; font-weight: 700; margin: 0 0 8px 0; }}
        .header .subtitle {{ color: #64748b; font-size: 1rem; margin: 0; }}
        .card {{
            background: rgba(255,255,255,0.95); border-radius: 15px;
            box-shadow: 0 8px 32px rgba(0,0,0,0.1); margin-bottom: 30px;
            transition: transform 0.2s ease;
        }}
        .card:hover {{ transform: translateY(-2px); }}
        .card-header {{
            background: linear-gradient(135deg, #075985 0%, #0284c7 100%);
            color: white; border-radius: 15px 15px 0 0 !important;
            padding: 18px 25px; font-weight: 600; font-size: 1.05rem;
        }}
        .summary-grid {{
            display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px; margin-bottom: 30px;
        }}
        .summary-card {{
            background: rgba(255,255,255,0.95); border-radius: 15px;
            padding: 22px; text-align: center;
            box-shadow: 0 8px 32px rgba(0,0,0,0.1);
        }}
        .summary-card .tool-label {{
            font-size: 0.8rem; text-transform: uppercase;
            letter-spacing: 1px; margin-bottom: 6px;
        }}
        .summary-card .value {{ font-size: 1.8rem; font-weight: 700; }}
        .summary-card .metric-name {{ color: #64748b; font-size: 0.85rem; }}
        .boltz-label {{ color: var(--boltz); }}
        .chai-label  {{ color: var(--chai);  }}
        .winner-badge {{
            display: inline-block; padding: 3px 10px; border-radius: 12px;
            font-size: 0.75rem; font-weight: 600; color: white;
            background: var(--excellent);
        }}
        .plot-container {{ padding: 20px; }}
        .table th {{ background: #f8fafc; font-weight: 600; color: #075985; }}
        .table tbody tr:hover {{ background: rgba(7,89,133,0.05); }}
    </style>
</head>
<body>
<div class="main-container">

    <!-- Header -->
    <div class="header">
        <h1><i class="fas fa-balance-scale"></i> Boltz-2 vs Chai-1 Structure Prediction</h1>
        <p class="subtitle">
            <strong>Generated:</strong> {timestamp} &nbsp;|&nbsp;
            <strong>Metrics on unified 0–1 scale</strong> (lower is better for PAE/PDE in Å)
        </p>
    </div>

    {mock_banner}

    <!-- Top metric summary cards -->
    <div class="summary-grid">
        <div class="summary-card">
            <div class="tool-label boltz-label">Boltz-2</div>
            <div class="value boltz-label">{b_best.get('plddt', 0):.3f}</div>
            <div class="metric-name">pLDDT (best model)</div>
        </div>
        <div class="summary-card">
            <div class="tool-label chai-label">Chai-1</div>
            <div class="value chai-label">{c_best.get('plddt', 0):.3f}</div>
            <div class="metric-name">pLDDT (best model)</div>
        </div>
        <div class="summary-card">
            <div class="tool-label boltz-label">Boltz-2</div>
            <div class="value boltz-label">{b_best.get('ptm', 0):.3f}</div>
            <div class="metric-name">pTM (best model)</div>
        </div>
        <div class="summary-card">
            <div class="tool-label chai-label">Chai-1</div>
            <div class="value chai-label">{c_best.get('ptm', 0):.3f}</div>
            <div class="metric-name">pTM (best model)</div>
        </div>
        <div class="summary-card">
            <div class="tool-label" style="color:#475569;">PAE Winner</div>
            <div class="value" style="font-size:1.1rem;color:#0f766e;">
                {winner(b_best.get('pae_mean'), c_best.get('pae_mean'), lower_is_better=True)}
            </div>
            <div class="metric-name">Lower PAE = better</div>
        </div>
    </div>

    <!-- Side-by-side metric bar chart -->
    <div class="card">
        <div class="card-header"><i class="fas fa-chart-bar"></i> Metric Comparison (Best Models)</div>
        <div class="plot-container" id="metricChart"></div>
    </div>

    <!-- PAE mean comparison -->
    <div class="card">
        <div class="card-header"><i class="fas fa-th"></i> PAE &amp; PDE Comparison (Å)</div>
        <div class="plot-container" id="errorChart"></div>
    </div>

    <!-- Detailed table -->
    <div class="card">
        <div class="card-header"><i class="fas fa-table"></i> All Models — Detailed Metrics</div>
        <div class="card-body">
            <div class="table-responsive">
                <table class="table table-striped table-hover">
                    <thead>
                        <tr>
                            <th>Tool</th><th>Sample</th>
                            <th>pLDDT (0–1)</th><th>pTM (0–1)</th><th>ipTM (0–1)</th>
                            <th>PAE mean (Å)</th><th>PDE mean (Å)</th>
                        </tr>
                    </thead>
                    <tbody>{table_rows}</tbody>
                </table>
            </div>
        </div>
    </div>

</div>

<script>
// ── Metric comparison chart ──
Plotly.newPlot('metricChart', [
    {{
        name: 'Boltz-2',
        type: 'bar',
        x: {json.dumps(labels)},
        y: {json.dumps(boltz_vals)},
        marker: {{ color: '#0284c7' }},
        text: {json.dumps([str(v) for v in boltz_vals])},
        textposition: 'outside'
    }},
    {{
        name: 'Chai-1',
        type: 'bar',
        x: {json.dumps(labels)},
        y: {json.dumps(chai_vals)},
        marker: {{ color: '#7c3aed' }},
        text: {json.dumps([str(v) for v in chai_vals])},
        textposition: 'outside'
    }}
], {{
    barmode: 'group',
    yaxis: {{ title: 'Score (0–1)', range: [0, 1.05] }},
    xaxis: {{ title: 'Metric' }},
    template: 'plotly_white',
    legend: {{ orientation: 'h', y: -0.2 }},
    height: 380,
    margin: {{ t: 30, b: 80 }}
}});

// ── PAE / PDE comparison chart ──
Plotly.newPlot('errorChart', [
    {{
        name: 'Boltz-2',
        type: 'bar',
        x: ['PAE mean', 'PDE mean'],
        y: [{round(b_best.get('pae_mean') or 0, 4)}, {round(b_best.get('pde_mean') or 0, 4)}],
        marker: {{ color: '#0284c7' }}
    }},
    {{
        name: 'Chai-1',
        type: 'bar',
        x: ['PAE mean', 'PDE mean'],
        y: [{round(c_best.get('pae_mean') or 0, 4)}, {round(c_best.get('pde_mean') or 0, 4)}],
        marker: {{ color: '#7c3aed' }}
    }}
], {{
    barmode: 'group',
    yaxis: {{ title: 'Error (Å)' }},
    template: 'plotly_white',
    legend: {{ orientation: 'h', y: -0.2 }},
    height: 350,
    margin: {{ t: 30, b: 80 }},
    annotations: [{{
        text: 'Lower is better',
        showarrow: false, x: 0.5, y: 1.08,
        xref: 'paper', yref: 'paper',
        font: {{ color: '#64748b', size: 12 }}
    }}]
}});
</script>
</body>
</html>"""


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Generate Boltz-2 vs Chai-1 comparison report')
    parser.add_argument('--boltz-summary', default='boltz_summary.json')
    parser.add_argument('--chai-summary',  default='chai_summary.json')
    parser.add_argument('--output',        default='comparison_report.html')
    args = parser.parse_args()

    print("\nLoading summaries:")
    boltz_data = load_summary(args.boltz_summary)
    chai_data  = load_summary(args.chai_summary)

    print("\nGenerating report...")
    html = generate_html(boltz_data, chai_data)

    with open(args.output, 'w') as f:
        f.write(html)

    print(f"  Report written to: {args.output}")
    print("\nNode 05 completed ✓")


if __name__ == '__main__':
    main()
