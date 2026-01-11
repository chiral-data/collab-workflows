#!/usr/bin/env python3
"""
Visualize final report with top compounds - Node 7

Usage:
    python visualize.py [top_compounds.sdf] [top_compounds_report.csv] [summary.txt] [output.html]

Default:
    python visualize.py top_compounds.sdf top_compounds_report.csv summary.txt report_viz.html
"""

import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D

try:
    import pandas as pd
except ImportError:
    pd = None


def mol_to_svg(mol, width=300, height=220):
    """Convert RDKit molecule to SVG string."""
    if mol is None:
        return "<p>Invalid molecule</p>"
    if mol.GetNumConformers() == 0:
        AllChem.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def mol_to_molblock(mol):
    """Convert molecule to MOL block for 3Dmol.js."""
    if mol is None:
        return ""
    if mol.GetNumConformers() == 0:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(mol)
    return Chem.MolToMolBlock(mol)


def get_score(mol):
    """Get binding score from molecule properties."""
    for prop in ['score', 'minimizedAffinity', 'Binding_Energy']:
        if mol.HasProp(prop):
            try:
                return float(mol.GetProp(prop))
            except:
                pass
    return None


def get_status(score):
    """Classify binding energy."""
    if score is None:
        return "", ""
    if score < -8:
        return "EXCELLENT", "excellent"
    elif score < -5:
        return "GOOD", "good"
    elif score < 0:
        return "FAIR", "fair"
    else:
        return "POOR", "poor"


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>FEGrow Report - Top Compounds</title>
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: #1a1a2e; color: #eee; padding: 20px;
        }}
        .header {{ text-align: center; padding: 24px; background: linear-gradient(135deg, #16213e, #1a1a2e); border-radius: 12px; margin-bottom: 20px; border: 1px solid #ea7d3d; }}
        .header h1 {{ color: #ea7d3d; font-size: 32px; }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .stats {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 16px; margin-bottom: 20px; }}
        .stat {{ background: #16213e; padding: 24px; border-radius: 12px; text-align: center; border: 1px solid #2a2a4a; }}
        .stat .value {{ font-size: 32px; color: #ea7d3d; font-weight: bold; }}
        .stat .label {{ color: #888; font-size: 14px; margin-top: 8px; }}
        .tabs {{ display: flex; gap: 8px; margin-bottom: 20px; }}
        .tab {{ padding: 12px 24px; background: #16213e; border: 1px solid #2a2a4a; border-radius: 8px; cursor: pointer; color: #888; }}
        .tab:hover, .tab.active {{ background: #ea7d3d; color: white; border-color: #ea7d3d; }}
        .tab-content {{ display: none; }}
        .tab-content.active {{ display: block; }}
        .grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(340px, 1fr)); gap: 16px; }}
        .card {{ background: #16213e; border-radius: 12px; overflow: hidden; border: 1px solid #2a2a4a; transition: all 0.2s; }}
        .card:hover {{ transform: translateY(-4px); box-shadow: 0 8px 32px rgba(234,125,61,0.2); }}
        .card-header {{ background: #1a1a2e; padding: 12px 16px; display: flex; justify-content: space-between; align-items: center; }}
        .card-header h3 {{ color: #ea7d3d; font-size: 16px; }}
        .rank {{ background: #ea7d3d; color: white; width: 32px; height: 32px; border-radius: 50%; display: flex; align-items: center; justify-content: center; font-weight: bold; }}
        .badge {{ padding: 4px 12px; border-radius: 20px; font-size: 12px; font-weight: 600; }}
        .badge.excellent {{ background: #22c55e; color: white; }}
        .badge.good {{ background: #3b82f6; color: white; }}
        .badge.fair {{ background: #eab308; color: black; }}
        .badge.poor {{ background: #ef4444; color: white; }}
        .card-body {{ padding: 16px; }}
        .mol-2d {{ text-align: center; background: white; border-radius: 8px; padding: 8px; margin-bottom: 14px; }}
        .mol-3d {{ width: 100%; height: 240px; border-radius: 8px; margin-bottom: 14px; }}
        .props {{ font-size: 13px; }}
        .props div {{ display: flex; justify-content: space-between; padding: 6px 0; border-bottom: 1px solid #2a2a4a; }}
        .props div:last-child {{ border-bottom: none; }}
        .props .label {{ color: #888; }}
        .props .value {{ color: #fff; font-weight: 500; }}
        table.data-table {{ width: 100%; border-collapse: collapse; background: #16213e; border-radius: 12px; overflow: hidden; }}
        table.data-table th {{ background: #ea7d3d; color: white; padding: 14px 16px; text-align: left; font-weight: 600; }}
        table.data-table td {{ padding: 12px 16px; border-bottom: 1px solid #2a2a4a; }}
        table.data-table tr:hover {{ background: #1a1a2e; }}
        .status-excellent {{ color: #22c55e; font-weight: 600; }}
        .status-good {{ color: #3b82f6; font-weight: 600; }}
        .status-fair {{ color: #eab308; font-weight: 600; }}
        .status-poor {{ color: #ef4444; font-weight: 600; }}
        .summary-box {{ background: #16213e; padding: 24px; border-radius: 12px; border: 1px solid #2a2a4a; }}
        .summary-box h2 {{ color: #ea7d3d; margin-bottom: 16px; }}
        .summary-box pre {{ background: #1a1a2e; padding: 20px; border-radius: 8px; overflow-x: auto; font-family: 'Monaco', 'Consolas', monospace; font-size: 13px; line-height: 1.7; white-space: pre-wrap; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>FEGrow Results Report</h1>
            <p>Top Compounds from Active Learning Docking</p>
        </div>
        <div class="stats">
            <div class="stat"><div class="value">{num_compounds}</div><div class="label">Top Compounds</div></div>
            <div class="stat"><div class="value">{best_energy}</div><div class="label">Best Energy (kcal/mol)</div></div>
            <div class="stat"><div class="value">{avg_energy}</div><div class="label">Avg Energy (kcal/mol)</div></div>
            <div class="stat"><div class="value">{worst_energy}</div><div class="label">Worst Energy (kcal/mol)</div></div>
        </div>
        <div class="tabs">
            <div class="tab active" data-tab="molecules">Molecules</div>
            <div class="tab" data-tab="table">Data Table</div>
            <div class="tab" data-tab="summary">Summary</div>
        </div>
        <div class="tab-content active" data-tab="molecules">
            <div class="grid">{cards}</div>
        </div>
        <div class="tab-content" data-tab="table">{table}</div>
        <div class="tab-content" data-tab="summary">
            <div class="summary-box">
                <h2>Analysis Summary</h2>
                <pre>{summary}</pre>
            </div>
        </div>
    </div>
    <script>
        document.querySelectorAll('.mol-3d').forEach(function(container) {{
            const molData = container.dataset.mol;
            if (molData) {{
                const viewer = $3Dmol.createViewer(container, {{backgroundColor: '#16213e'}});
                viewer.addModel(molData.replace(/&#10;/g, '\\n').replace(/&quot;/g, '"'), 'sdf');
                viewer.setStyle({{}}, {{stick: {{colorscheme: 'orangeCarbon', radius: 0.12}}}});
                viewer.zoomTo();
                viewer.render();
            }}
        }});
        document.querySelectorAll('.tab').forEach(function(tab) {{
            tab.addEventListener('click', function() {{
                const tabId = this.dataset.tab;
                document.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
                document.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
                this.classList.add('active');
                document.querySelector('.tab-content[data-tab="' + tabId + '"]').classList.add('active');
            }});
        }});
    </script>
</body>
</html>
"""

CARD_TEMPLATE = """
<div class="card">
    <div class="card-header">
        <div style="display:flex;align-items:center;gap:12px">
            <div class="rank">{rank}</div>
            <h3>Compound #{rank}</h3>
        </div>
        {badge}
    </div>
    <div class="card-body">
        <div class="mol-2d">{svg_2d}</div>
        <div class="mol-3d" data-mol="{molblock}"></div>
        <div class="props">
            <div><span class="label">Binding Energy</span><span class="value">{score} kcal/mol</span></div>
            <div><span class="label">Molecular Weight</span><span class="value">{mw}</span></div>
            <div><span class="label">LogP</span><span class="value">{logp}</span></div>
        </div>
    </div>
</div>
"""


def main():
    sdf_file = sys.argv[1] if len(sys.argv) > 1 else "top_compounds.sdf"
    csv_file = sys.argv[2] if len(sys.argv) > 2 else "top_compounds_report.csv"
    txt_file = sys.argv[3] if len(sys.argv) > 3 else "summary.txt"
    output_file = sys.argv[4] if len(sys.argv) > 4 else "report_viz.html"

    # Load SDF
    molecules = []
    if Path(sdf_file).exists():
        supplier = Chem.SDMolSupplier(sdf_file)
        molecules = [mol for mol in supplier if mol is not None]
        print(f"Loaded {len(molecules)} molecules from {sdf_file}")
    else:
        print(f"Warning: {sdf_file} not found")

    # Load CSV
    csv_data = None
    if Path(csv_file).exists() and pd is not None:
        csv_data = pd.read_csv(csv_file)
        print(f"Loaded {len(csv_data)} rows from {csv_file}")
    elif Path(csv_file).exists():
        print("Warning: pandas not available, CSV skipped")

    # Load summary
    summary_text = "No summary available"
    if Path(txt_file).exists():
        with open(txt_file, 'r') as f:
            summary_text = f.read()
        print(f"Loaded summary from {txt_file}")

    # Calculate statistics
    energies = [get_score(mol) for mol in molecules if get_score(mol) is not None]
    best_energy = f"{min(energies):.2f}" if energies else "N/A"
    worst_energy = f"{max(energies):.2f}" if energies else "N/A"
    avg_energy = f"{sum(energies)/len(energies):.2f}" if energies else "N/A"

    # Generate molecule cards
    cards = []
    for i, mol in enumerate(molecules):
        score = get_score(mol)
        status_text, status_class = get_status(score)

        svg_2d = mol_to_svg(mol)
        molblock = mol_to_molblock(mol).replace('"', '&quot;').replace('\n', '&#10;')

        badge = f'<span class="badge {status_class}">{status_text}</span>' if status_text else ""

        card = CARD_TEMPLATE.format(
            rank=i + 1,
            badge=badge,
            svg_2d=svg_2d,
            molblock=molblock,
            score=f"{score:.2f}" if score else "N/A",
            mw=round(Descriptors.MolWt(mol), 1),
            logp=round(Descriptors.MolLogP(mol), 2)
        )
        cards.append(card)

    # Generate data table
    table_html = "<p>No CSV data available</p>"
    if csv_data is not None and not csv_data.empty:
        table_html = '<table class="data-table"><thead><tr>'
        for col in csv_data.columns:
            table_html += f'<th>{col}</th>'
        table_html += '</tr></thead><tbody>'

        for _, row in csv_data.iterrows():
            table_html += '<tr>'
            for col in csv_data.columns:
                val = row[col]
                if col == 'Status' and val in ['EXCELLENT', 'GOOD', 'FAIR', 'POOR']:
                    table_html += f'<td class="status-{str(val).lower()}">{val}</td>'
                else:
                    table_html += f'<td>{val}</td>'
            table_html += '</tr>'
        table_html += '</tbody></table>'

    html = HTML_TEMPLATE.format(
        num_compounds=len(molecules),
        best_energy=best_energy,
        avg_energy=avg_energy,
        worst_energy=worst_energy,
        cards="\n".join(cards),
        table=table_html,
        summary=summary_text
    )

    with open(output_file, 'w') as f:
        f.write(html)

    print(f"Generated: {output_file}")


if __name__ == '__main__':
    main()
