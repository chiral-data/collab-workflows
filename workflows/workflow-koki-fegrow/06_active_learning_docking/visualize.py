#!/usr/bin/env python3
"""
Visualize evaluated molecules from active learning docking - Node 6

Usage:
    python visualize.py [chemspace_evaluated.sdf] [output.html]

Default:
    python visualize.py chemspace_evaluated.sdf evaluated_viz.html
"""

import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D


def mol_to_svg(mol, width=280, height=200):
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
    <title>Evaluated Molecules</title>
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: #1a1a2e; color: #eee; padding: 20px;
        }}
        .header {{ text-align: center; padding: 20px; background: #16213e; border-radius: 12px; margin-bottom: 20px; }}
        .header h1 {{ color: #ea7d3d; }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .stats {{ display: grid; grid-template-columns: repeat(5, 1fr); gap: 16px; margin-bottom: 20px; }}
        .stat {{ background: #16213e; padding: 20px; border-radius: 12px; text-align: center; border: 1px solid #2a2a4a; }}
        .stat .value {{ font-size: 28px; color: #ea7d3d; font-weight: bold; }}
        .stat .label {{ color: #888; font-size: 13px; margin-top: 4px; }}
        .grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(320px, 1fr)); gap: 16px; }}
        .card {{ background: #16213e; border-radius: 12px; overflow: hidden; border: 1px solid #2a2a4a; }}
        .card:hover {{ transform: translateY(-2px); box-shadow: 0 8px 24px rgba(234,125,61,0.2); }}
        .card-header {{ background: #1a1a2e; padding: 10px 14px; display: flex; justify-content: space-between; align-items: center; }}
        .card-header h3 {{ color: #ea7d3d; font-size: 14px; }}
        .badge {{ padding: 4px 10px; border-radius: 20px; font-size: 11px; font-weight: 600; }}
        .badge.excellent {{ background: #22c55e; color: white; }}
        .badge.good {{ background: #3b82f6; color: white; }}
        .badge.fair {{ background: #eab308; color: black; }}
        .badge.poor {{ background: #ef4444; color: white; }}
        .card-body {{ padding: 14px; }}
        .mol-2d {{ text-align: center; background: white; border-radius: 8px; padding: 8px; margin-bottom: 12px; }}
        .mol-3d {{ width: 100%; height: 220px; border-radius: 8px; margin-bottom: 12px; }}
        .props {{ font-size: 12px; }}
        .props div {{ display: flex; justify-content: space-between; padding: 4px 0; border-bottom: 1px solid #2a2a4a; }}
        .props div:last-child {{ border-bottom: none; }}
        .props .label {{ color: #888; }}
        .props .value {{ color: #fff; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>Evaluated Molecules</h1>
            <p>Node 6: Active Learning Docking Results</p>
        </div>
        <div class="stats">
            <div class="stat"><div class="value">{total}</div><div class="label">Total Molecules</div></div>
            <div class="stat"><div class="value" style="color:#22c55e">{excellent}</div><div class="label">Excellent (&lt;-8)</div></div>
            <div class="stat"><div class="value" style="color:#3b82f6">{good}</div><div class="label">Good (-5 to -8)</div></div>
            <div class="stat"><div class="value" style="color:#eab308">{fair}</div><div class="label">Fair (0 to -5)</div></div>
            <div class="stat"><div class="value" style="color:#ef4444">{poor}</div><div class="label">Poor (&gt;0)</div></div>
        </div>
        <div class="grid">
            {cards}
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
    </script>
</body>
</html>
"""

CARD_TEMPLATE = """
<div class="card">
    <div class="card-header">
        <h3>Molecule {index}</h3>
        {badge}
    </div>
    <div class="card-body">
        <div class="mol-2d">{svg_2d}</div>
        <div class="mol-3d" data-mol="{molblock}"></div>
        <div class="props">
            <div><span class="label">Score</span><span class="value">{score}</span></div>
            <div><span class="label">MW</span><span class="value">{mw}</span></div>
            <div><span class="label">LogP</span><span class="value">{logp}</span></div>
        </div>
    </div>
</div>
"""


def main():
    input_file = sys.argv[1] if len(sys.argv) > 1 else "chemspace_evaluated.sdf"
    output_file = sys.argv[2] if len(sys.argv) > 2 else "evaluated_viz.html"

    if not Path(input_file).exists():
        print(f"Error: {input_file} not found")
        sys.exit(1)

    supplier = Chem.SDMolSupplier(input_file)
    molecules = [(i, mol) for i, mol in enumerate(supplier) if mol is not None]

    if not molecules:
        print("Error: No valid molecules found")
        sys.exit(1)

    print(f"Loaded {len(molecules)} molecules from {input_file}")

    # Sort by score (best first)
    molecules.sort(key=lambda x: get_score(x[1]) or 0)

    # Count by status
    counts = {"excellent": 0, "good": 0, "fair": 0, "poor": 0}
    for _, mol in molecules:
        score = get_score(mol)
        _, status_class = get_status(score)
        if status_class:
            counts[status_class] += 1

    # Generate cards
    cards = []
    for idx, mol in molecules:
        score = get_score(mol)
        status_text, status_class = get_status(score)

        svg_2d = mol_to_svg(mol)
        molblock = mol_to_molblock(mol).replace('"', '&quot;').replace('\n', '&#10;')

        badge = f'<span class="badge {status_class}">{status_text}</span>' if status_text else ""

        card = CARD_TEMPLATE.format(
            index=idx + 1,
            badge=badge,
            svg_2d=svg_2d,
            molblock=molblock,
            score=f"{score:.2f}" if score else "N/A",
            mw=round(Descriptors.MolWt(mol), 1),
            logp=round(Descriptors.MolLogP(mol), 2)
        )
        cards.append(card)

    html = HTML_TEMPLATE.format(
        total=len(molecules),
        excellent=counts["excellent"],
        good=counts["good"],
        fair=counts["fair"],
        poor=counts["poor"],
        cards="\n".join(cards)
    )

    with open(output_file, 'w') as f:
        f.write(html)

    print(f"Generated: {output_file}")


if __name__ == '__main__':
    main()
