#!/usr/bin/env python3
"""
Visualize ligand SDF file - Node 1

Usage:
    python visualize.py [input.sdf] [output.html]

Default:
    python visualize.py ligand.sdf ligand_viz.html
"""

import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D


def mol_to_svg(mol, width=400, height=300):
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


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Ligand Visualization</title>
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: #1a1a2e; color: #eee; padding: 20px;
        }}
        .header {{ text-align: center; padding: 20px; background: #16213e; border-radius: 12px; margin-bottom: 20px; }}
        .header h1 {{ color: #ea7d3d; }}
        .container {{ max-width: 900px; margin: 0 auto; }}
        .card {{ background: #16213e; border-radius: 12px; padding: 20px; margin-bottom: 20px; }}
        .mol-2d {{ text-align: center; background: white; border-radius: 8px; padding: 10px; margin-bottom: 20px; }}
        .mol-3d {{ width: 100%; height: 400px; border-radius: 8px; }}
        .props {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 16px; }}
        .prop {{ background: #1a1a2e; padding: 16px; border-radius: 8px; text-align: center; }}
        .prop .value {{ font-size: 24px; color: #ea7d3d; font-weight: bold; }}
        .prop .label {{ color: #888; font-size: 14px; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>Ligand Visualization</h1>
            <p>{filename}</p>
        </div>
        <div class="card">
            <div class="mol-2d">{svg_2d}</div>
            <div class="mol-3d" id="viewer"></div>
        </div>
        <div class="props">
            <div class="prop"><div class="value">{num_atoms}</div><div class="label">Atoms</div></div>
            <div class="prop"><div class="value">{mw}</div><div class="label">Mol. Weight</div></div>
            <div class="prop"><div class="value">{logp}</div><div class="label">LogP</div></div>
        </div>
        <div class="card">
            <p><strong>SMILES:</strong> {smiles}</p>
        </div>
    </div>
    <script>
        const viewer = $3Dmol.createViewer('viewer', {{backgroundColor: '#16213e'}});
        viewer.addModel(`{molblock}`, 'sdf');
        viewer.setStyle({{}}, {{stick: {{colorscheme: 'orangeCarbon', radius: 0.15}}}});
        viewer.zoomTo();
        viewer.render();
    </script>
</body>
</html>
"""


def main():
    input_file = sys.argv[1] if len(sys.argv) > 1 else "ligand.sdf"
    output_file = sys.argv[2] if len(sys.argv) > 2 else "ligand_viz.html"

    if not Path(input_file).exists():
        print(f"Error: {input_file} not found")
        sys.exit(1)

    supplier = Chem.SDMolSupplier(input_file)
    mol = next((m for m in supplier if m is not None), None)

    if mol is None:
        print("Error: No valid molecule found")
        sys.exit(1)

    svg_2d = mol_to_svg(mol)
    molblock = mol_to_molblock(mol).replace('`', '\\`').replace('\n', '\\n')

    html = HTML_TEMPLATE.format(
        filename=Path(input_file).name,
        svg_2d=svg_2d,
        num_atoms=mol.GetNumAtoms(),
        mw=round(Descriptors.MolWt(mol), 2),
        logp=round(Descriptors.MolLogP(mol), 2),
        smiles=Chem.MolToSmiles(mol),
        molblock=molblock
    )

    with open(output_file, 'w') as f:
        f.write(html)

    print(f"Generated: {output_file}")


if __name__ == '__main__':
    main()
