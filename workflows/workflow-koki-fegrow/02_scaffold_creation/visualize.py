#!/usr/bin/env python3
"""
Visualize scaffold from pickle file - Node 2

Usage:
    python visualize.py [scaffold.pkl] [output.html]

Default:
    python visualize.py scaffold.pkl scaffold_viz.html
"""

import sys
import pickle
from pathlib import Path

try:
    import cloudpickle
except ImportError:
    cloudpickle = None

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D


def mol_to_svg(mol, width=500, height=400):
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


def load_scaffold(pkl_path):
    """Load scaffold from pickle file and extract RDKit molecule."""
    try:
        with open(pkl_path, 'rb') as f:
            scaffold = pickle.load(f)
    except:
        if cloudpickle:
            with open(pkl_path, 'rb') as f:
                scaffold = cloudpickle.load(f)
        else:
            raise

    # Extract RDKit molecule from FEGrow RMol
    if hasattr(scaffold, 'mol'):
        return scaffold.mol
    elif hasattr(scaffold, 'GetMol'):
        return scaffold.GetMol()
    elif isinstance(scaffold, Chem.Mol):
        return scaffold
    else:
        raise ValueError(f"Cannot extract molecule from: {type(scaffold)}")


def find_attachment_point(mol):
    """Find attachment point (atom with atomic number 0)."""
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            return atom.GetIdx()
    return None


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Scaffold Visualization</title>
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: #1a1a2e; color: #eee; padding: 20px;
        }}
        .header {{ text-align: center; padding: 20px; background: #16213e; border-radius: 12px; margin-bottom: 20px; }}
        .header h1 {{ color: #ea7d3d; }}
        .container {{ max-width: 1000px; margin: 0 auto; }}
        .grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
        .card {{ background: #16213e; border-radius: 12px; padding: 20px; }}
        .mol-2d {{ text-align: center; background: white; border-radius: 8px; padding: 10px; }}
        .mol-3d {{ width: 100%; height: 400px; border-radius: 8px; }}
        .props {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 16px; margin-bottom: 20px; }}
        .prop {{ background: #16213e; padding: 16px; border-radius: 8px; text-align: center; border: 1px solid #2a2a4a; }}
        .prop .value {{ font-size: 24px; color: #ea7d3d; font-weight: bold; }}
        .prop .label {{ color: #888; font-size: 14px; }}
        .info {{ background: #16213e; padding: 20px; border-radius: 12px; margin-top: 20px; }}
        .info h3 {{ color: #ea7d3d; margin-bottom: 12px; }}
        .highlight {{ color: #22c55e; font-weight: bold; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>Scaffold Visualization</h1>
            <p>Node 2: Scaffold Creation</p>
        </div>
        <div class="props">
            <div class="prop"><div class="value">{num_atoms}</div><div class="label">Atoms</div></div>
            <div class="prop"><div class="value">{num_bonds}</div><div class="label">Bonds</div></div>
            <div class="prop"><div class="value">{attachment}</div><div class="label">Attachment Point</div></div>
            <div class="prop"><div class="value">{mw}</div><div class="label">Mol. Weight</div></div>
        </div>
        <div class="grid">
            <div class="card">
                <h3 style="color: #ea7d3d; margin-bottom: 16px;">2D Structure</h3>
                <div class="mol-2d">{svg_2d}</div>
            </div>
            <div class="card">
                <h3 style="color: #ea7d3d; margin-bottom: 16px;">3D Structure</h3>
                <div class="mol-3d" id="viewer"></div>
            </div>
        </div>
        <div class="info">
            <h3>Scaffold Information</h3>
            <p>The scaffold has an attachment point at atom index <span class="highlight">{attachment}</span>.</p>
            <p>This is where R-groups and linkers will be attached during chemical space generation.</p>
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
    input_file = sys.argv[1] if len(sys.argv) > 1 else "scaffold.pkl"
    output_file = sys.argv[2] if len(sys.argv) > 2 else "scaffold_viz.html"

    if not Path(input_file).exists():
        print(f"Error: {input_file} not found")
        sys.exit(1)

    mol = load_scaffold(input_file)
    attachment = find_attachment_point(mol)

    svg_2d = mol_to_svg(mol)
    molblock = mol_to_molblock(mol).replace('`', '\\`').replace('\n', '\\n')

    html = HTML_TEMPLATE.format(
        svg_2d=svg_2d,
        num_atoms=mol.GetNumAtoms(),
        num_bonds=mol.GetNumBonds(),
        attachment=attachment if attachment is not None else "N/A",
        mw=round(Descriptors.MolWt(mol), 2),
        smiles=Chem.MolToSmiles(mol),
        molblock=molblock
    )

    with open(output_file, 'w') as f:
        f.write(html)

    print(f"Generated: {output_file}")


if __name__ == '__main__':
    main()
