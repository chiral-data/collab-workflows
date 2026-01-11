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


def mol_to_svg(mol, width=800, height=600):
    """Convert RDKit molecule to SVG with atom indices and explicit hydrogens."""
    if mol is None or mol.GetNumAtoms() == 0:
        return "<p>Invalid molecule</p>"

    # Add explicit hydrogens
    mol_with_h = Chem.AddHs(mol)

    # Generate 2D coordinates
    AllChem.Compute2DCoords(mol_with_h)

    # Highlight all scaffold atoms (non-H) and attachment points specially
    highlight_atoms = []
    highlight_colors = {}
    for atom in mol_with_h.GetAtoms():
        idx = atom.GetIdx()
        if atom.GetAtomicNum() == 0:  # Dummy atom (attachment point) - orange
            highlight_atoms.append(idx)
            highlight_colors[idx] = (0.92, 0.49, 0.24)  # Orange (#ea7d3d)
        elif atom.GetAtomicNum() != 1:  # Non-hydrogen (scaffold) - light blue
            highlight_atoms.append(idx)
            highlight_colors[idx] = (0.7, 0.85, 1.0)  # Light blue

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    opts.addAtomIndices = True  # Show atom index numbers

    drawer.DrawMolecule(mol_with_h, highlightAtoms=highlight_atoms, highlightAtomColors=highlight_colors)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


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


def find_attachment_points(mol):
    """Find all attachment points (atoms with atomic number 0)."""
    points = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 0:
            points.append(atom.GetIdx())
    return points


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Scaffold Visualization</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: #1a1a2e; color: #eee; padding: 20px;
        }}
        .header {{ text-align: center; padding: 20px; background: #16213e; border-radius: 12px; margin-bottom: 20px; }}
        .header h1 {{ color: #ea7d3d; }}
        .container {{ max-width: 1000px; margin: 0 auto; }}
        .card {{ background: #16213e; border-radius: 12px; padding: 20px; margin-bottom: 20px; }}
        .mol-2d {{ text-align: center; background: white; border-radius: 8px; padding: 10px; }}
        .props {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 16px; margin-bottom: 20px; }}
        .prop {{ background: #16213e; padding: 16px; border-radius: 8px; text-align: center; border: 1px solid #2a2a4a; }}
        .prop .value {{ font-size: 24px; color: #ea7d3d; font-weight: bold; }}
        .prop .label {{ color: #888; font-size: 14px; }}
        .info {{ background: #16213e; padding: 20px; border-radius: 12px; }}
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
            <div class="prop"><div class="value">{num_atoms}</div><div class="label">Heavy Atoms</div></div>
            <div class="prop"><div class="value">{num_atoms_with_h}</div><div class="label">Total Atoms (with H)</div></div>
            <div class="prop"><div class="value">{num_bonds}</div><div class="label">Bonds</div></div>
            <div class="prop"><div class="value">{mw}</div><div class="label">Mol. Weight</div></div>
        </div>
        <div class="card">
            <h3 style="color: #ea7d3d; margin-bottom: 16px;">2D Structure (with Atom Indices)</h3>
            <div class="mol-2d">{svg_2d}</div>
        </div>
        <div class="info">
            <h3>Scaffold Information</h3>
            <p>Attachment point(s) at atom index: <span class="highlight">{attachments}</span></p>
            <p>Attachment points are highlighted in orange. These are where R-groups will be attached.</p>
            <p><strong>SMILES:</strong> {smiles}</p>
        </div>
    </div>
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
    attachments = find_attachment_points(mol)
    mol_with_h = Chem.AddHs(mol)

    svg_2d = mol_to_svg(mol)

    html = HTML_TEMPLATE.format(
        svg_2d=svg_2d,
        num_atoms=mol.GetNumAtoms(),
        num_atoms_with_h=mol_with_h.GetNumAtoms(),
        num_bonds=mol.GetNumBonds(),
        attachments=', '.join(str(a) for a in attachments) if attachments else "N/A",
        mw=round(Descriptors.MolWt(mol), 2),
        smiles=Chem.MolToSmiles(mol)
    )

    with open(output_file, 'w') as f:
        f.write(html)

    print(f"Generated: {output_file}")


if __name__ == '__main__':
    main()
