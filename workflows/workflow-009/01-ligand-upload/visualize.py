#!/usr/bin/env python3
"""
Visualize ligand SMILES file - Node 1

Usage:
    python visualize.py [input.smi] [output.html]

Default:
    python visualize.py ligand.smi ligand_viz.html
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


def mol_to_svg_with_h(smiles, width=800, height=600):
    """Convert SMILES to SVG with explicit H atoms and atom map numbers."""
    # Load SMILES without sanitization to preserve explicit H atoms with map numbers
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return "<p>Invalid molecule</p>"
    # Sanitize without adjusting/removing Hs
    Chem.SanitizeMol(
        mol,
        sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL
        ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS,
    )

    AllChem.Compute2DCoords(mol)

    # Highlight all heavy atoms (non-H) in light blue
    highlight_atoms = []
    highlight_colors = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if atom.GetAtomicNum() != 1:  # Non-hydrogen
            highlight_atoms.append(idx)
            highlight_colors[idx] = (0.7, 0.85, 1.0)  # Light blue

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    opts = drawer.drawOptions()
    # Set custom labels for ALL atoms to show their map numbers (forces H atoms to be visible)
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        map_num = atom.GetAtomMapNum()
        symbol = atom.GetSymbol()
        if map_num > 0:
            opts.atomLabels[idx] = f"{symbol}:{map_num}"
        else:
            opts.atomLabels[idx] = symbol

    drawer.DrawMolecule(
        mol, highlightAtoms=highlight_atoms, highlightAtomColors=highlight_colors
    )
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def mol_to_molblock(mol):
    """Convert molecule to MOL block for 3Dmol.js."""
    if mol is None:
        return ""
    if mol.GetNumConformers() == 0:
        # Only add Hs if not already present
        mol_no_h = Chem.RemoveHs(mol)
        if mol.GetNumAtoms() == mol_no_h.GetNumAtoms():
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
        .mol-3d {{ width: 100%; height: 400px; border-radius: 8px; position: relative; }}
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
            <h3 style="color: #ea7d3d; margin-bottom: 16px;">2D Structure</h3>
            <div class="mol-2d">{svg_2d}</div>
            <div class="mol-3d" id="viewer"></div>
        </div>
        <div class="card">
            <h3 style="color: #ea7d3d; margin-bottom: 16px;">2D Structure (with H Atoms and Atom IDs)</h3>
            <p style="color: #888; margin-bottom: 16px;">Use these atom IDs when selecting attachment points in Node 2</p>
            <div class="mol-2d">{svg_2d_with_h}</div>
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
    input_file = sys.argv[1] if len(sys.argv) > 1 else "ligand.smi"
    output_file = sys.argv[2] if len(sys.argv) > 2 else "ligand_viz.html"

    if not Path(input_file).exists():
        print(f"Error: {input_file} not found")
        sys.exit(1)

    with open(input_file, "r") as f:
        smiles = f.read().strip()

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Error: No valid molecule found")
        sys.exit(1)

    svg_2d = mol_to_svg(mol)
    svg_2d_with_h = mol_to_svg_with_h(smiles)
    molblock = mol_to_molblock(mol).replace("`", "\\`").replace("\n", "\\n")

    html = HTML_TEMPLATE.format(
        filename=Path(input_file).name,
        svg_2d=svg_2d,
        svg_2d_with_h=svg_2d_with_h,
        num_atoms=mol.GetNumAtoms(),
        mw=round(Descriptors.MolWt(mol), 2),
        logp=round(Descriptors.MolLogP(mol), 2),
        smiles=smiles,
        molblock=molblock,
    )

    with open(output_file, "w") as f:
        f.write(html)

    print(f"Generated: {output_file}")


if __name__ == "__main__":
    main()
