#!/usr/bin/env python3
"""
Visualize chemical space summary with 2D molecule grid - Node 5

Features:
- Grid view of all generated molecules
- Click to enlarge molecule view
- Highlights R-groups/linkers (parts not in scaffold) in orange

Usage:
    python visualize.py [chemspace.pkl] [output.html]

Default:
    python visualize.py chemspace.pkl chemspace_viz.html
"""

import pickle
import sys
from pathlib import Path

try:
    import cloudpickle
except ImportError:
    cloudpickle = None

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D


def load_chemspace(pkl_path):
    """Load chemical space from pickle file."""
    try:
        if cloudpickle:
            with open(pkl_path, "rb") as f:
                return cloudpickle.load(f)
        else:
            with open(pkl_path, "rb") as f:
                return pickle.load(f)
    except:
        with open(pkl_path, "rb") as f:
            return pickle.load(f)


def remove_dummy_atoms(mol):
    """Remove dummy atoms (atomic number 0) from molecule."""
    if mol is None:
        return None
    dummy_indices = [
        atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0
    ]
    if not dummy_indices:
        return mol
    edit_mol = Chem.RWMol(mol)
    for idx in sorted(dummy_indices, reverse=True):
        edit_mol.RemoveAtom(idx)
    return edit_mol.GetMol()


def get_scaffold_match(mol, scaffold_mol):
    """Get atom indices that match the scaffold (to identify R-groups).

    Uses substructure matching to find scaffold atoms in the molecule.
    Atoms NOT in the match are the R-groups/linkers added at attachment point.
    """
    if mol is None or scaffold_mol is None:
        return []

    # Clean both molecules (remove dummy atoms and Hs for matching)
    clean_mol = remove_dummy_atoms(mol)
    clean_scaffold = remove_dummy_atoms(scaffold_mol)

    # Remove Hs for cleaner matching
    if clean_mol is not None:
        clean_mol = Chem.RemoveHs(clean_mol)
    if clean_scaffold is not None:
        clean_scaffold = Chem.RemoveHs(clean_scaffold)

    if clean_mol is None or clean_scaffold is None:
        return []
    if clean_scaffold.GetNumAtoms() == 0:
        return []

    # Try substructure match
    match = clean_mol.GetSubstructMatch(clean_scaffold)
    if match:
        return list(match)

    # If exact match fails, try with more flexible matching
    # This handles cases where bond orders might differ
    from rdkit.Chem import AllChem

    params = AllChem.AdjustQueryParameters()
    params.makeBondsGeneric = True
    params.makeAtomsGeneric = False

    try:
        query = AllChem.AdjustQueryProperties(clean_scaffold, params)
        match = clean_mol.GetSubstructMatch(query)
        if match:
            return list(match)
    except:
        pass

    # Fallback: try matching by SMILES core pattern
    # Get the core ring system from scaffold
    try:
        from rdkit.Chem.Scaffolds import MurckoScaffold

        core = MurckoScaffold.GetScaffoldForMol(clean_scaffold)
        if core and core.GetNumAtoms() > 0:
            match = clean_mol.GetSubstructMatch(core)
            if match:
                return list(match)
    except:
        pass

    return []


def mol_to_svg(mol, scaffold_mol=None, width=200, height=150, highlight_rgroups=True):
    """Convert RDKit molecule to SVG string with optional R-group highlighting."""
    if mol is None:
        return "<p>Invalid</p>"
    try:
        clean_mol = remove_dummy_atoms(mol)
        if clean_mol is None or clean_mol.GetNumAtoms() == 0:
            return "<p>Empty</p>"

        # Remove explicit hydrogens for cleaner visualization
        clean_mol = Chem.RemoveHs(clean_mol)

        AllChem.Compute2DCoords(clean_mol)

        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        # Don't show H atoms
        opts = drawer.drawOptions()
        opts.addAtomIndices = False

        # Highlight R-groups (atoms NOT in scaffold) in orange
        if highlight_rgroups and scaffold_mol is not None:
            # Use the clean mol (without Hs) for matching
            scaffold_atoms = get_scaffold_match(clean_mol, scaffold_mol)
            if scaffold_atoms:
                # All atoms not in scaffold are R-groups
                all_atoms = set(range(clean_mol.GetNumAtoms()))
                rgroup_atoms = list(all_atoms - set(scaffold_atoms))

                if rgroup_atoms:
                    # Orange highlight for R-groups
                    highlight_colors = {
                        idx: (0.92, 0.49, 0.24) for idx in rgroup_atoms
                    }  # #ea7d3d
                    drawer.DrawMolecule(
                        clean_mol,
                        highlightAtoms=rgroup_atoms,
                        highlightAtomColors=highlight_colors,
                    )
                else:
                    drawer.DrawMolecule(clean_mol)
            else:
                drawer.DrawMolecule(clean_mol)
        else:
            drawer.DrawMolecule(clean_mol)

        drawer.FinishDrawing()
        return drawer.GetDrawingText()
    except Exception as e:
        return f"<p>Error: {str(e)[:20]}</p>"


def mol_to_svg_large(mol, scaffold_mol=None, width=700, height=550):
    """Generate larger SVG for modal view."""
    return mol_to_svg(mol, scaffold_mol, width, height, highlight_rgroups=True)


def extract_molecules(chemspace, max_molecules=100):
    """Extract RDKit molecules from chemical space."""
    molecules = []

    if hasattr(chemspace, "df") and chemspace.df is not None:
        df = chemspace.df
        mol_cols = [
            col for col in df.columns if "mol" in col.lower() or "rmol" in col.lower()
        ]
        if mol_cols:
            for idx, row in df.head(max_molecules).iterrows():
                rmol = row[mol_cols[0]]
                if rmol is not None:
                    if hasattr(rmol, "mol"):
                        molecules.append((idx, rmol.mol))
                    elif hasattr(rmol, "GetMol"):
                        molecules.append((idx, rmol.GetMol()))
                    elif isinstance(rmol, Chem.Mol):
                        molecules.append((idx, rmol))

    if not molecules and hasattr(chemspace, "__iter__"):
        try:
            for idx, item in enumerate(chemspace):
                if idx >= max_molecules:
                    break
                if hasattr(item, "mol"):
                    molecules.append((idx, item.mol))
                elif hasattr(item, "GetMol"):
                    molecules.append((idx, item.GetMol()))
                elif isinstance(item, Chem.Mol):
                    molecules.append((idx, item))
        except:
            pass

    return molecules


def get_scaffold_mol(chemspace):
    """Extract scaffold molecule from chemspace or external files.

    Priority:
    1. chemspace.scaffold attribute
    2. scaffold.pkl from Node 2 (has attachment point info)
    3. ligand.sdf from Node 1 (original ligand)
    """
    # First try from chemspace object
    if hasattr(chemspace, "scaffold") and chemspace.scaffold is not None:
        scaffold = chemspace.scaffold
        if hasattr(scaffold, "mol"):
            print("  - Scaffold source: chemspace.scaffold")
            return scaffold.mol
        elif hasattr(scaffold, "GetMol"):
            print("  - Scaffold source: chemspace.scaffold")
            return scaffold.GetMol()
        elif isinstance(scaffold, Chem.Mol):
            print("  - Scaffold source: chemspace.scaffold")
            return scaffold

    # Try loading from scaffold.pkl (Node 2 output - preferred, has attachment point)
    scaffold_paths = [
        Path("scaffold.pkl"),
    ]
    for scaffold_path in scaffold_paths:
        if scaffold_path.exists():
            try:
                if cloudpickle:
                    with open(scaffold_path, "rb") as f:
                        scaffold = cloudpickle.load(f)
                else:
                    with open(scaffold_path, "rb") as f:
                        scaffold = pickle.load(f)
                if hasattr(scaffold, "mol"):
                    print(f"  - Scaffold source: {scaffold_path}")
                    return scaffold.mol
                elif hasattr(scaffold, "GetMol"):
                    print(f"  - Scaffold source: {scaffold_path}")
                    return scaffold.GetMol()
                elif isinstance(scaffold, Chem.Mol):
                    print(f"  - Scaffold source: {scaffold_path}")
                    return scaffold
            except Exception as e:
                print(f"  - Failed to load {scaffold_path}: {e}")

    # Fallback: Try loading from ligand.smi (Node 1 output - original ligand as SMILES)
    ligand_smi_paths = [
        Path("ligand.smi"),
    ]
    for ligand_path in ligand_smi_paths:
        if ligand_path.exists():
            try:
                with open(ligand_path, "r") as f:
                    smiles = f.read().strip()
                mol = Chem.MolFromSmiles(smiles)
                if mol is not None:
                    print(f"  - Scaffold source: {ligand_path}")
                    return mol
            except Exception as e:
                print(f"  - Failed to load {ligand_path}: {e}")

    print("  - No scaffold found")
    return None


def generate_molecule_cards(molecules, scaffold_mol=None):
    """Generate HTML cards for each molecule."""
    if not molecules:
        return '<p style="text-align: center; color: #888;">No molecules to display</p>'

    cards = []
    large_svgs = []

    for idx, mol in molecules:
        # Small SVG for card
        svg_small = mol_to_svg(mol, scaffold_mol, width=200, height=150)

        # Large SVG for modal
        svg_large = mol_to_svg_large(mol, scaffold_mol)
        large_svgs.append((idx, svg_large))

        # Get SMILES
        smiles = ""
        smiles_full = ""
        try:
            clean_mol = remove_dummy_atoms(mol)
            if clean_mol:
                smiles_full = Chem.MolToSmiles(clean_mol)
                smiles = (
                    smiles_full if len(smiles_full) <= 30 else smiles_full[:27] + "..."
                )
        except:
            smiles = "N/A"
            smiles_full = "N/A"

        card = f'''
        <div class="mol-card" onclick="showModal({idx})">
            <div class="mol-idx">#{idx}</div>
            <div class="mol-svg">{svg_small}</div>
            <div class="mol-smiles" title="{smiles_full}">{smiles}</div>
        </div>
        '''
        cards.append(card)

    # Generate hidden divs for large SVGs
    modal_data = []
    for idx, svg_large in large_svgs:
        mol = next((m for i, m in molecules if i == idx), None)
        smiles_full = ""
        try:
            if mol:
                clean_mol = remove_dummy_atoms(mol)
                if clean_mol:
                    smiles_full = Chem.MolToSmiles(clean_mol)
        except:
            smiles_full = "N/A"

        modal_data.append(f"""
        <div id="modal-content-{idx}" style="display:none;">
            <div class="modal-svg">{svg_large}</div>
            <div class="modal-smiles">{smiles_full}</div>
        </div>
        """)

    return "\n".join(cards), "\n".join(modal_data)


HTML_TEMPLATE = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Chemical Space Summary</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: #1a1a2e; color: #eee; padding: 20px;
        }}
        .header {{ text-align: center; padding: 20px; background: #16213e; border-radius: 12px; margin-bottom: 20px; }}
        .header h1 {{ color: #ea7d3d; }}
        .container {{ max-width: 1400px; margin: 0 auto; }}
        .stats {{ display: grid; grid-template-columns: repeat(3, 1fr); gap: 20px; margin-bottom: 20px; }}
        .stat {{ background: #16213e; padding: 30px; border-radius: 12px; text-align: center; border: 1px solid #2a2a4a; }}
        .stat .value {{ font-size: 48px; color: #ea7d3d; font-weight: bold; }}
        .stat .label {{ color: #888; font-size: 16px; margin-top: 8px; }}
        .info {{ background: #16213e; padding: 24px; border-radius: 12px; border: 1px solid #2a2a4a; }}
        .info h3 {{ color: #ea7d3d; margin-bottom: 16px; }}
        .info p {{ line-height: 1.8; color: #ccc; }}
        .next-step {{ background: #16213e; padding: 24px; border-radius: 12px; margin-top: 20px; border-left: 4px solid #ea7d3d; }}
        .next-step h4 {{ color: #ea7d3d; margin-bottom: 12px; }}

        /* Legend */
        .legend {{
            background: #16213e;
            padding: 16px 24px;
            border-radius: 12px;
            margin-bottom: 20px;
            display: flex;
            align-items: center;
            gap: 24px;
            border: 1px solid #2a2a4a;
        }}
        .legend-title {{ color: #ea7d3d; font-weight: bold; }}
        .legend-item {{ display: flex; align-items: center; gap: 8px; }}
        .legend-color {{
            width: 20px;
            height: 20px;
            border-radius: 4px;
        }}
        .legend-scaffold {{ background: #333; border: 2px solid #666; }}
        .legend-rgroup {{ background: #ea7d3d; }}

        /* Molecule Grid */
        .molecules-section {{ margin-top: 30px; }}
        .molecules-section h3 {{ color: #ea7d3d; margin-bottom: 20px; }}
        .mol-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(220px, 1fr));
            gap: 16px;
        }}
        .mol-card {{
            background: #16213e;
            border-radius: 12px;
            padding: 12px;
            border: 1px solid #2a2a4a;
            text-align: center;
            transition: transform 0.2s, border-color 0.2s;
            cursor: pointer;
        }}
        .mol-card:hover {{
            transform: translateY(-4px);
            border-color: #ea7d3d;
        }}
        .mol-idx {{
            color: #ea7d3d;
            font-weight: bold;
            margin-bottom: 8px;
        }}
        .mol-svg {{
            background: white;
            border-radius: 8px;
            padding: 8px;
            display: flex;
            justify-content: center;
            align-items: center;
            min-height: 150px;
        }}
        .mol-svg svg {{
            max-width: 100%;
            height: auto;
        }}
        .mol-smiles {{
            font-size: 11px;
            color: #888;
            margin-top: 8px;
            word-break: break-all;
            max-height: 2.4em;
            overflow: hidden;
        }}

        /* Modal */
        .modal {{
            display: none;
            position: fixed;
            z-index: 1000;
            left: 0;
            top: 0;
            width: 100%;
            height: 100%;
            background-color: rgba(0,0,0,0.8);
            justify-content: center;
            align-items: center;
        }}
        .modal.active {{
            display: flex;
        }}
        .modal-content {{
            background: #16213e;
            padding: 30px;
            border-radius: 16px;
            width: 85vw;
            max-width: 900px;
            height: 85vh;
            max-height: 85vh;
            overflow: auto;
            position: relative;
            border: 2px solid #ea7d3d;
            display: flex;
            flex-direction: column;
        }}
        .modal-close {{
            position: absolute;
            top: 10px;
            right: 15px;
            font-size: 28px;
            color: #ea7d3d;
            cursor: pointer;
            background: none;
            border: none;
            z-index: 10;
        }}
        .modal-close:hover {{
            color: #fff;
        }}
        .modal-title {{
            color: #ea7d3d;
            margin-bottom: 20px;
            text-align: center;
        }}
        .modal-svg {{
            background: white;
            border-radius: 12px;
            padding: 20px;
            display: flex;
            justify-content: center;
            align-items: center;
            margin-bottom: 20px;
            flex: 1;
            min-height: 0;
            overflow: hidden;
        }}
        .modal-svg svg {{
            max-width: 100%;
            max-height: 100%;
            height: auto;
            width: auto;
        }}
        .modal-smiles {{
            font-family: monospace;
            font-size: 12px;
            color: #ccc;
            background: #1a1a2e;
            padding: 12px;
            border-radius: 8px;
            word-break: break-all;
        }}
        /* Navigation arrows */
        .modal-nav {{
            position: absolute;
            top: 50%;
            transform: translateY(-50%);
            font-size: 36px;
            color: #ea7d3d;
            cursor: pointer;
            background: rgba(22, 33, 62, 0.9);
            border: none;
            padding: 10px 15px;
            border-radius: 8px;
            transition: background 0.2s, color 0.2s;
            z-index: 10;
        }}
        .modal-nav:hover {{
            background: #ea7d3d;
            color: #fff;
        }}
        .modal-nav.prev {{
            left: -60px;
        }}
        .modal-nav.next {{
            right: -60px;
        }}
        .modal-wrapper {{
            position: relative;
            display: flex;
            align-items: center;
        }}
        .modal-counter {{
            text-align: center;
            color: #888;
            font-size: 14px;
            margin-top: 10px;
        }}
        @media (max-width: 768px) {{
            .modal-nav.prev {{ left: 10px; }}
            .modal-nav.next {{ right: 10px; }}
            .modal-nav {{ font-size: 24px; padding: 8px 12px; }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>Chemical Space Summary</h1>
            <p>Node 5: Chemical Space Creation</p>
        </div>
        <div class="stats">
            <div class="stat">
                <div class="value">{num_molecules}</div>
                <div class="label">Generated Molecules</div>
            </div>
            <div class="stat">
                <div class="value">{has_scaffold}</div>
                <div class="label">Scaffold Loaded</div>
            </div>
            <div class="stat">
                <div class="value">{has_protein}</div>
                <div class="label">Protein Loaded</div>
            </div>
        </div>

        <div class="legend">
            <span class="legend-title">Legend:</span>
            <div class="legend-item">
                <div class="legend-color legend-scaffold"></div>
                <span>Scaffold atoms</span>
            </div>
            <div class="legend-item">
                <div class="legend-color legend-rgroup"></div>
                <span>R-groups / Linkers (highlighted)</span>
            </div>
            <span style="color: #888; margin-left: auto;">Click molecule for larger view</span>
        </div>

        <div class="molecules-section">
            <h3>Generated Molecules {display_note}</h3>
            <div class="mol-grid">
                {molecule_cards}
            </div>
        </div>

        <div class="info">
            <h3>Chemical Space Information</h3>
            <p>The chemical space has been generated by combining the scaffold with various linkers and R-groups.</p>
            <p>Total of <strong>{num_molecules}</strong> unique molecules have been created in this combinatorial library.</p>
            <p><strong>Orange highlighted atoms</strong> indicate the R-groups and linkers that were added to the scaffold.</p>
        </div>
    </div>

    <!-- Modal for enlarged view -->
    <div class="modal" id="mol-modal" onclick="closeModal(event)">
        <div class="modal-wrapper">
            <button class="modal-nav prev" onclick="event.stopPropagation(); navigateModal(-1)">&#8249;</button>
            <div class="modal-content" onclick="event.stopPropagation()">
                <button class="modal-close" onclick="closeModal()">&times;</button>
                <h3 class="modal-title" id="modal-title">Molecule #0</h3>
                <div id="modal-body"></div>
                <div class="modal-counter" id="modal-counter">1 / 100</div>
            </div>
            <button class="modal-nav next" onclick="event.stopPropagation(); navigateModal(1)">&#8250;</button>
        </div>
    </div>

    <!-- Hidden modal content -->
    <div style="display:none;">
        {modal_data}
    </div>

    <script>
        // Store molecule indices for navigation
        const molIndices = {mol_indices_json};
        let currentMolPosition = 0;

        function showModal(idx) {{
            const modal = document.getElementById('mol-modal');
            const title = document.getElementById('modal-title');
            const body = document.getElementById('modal-body');
            const counter = document.getElementById('modal-counter');
            const content = document.getElementById('modal-content-' + idx);

            if (content) {{
                // Find position in array
                currentMolPosition = molIndices.indexOf(idx);
                if (currentMolPosition === -1) currentMolPosition = 0;

                title.textContent = 'Molecule #' + idx;
                body.innerHTML = content.innerHTML;
                counter.textContent = (currentMolPosition + 1) + ' / ' + molIndices.length;
                modal.classList.add('active');
            }}
        }}

        function navigateModal(direction) {{
            currentMolPosition += direction;
            // Wrap around
            if (currentMolPosition < 0) currentMolPosition = molIndices.length - 1;
            if (currentMolPosition >= molIndices.length) currentMolPosition = 0;

            const idx = molIndices[currentMolPosition];
            const title = document.getElementById('modal-title');
            const body = document.getElementById('modal-body');
            const counter = document.getElementById('modal-counter');
            const content = document.getElementById('modal-content-' + idx);

            if (content) {{
                title.textContent = 'Molecule #' + idx;
                body.innerHTML = content.innerHTML;
                counter.textContent = (currentMolPosition + 1) + ' / ' + molIndices.length;
            }}
        }}

        function closeModal(event) {{
            if (event && event.target !== document.getElementById('mol-modal')) return;
            document.getElementById('mol-modal').classList.remove('active');
        }}

        // Keyboard navigation
        document.addEventListener('keydown', function(e) {{
            const modal = document.getElementById('mol-modal');
            if (!modal.classList.contains('active')) return;

            if (e.key === 'Escape') {{
                closeModal();
            }} else if (e.key === 'ArrowLeft') {{
                navigateModal(-1);
            }} else if (e.key === 'ArrowRight') {{
                navigateModal(1);
            }}
        }});
    </script>
</body>
</html>
"""


def main():
    input_file = sys.argv[1] if len(sys.argv) > 1 else "chemspace.pkl"
    output_file = sys.argv[2] if len(sys.argv) > 2 else "chemspace_viz.html"

    if not Path(input_file).exists():
        print(f"Error: {input_file} not found")
        sys.exit(1)

    chemspace = load_chemspace(input_file)

    # Extract information
    num_molecules = 0
    has_scaffold = "No"
    has_protein = "No"

    if hasattr(chemspace, "df") and chemspace.df is not None:
        num_molecules = len(chemspace.df)

    if hasattr(chemspace, "scaffold") and chemspace.scaffold is not None:
        has_scaffold = "Yes"

    if hasattr(chemspace, "protein") or hasattr(chemspace, "receptor"):
        has_protein = "Yes"

    # Get scaffold for highlighting
    scaffold_mol = get_scaffold_mol(chemspace)

    # Extract and visualize molecules (limit to 100 for performance)
    max_display = 100
    molecules = extract_molecules(chemspace, max_molecules=max_display)
    molecule_cards, modal_data = generate_molecule_cards(molecules, scaffold_mol)

    display_note = ""
    if num_molecules > max_display:
        display_note = f"(showing first {max_display} of {num_molecules})"

    # Get molecule indices for navigation
    import json

    mol_indices = [idx for idx, _ in molecules]

    html = HTML_TEMPLATE.format(
        num_molecules=num_molecules,
        has_scaffold=has_scaffold,
        has_protein=has_protein,
        molecule_cards=molecule_cards,
        modal_data=modal_data,
        display_note=display_note,
        mol_indices_json=json.dumps(mol_indices),
    )

    with open(output_file, "w") as f:
        f.write(html)

    print(f"Generated: {output_file}")
    print(f"  - Molecules: {num_molecules}")
    print(f"  - Displayed: {len(molecules)}")
    print(
        f"  - Scaffold: {has_scaffold} (R-group highlighting: {'enabled' if scaffold_mol else 'disabled'})"
    )
    print(f"  - Protein: {has_protein}")


if __name__ == "__main__":
    main()
