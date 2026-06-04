#!/usr/bin/env python3
"""
Node 02: Prepare Structures

1. Receptor: PDBFixer cleans → OpenBabel converts to PDBQT
2. Pocket:   Parse native ligand 3D coords → centroid + extent → pocket_config.txt
3. Native ligand: OpenBabel → native_ligand.pdbqt + native_ligand.mol2
4. Screening library: SMILES → OpenBabel 3D → xTB optimise → Meeko PDBQT
"""

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


# ── Receptor preparation ──────────────────────────────────────────────────────

def prepare_receptor(pdb_path, out_pdbqt):
    """PDBFixer → fixed.pdb, then OpenBabel → PDBQT."""
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile

    print("  Running PDBFixer ...", flush=True)
    fixer = PDBFixer(filename=str(pdb_path))
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)

    fixed_pdb = Path("receptor_fixed.pdb")
    with open(fixed_pdb, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    print(f"  Wrote {fixed_pdb}", flush=True)

    print("  Converting to PDBQT with OpenBabel ...", flush=True)
    _run(["obabel", str(fixed_pdb), "-O", str(out_pdbqt), "-xr", "-xh"],
         "receptor PDBQT conversion")
    print(f"  Wrote {out_pdbqt}", flush=True)


# ── Pocket grid box ───────────────────────────────────────────────────────────

def _coords_from_sdf(sdf_path):
    """Return Nx3 numpy array of heavy-atom coords from SDF, or None if 2D/unreadable."""
    import numpy as np
    from rdkit import Chem

    mol = Chem.MolFromMolFile(str(sdf_path), removeHs=True, sanitize=False)
    if mol is None or mol.GetNumConformers() == 0:
        return None
    conf = mol.GetConformer()
    coords = np.array([
        [conf.GetAtomPosition(i).x,
         conf.GetAtomPosition(i).y,
         conf.GetAtomPosition(i).z]
        for i in range(mol.GetNumAtoms())
    ])
    if np.std(coords[:, 2]) < 0.01:   # flat → 2D file
        return None
    return coords


def _coords_from_receptor_hetatm(receptor_pdb, ligand_res_name):
    """Extract HETATM coords for ligand_res_name from the receptor PDB (crystal binding pose)."""
    import numpy as np
    coords = []
    for line in Path(receptor_pdb).read_text().splitlines():
        if not line.startswith("HETATM"):
            continue
        res = line[17:20].strip()
        if res != ligand_res_name:
            continue
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            coords.append((x, y, z))
        except ValueError:
            continue
    return np.array(coords) if coords else None


def compute_pocket_box(native_sdf_path, padding, receptor_pdb="receptor.pdb"):
    """Return (center_x, center_y, center_z, size_x, size_y, size_z).

    Priority for coords:
    1. native_ligand.sdf if it has 3D coordinates
    2. HETATM records in receptor.pdb (crystal binding pose — preferred for 1OKL)
    3. OpenBabel --gen3d fallback (last resort; pose may not match binding site)
    """
    import numpy as np

    coords = _coords_from_sdf(native_sdf_path)

    if coords is None:
        print("  native_ligand.sdf is 2D — extracting crystal coords from receptor.pdb", flush=True)
        # Infer ligand residue name from validation_report.json if available
        lig_id = None
        report_path = Path("validation_report.json")
        if report_path.exists():
            import json
            lig_id = json.load(open(report_path)).get("native_ligand_id")
        if lig_id:
            coords = _coords_from_receptor_hetatm(receptor_pdb, lig_id)
            if coords is not None and len(coords) > 0:
                print(f"  Extracted {len(coords)} atoms for ligand {lig_id} from receptor.pdb", flush=True)

    if coords is None or len(coords) == 0:
        print("  WARNING: falling back to OpenBabel 3D generation for pocket box", flush=True)
        tmp_sdf = Path("native_ligand_3d.sdf")
        _run(["obabel", str(native_sdf_path), "-O", str(tmp_sdf), "--gen3d"],
             "native ligand 3D generation")
        coords = _coords_from_sdf(tmp_sdf)
        if coords is None:
            raise ValueError("Could not obtain 3D coordinates for the native ligand")

    center = coords.mean(axis=0)
    extent = coords.max(axis=0) - coords.min(axis=0) + 2 * padding

    cx, cy, cz = float(center[0]), float(center[1]), float(center[2])
    sx, sy, sz = float(extent[0]), float(extent[1]), float(extent[2])
    print(f"  Pocket center: ({cx:.2f}, {cy:.2f}, {cz:.2f})", flush=True)
    print(f"  Box size:      ({sx:.2f}, {sy:.2f}, {sz:.2f})", flush=True)
    return cx, cy, cz, sx, sy, sz


def write_pocket_config(cx, cy, cz, sx, sy, sz, out_path):
    lines = [
        f"center_x = {cx:.3f}",
        f"center_y = {cy:.3f}",
        f"center_z = {cz:.3f}",
        f"size_x = {sx:.3f}",
        f"size_y = {sy:.3f}",
        f"size_z = {sz:.3f}",
    ]
    Path(out_path).write_text("\n".join(lines) + "\n")
    print(f"  Wrote {out_path}", flush=True)


# ── Native ligand ─────────────────────────────────────────────────────────────

def prepare_native_ligand(native_sdf_path, out_pdbqt, out_mol2):
    print("  Converting native ligand to PDBQT ...", flush=True)
    _sdf_to_pdbqt(native_sdf_path, out_pdbqt)
    print("  Converting native ligand to MOL2 ...", flush=True)
    _run(["obabel", str(native_sdf_path), "-O", str(out_mol2)],
         "native ligand MOL2 conversion")
    print(f"  Wrote {out_pdbqt} and {out_mol2}", flush=True)


# ── Screening library ─────────────────────────────────────────────────────────

def generate_3d_from_smiles(smiles_path, out_sdf):
    """SMILES → 3D SDF via OpenBabel --gen3d."""
    print(f"  Generating 3D conformers from {smiles_path} ...", flush=True)
    _run(["obabel", str(smiles_path), "-ismi", "-osdf", "--gen3d",
          "-O", str(out_sdf)],
         "3D conformer generation")
    count = _count_sdf_mols(out_sdf)
    print(f"  Generated {count} 3D molecules → {out_sdf}", flush=True)
    return count


def xtb_optimize_sdf(input_sdf, output_sdf):
    """Optimise each molecule with GFN2-xTB; fall back to input geometry on failure."""
    from rdkit import Chem

    suppl = Chem.SDMolSupplier(str(input_sdf), removeHs=False, sanitize=False)
    writer = Chem.SDWriter(str(output_sdf))
    ok, failed = 0, 0

    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{i}"

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                mol_in = os.path.join(tmpdir, "input.sdf")
                w = Chem.SDWriter(mol_in)
                w.write(mol)
                w.close()

                res = subprocess.run(
                    ["xtb", "input.sdf", "--opt", "--gfn2", "-P", "1", "--silent"],
                    cwd=tmpdir, capture_output=True, text=True, timeout=180
                )
                opt_sdf = os.path.join(tmpdir, "xtbopt.sdf")
                if res.returncode == 0 and os.path.exists(opt_sdf):
                    opt_mol = Chem.MolFromMolFile(opt_sdf, removeHs=False, sanitize=False)
                    if opt_mol is not None:
                        opt_mol.SetProp("_Name", name)
                        writer.write(opt_mol)
                        print(f"    {name}: xTB optimised", flush=True)
                        ok += 1
                        continue
        except Exception as e:
            print(f"    {name}: xTB failed ({e}), using OpenBabel geometry", flush=True)

        writer.write(mol)
        failed += 1

    writer.close()
    print(f"  xTB: {ok} optimised, {failed} fell back to OpenBabel geometry", flush=True)


def screening_library_to_pdbqt(sdf_path, out_pdbqt):
    """Multi-mol SDF → single multi-model PDBQT with Meeko, fallback to obabel."""
    from rdkit import Chem

    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
    mols = [m for m in suppl if m is not None]
    if not mols:
        print("  WARNING: no molecules to convert to PDBQT", flush=True)
        return 0

    # Try Meeko for each molecule; fall back to obabel per-mol on failure
    pdbqt_blocks = []
    for mol in mols:
        block = _mol_to_pdbqt_meeko(mol)
        if block is None:
            block = _mol_to_pdbqt_obabel(mol)
        if block:
            pdbqt_blocks.append(block)

    Path(out_pdbqt).write_text("\n".join(pdbqt_blocks))
    print(f"  Wrote {len(pdbqt_blocks)} molecules to {out_pdbqt}", flush=True)
    return len(pdbqt_blocks)


def _mol_to_pdbqt_meeko(mol):
    """Return PDBQT string for mol via Meeko, or None on failure."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
        from meeko import MoleculePreparation

        mol = Chem.RWMol(mol)
        Chem.SanitizeMol(mol)
        if mol.GetNumConformers() == 0:
            AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        AllChem.MMFFOptimizeMolecule(mol)

        preparator = MoleculePreparation()
        mol_setups = preparator.prepare(mol)
        if not mol_setups:
            return None

        from meeko import PDBQTWriterLegacy
        pdbqt_string, _, _ = PDBQTWriterLegacy.write_string(mol_setups[0])
        return pdbqt_string
    except Exception:
        return None


def _mol_to_pdbqt_obabel(mol):
    """Return PDBQT string for mol via obabel, or None on failure."""
    from rdkit import Chem
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            in_sdf = os.path.join(tmpdir, "mol.sdf")
            out_pdbqt = os.path.join(tmpdir, "mol.pdbqt")
            w = Chem.SDWriter(in_sdf)
            w.write(mol)
            w.close()
            res = subprocess.run(
                ["obabel", in_sdf, "-O", out_pdbqt, "-xh"],
                capture_output=True, text=True, timeout=30
            )
            if res.returncode == 0 and os.path.exists(out_pdbqt):
                return Path(out_pdbqt).read_text()
    except Exception:
        pass
    return None


# ── Helpers ───────────────────────────────────────────────────────────────────

def _run(cmd, label):
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print(f"ERROR: {label} failed", flush=True)
        print(res.stderr, flush=True)
        sys.exit(1)


def _sdf_to_pdbqt(sdf_path, pdbqt_path):
    _run(["obabel", str(sdf_path), "-O", str(pdbqt_path), "-xh"],
         f"SDF→PDBQT for {sdf_path}")


def _count_sdf_mols(sdf_path):
    try:
        return Path(sdf_path).read_text().count("$$$$")
    except Exception:
        return 0


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Node 02: Prepare Structures")
    parser.add_argument("--receptor",      default="receptor.pdb")
    parser.add_argument("--native-ligand", default="native_ligand.sdf")
    parser.add_argument("--ligands",       default="validated_ligands.smiles")
    parser.add_argument("--box-padding",   type=float, default=5.0,
                        help="Angstroms to add on each side of the native ligand bounding box")
    args = parser.parse_args()

    receptor_path  = Path(args.receptor)
    native_lig     = Path(args.native_ligand)
    ligands_smiles = Path(args.ligands)

    for p in (receptor_path, native_lig, ligands_smiles):
        if not p.exists():
            print(f"ERROR: input not found: {p}", flush=True)
            sys.exit(1)

    report = {}

    # 1. Receptor
    print("\n--- Step 1: Prepare receptor ---", flush=True)
    prepare_receptor(receptor_path, "receptor.pdbqt")
    report["receptor_pdbqt"] = "receptor.pdbqt"

    # 2. Pocket box
    print("\n--- Step 2: Compute pocket grid box ---", flush=True)
    cx, cy, cz, sx, sy, sz = compute_pocket_box(native_lig, args.box_padding,
                                                receptor_pdb=str(receptor_path))
    write_pocket_config(cx, cy, cz, sx, sy, sz, "pocket_config.txt")
    report["pocket"] = {"center": [cx, cy, cz], "size": [sx, sy, sz],
                        "padding_A": args.box_padding}

    # 3. Native ligand
    print("\n--- Step 3: Prepare native ligand ---", flush=True)
    prepare_native_ligand(native_lig, "native_ligand.pdbqt", "native_ligand.mol2")
    report["native_ligand_pdbqt"] = "native_ligand.pdbqt"
    report["native_ligand_mol2"]  = "native_ligand.mol2"

    # 4. Screening library: SMILES → 3D → xTB → PDBQT
    print("\n--- Step 4: Prepare screening library ---", flush=True)
    raw_3d_sdf = Path("ligands_3d_raw.sdf")
    opt_3d_sdf = Path("ligands_3d_opt.sdf")

    n_raw = generate_3d_from_smiles(ligands_smiles, raw_3d_sdf)

    print("  Running xTB geometry optimisation ...", flush=True)
    xtb_optimize_sdf(raw_3d_sdf, opt_3d_sdf)

    n_pdbqt = screening_library_to_pdbqt(opt_3d_sdf, "optimized_screening_library.pdbqt")
    report["screening_library"] = {
        "smiles_count": n_raw,
        "pdbqt_count":  n_pdbqt,
        "pdbqt_file":   "optimized_screening_library.pdbqt",
    }

    # Summary
    with open("prepare_report.json", "w") as f:
        json.dump(report, f, indent=2)
    print("\nWrote prepare_report.json", flush=True)
    print(json.dumps(report, indent=2), flush=True)
    print("\nNode 02 completed", flush=True)


if __name__ == "__main__":
    main()
