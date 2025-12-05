#!/usr/bin/env python3

import os
import sys
import subprocess
from pathlib import Path
from pdbfixer import PDBFixer
from openmm.app import PDBFile


def main():

    input_dir = Path("./")
    output_dir = Path("./")

    # 1) Locate PDB file
    pdb_files = list(input_dir.glob("*.pdb"))
    if not pdb_files:
        print("ERROR: No .pdb file found in ./ — cannot prepare receptor")
        sys.exit(1)

    input_pdb = pdb_files[0]
    fixed_pdb = output_dir / "protein_fixed.pdb"
    output_pdbqt = output_dir / "receptor.pdbqt"

    print(f"[Node3] Input PDB: {input_pdb}")

    # 2) Clean with PDBFixer
    print(f"[Node3] Running PDBFixer on {input_pdb} ...")
    fixer = PDBFixer(filename=str(input_pdb))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens()

    with open(fixed_pdb, "w") as out_f:
        PDBFile.writeFile(fixer.topology, fixer.positions, out_f)

    print(f"[Node3] Cleaned PDB saved to: {fixed_pdb}")

    # 3) Convert with pythonsh
    print("[Node3] Converting to PDBQT with pythonsh prepare_receptor4.py ...")

    pythonsh = "/opt/mgltools_install/bin/pythonsh"
    prep_script = "/opt/mgltools_install/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor.py"

    cmd = [
        pythonsh,
        prep_script,
        "-r", str(fixed_pdb),
        "-o", str(output_pdbqt),
        "-A", "hydrogens"
        "-U", "nphs_lps_waters_nonstdres"
    ]

    try:
        subprocess.run(cmd, check=True)
    except Exception as e:
        print(f"ERROR running pythonsh: {e}")
        sys.exit(1)

    print(f"[Node3] Receptor PDBQT saved to: {output_pdbqt}")
    print("[Node3] Receptor preparation complete ✅")


if __name__ == "__main__":
    main()
