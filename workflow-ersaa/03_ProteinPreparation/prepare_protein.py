#!/usr/bin/env python3

import sys
import subprocess
from pathlib import Path
from pdbfixer import PDBFixer
from openmm.app import PDBFile


def main():

    workdir = Path("./")

    pdb_files = list(workdir.glob("*.pdb")) + list(workdir.glob("*.cif"))
    if not pdb_files:
        print("ERROR: No PDB or CIF file found")
        sys.exit(1)

    input_structure = pdb_files[0]
    fixed_pdb = workdir / "protein_fixed.pdb"
    receptor_pdbqt = workdir / "receptor.pdbqt"

    print(f"[Node3] Input structure: {input_structure}")

    # -----------------------
    # PDBFixer (NO hydrogens)
    # -----------------------
    fixer = PDBFixer(filename=str(input_structure))

    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()

    # Keep crystallographic waters if present
    fixer.removeHeterogens(False)

    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    # IMPORTANT: DO NOT add hydrogens here
    # AutoDockTools must do this

    with open(fixed_pdb, "w") as f:
        PDBFile.writeFile(
            fixer.topology,
            fixer.positions,
            f,
            keepIds=True
        )

    print(f"[Node3] Hydrogen-free PDB written: {fixed_pdb}")

    # -----------------------
    # AutoDockTools (Vina)
    # -----------------------
    pythonsh = "/opt/mgltools_install/bin/pythonsh"
    prep_script = (
        "/opt/mgltools_install/MGLToolsPckgs/"
        "AutoDockTools/Utilities24/prepare_receptor4.py"
    )

    cmd = [
        pythonsh,
        prep_script,
        "-r", str(fixed_pdb),
        "-o", str(receptor_pdbqt),
        "-A", "checkhydrogens",
        "-U", "nphs_lps_waters_nonstdres"
    ]

    subprocess.run(cmd, check=True)

    print(f"[Node3] Receptor PDBQT saved: {receptor_pdbqt}")
    print("[Node3] Protein preparation complete ✅")


if __name__ == "__main__":
    main()
