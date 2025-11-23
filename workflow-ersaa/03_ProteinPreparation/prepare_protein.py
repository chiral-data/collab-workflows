#!/usr/bin/env python3

import sys
import os
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import subprocess
from pathlib import Path

def prepare_receptor(input_pdb: Path, output_pdbqt: Path):
    print(f"[Node3] Fixing protein: {input_pdb}")

    fixer = PDBFixer(filename=str(input_pdb))
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)

    fixed = Path("protein_fixed.pdb")
    with open(fixed, "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)

    print("[Node3] Converting to PDBQT...")
    cmd = [
        "python2", "/opt/mgltools/bin/prepare_receptor4.py",
        "-r", str(fixed),
        "-o", str(output_pdbqt),
        "-A", "hydrogens"
    ]
    subprocess.run(cmd, check=True)

    print(f"[Node3] Saved receptor to {output_pdbqt}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python prepare_receptor.py <INPUT_PDB> <OUTPUT_PDBQT>")
        sys.exit(2)

    input_pdb = Path(sys.argv[1])
    output_pdbqt = Path(sys.argv[2])
    output_pdbqt.parent.mkdir(parents=True, exist_ok=True)

    prepare_receptor(input_pdb, output_pdbqt)
