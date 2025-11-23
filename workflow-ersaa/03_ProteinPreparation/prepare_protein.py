import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile
import os
import subprocess

def prepare_protein(input_pdb, output_pdbqt):
    print(f"Fixing protein: {input_pdb}")
    fixer = PDBFixer(filename=input_pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.4)
    fixed_pdb = "protein_fixed.pdb"
    PDBFile.writeFile(fixer.topology, fixer.positions, open(fixed_pdb, 'w'))

    print("Converting protein to PDBQT format...")
    cmd = [
        "python2", "/opt/mgltools/bin/prepare_receptor4.py",
        "-r", fixed_pdb,
        "-o", output_pdbqt,
        "-A", "hydrogens"
    ]
    subprocess.run(cmd, check=True)
    print(f"Protein prepared and saved to {output_pdbqt}")

if __name__ == "__main__":
    input_pdb = sys.argv[1]
    output_pdbqt = sys.argv[2]
    os.makedirs(os.path.dirname(output_pdbqt), exist_ok=True)
    prepare_protein(input_pdb, output_pdbqt)
