#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract the 2TK ligand from a PDB file and save it as 2TK.sdf
(Preserves original coordinates from the PDB)
"""

from Bio.PDB import PDBParser, PDBIO, Select
from rdkit import Chem
import os


class Ligand2TKSelect(Select):
    """Select only 2TK residues"""

    def accept_residue(self, residue):
        return residue.get_resname() == "2TK"


def extract_2TK_from_pdb(pdb_path, output_sdf="2TK.sdf"):
    """Extract the 2TK ligand from a PDB file and save as SDF (keeping original conformation)"""
    if not os.path.exists(pdb_path):
        print(f"File not found: {pdb_path}")
        return

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    # Save only the 2TK residue to a temporary PDB
    temp_pdb = "2TK_temp.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(temp_pdb, Ligand2TKSelect())

    # Read molecule from the PDB without generating new coordinates
    mol = Chem.MolFromPDBFile(temp_pdb, removeHs=False)
    if mol is None:
        print("RDKit could not parse the extracted 2TK PDB.")
        return

    # Save to SDF (preserves PDB coordinates)
    writer = Chem.SDWriter(output_sdf)
    writer.write(mol)
    writer.close()

    os.remove(temp_pdb)
    print(f"Extracted 2TK ligand and saved as: {output_sdf}")


if __name__ == "__main__":
    extract_2TK_from_pdb("4OHU_A_NAD.pdb")
