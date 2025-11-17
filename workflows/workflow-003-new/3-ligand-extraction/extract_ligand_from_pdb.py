#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract the ligand from a PDB file and save it as ligand_name.sdf
(Preserves original coordinates from the PDB)
"""

import os

from Bio.PDB import PDBIO, PDBParser, Select
from rdkit import Chem


class LigandSelect(Select):
    """Select only residues"""

    def accept_residue(self, residue):
        ligand_name = os.getenv("PARAM_LIGAND_NAME")
        return residue.get_resname() == ligand_name
        # return residue.get_resname() == "2TK"


def extract_ligand_from_pdb(pdb_path, ligand_name):
    """Extract the ligand from a PDB file and save as SDF (keeping original conformation)"""
    if not os.path.exists(pdb_path):
        print(f"File not found: {pdb_path}")
        return

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)

    # Save only the ligand residue to a temporary PDB
    temp_pdb = f"{ligand_name}_temp.pdb"
    io = PDBIO()
    io.set_structure(structure)
    io.save(temp_pdb, LigandSelect())

    # Read molecule from the PDB without generating new coordinates
    mol = Chem.MolFromPDBFile(temp_pdb, removeHs=False)
    if mol is None:
        print(f"RDKit could not parse the extracted {ligand_name} PDB.")
        return

    # Save to SDF (preserves PDB coordinates)
    writer = Chem.SDWriter(f"{ligand_name}.sdf")
    writer.write(mol)
    writer.close()

    os.remove(temp_pdb)
    print(f"Extracted {ligand_name} ligand and saved as: {ligand_name}.sdf")


if __name__ == "__main__":
    # pdb_id = "4OHU"
    pdb_id = os.getenv("PARAM_PDB_ID")
    # ligand_name = "2TK"
    ligand_name = os.getenv("PARAM_LIGAND_NAME")
    extract_ligand_from_pdb(f"{pdb_id}_A_NAD.pdb", ligand_name)
