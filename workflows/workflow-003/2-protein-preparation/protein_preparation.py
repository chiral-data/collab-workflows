#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

from Bio.PDB import PDBIO, PDBParser, Select

# ------------------------------------------------------------------------------
# 2. Extract Chain A + NAD
# ------------------------------------------------------------------------------


def extract_chain_a_with_nad(pdb_file, output_dir="outputs"):
    """Extract Chain A with NAD cofactor (exclude 2TK)"""
    print("\n=== Extracting Chain A with NAD cofactor ===")
    os.makedirs(output_dir, exist_ok=True)

    cofactor_name = "NAD"
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    chains = [chain.id for model in structure for chain in model]
    unique_chains = sorted(set(chains))
    print(f"Chains present: {', '.join(unique_chains)}")

    if "A" not in unique_chains:
        print("⚠ Chain A not found. Using original file.")
        return pdb_file

    basename = os.path.basename(pdb_file)
    output_pdb_file = f"{output_dir}/{os.path.splitext(basename)[0]}_A_NAD.pdb"

    class AChainAndNADSelect(Select):
        def accept_chain(self, chain):
            return chain.id == "A"

        def accept_residue(self, residue):
            if residue.get_parent().id == "A":
                return True
            return residue.get_resname() == cofactor_name

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb_file, AChainAndNADSelect())
    print(f"✅ Extracted: {output_pdb_file}")

    return output_pdb_file


# ------------------------------------------------------------------------------
# Entry point
# ------------------------------------------------------------------------------

if __name__ == "__main__":
    import shutil
    # pdb_id = "4OHU"
    pdb_id = os.getenv("PARAM_PDB_ID")
    # Extract Chain A with NAD (exclude 2TK)
    pdb_file = f"inputs/{pdb_id}.pdb"
    output_pdb_file = extract_chain_a_with_nad(pdb_file)

    # Also copy original PDB to outputs for downstream nodes
    os.makedirs("outputs", exist_ok=True)
    shutil.copy(pdb_file, f"outputs/{pdb_id}.pdb")
    print(f"✅ Copied original: outputs/{pdb_id}.pdb")
