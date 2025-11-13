#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 7: Biopython - Chain extraction - Extract chains
Input: 5Y7J.pdb
Output: 5Y7J_chain.pdb (PDB file with extracted chains)
"""

import os
from Bio.PDB import PDBIO, PDBParser, Select

DEFAULT_PDB_ID = "5Y7J"
LIGAND_NAME = "8OL"

def main():
    """Main execution function"""
    print("=== Node 7: Biopython - Chain extraction ===")
    
    # Get PDB ID
    pdb_id = os.environ.get("PDB_ID", DEFAULT_PDB_ID)
    input_pdb_file = f"./{pdb_id}.pdb"
    output_pdb_file = f"./{pdb_id}_chain.pdb"
    
    print(f"Input file: {input_pdb_file}")
    print(f"Output file: {output_pdb_file}")
    
    if not os.path.exists(input_pdb_file):
        print(f"❌ Error: {input_pdb_file} not found.")
        print("   Please run Node 2 (node_02_download_pdb.py) first.")
        exit(1)
    
    # Read PDB file
    parser = PDBParser()
    structure = parser.get_structure("protein", input_pdb_file)
    
    # Check chains
    chains = [chain.id for model in structure for chain in model]
    unique_chains = sorted(list(set(chains)))
    print(f"Chains present in PDB file: {', '.join(unique_chains)}")
    
    # Check if both chain A and B exist
    if "A" in unique_chains and "B" in unique_chains:
        print("Chains A and B detected. Extracting these chains.")
        
        # Class to select AB chains and ligand
        class ABChainAndLigandSelect(Select):
            def accept_chain(self, chain):
                return chain.id in ["A", "B"]
            
            def accept_residue(self, residue):
                return (
                    residue.get_parent().id in ["A", "B"]
                    or residue.get_resname() == LIGAND_NAME
                )
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(output_pdb_file, ABChainAndLigandSelect())
        print(f"✅ Chain extraction complete: {output_pdb_file}")
    else:
        print("Chains A and B not detected. Using original file as is.")
        import shutil
        shutil.copy(input_pdb_file, output_pdb_file)
        print(f"✅ File copied: {output_pdb_file}")


if __name__ == "__main__":
    main()

