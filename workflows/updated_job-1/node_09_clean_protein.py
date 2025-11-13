#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 9: OpenMM PDBFixer - Clean protein - Clean up protein structure
Input: 5Y7J_chain.pdb
Output: 5Y7J_clean.pdb (cleaned up PDB file)
"""

import os
from openmm.app import PDBFile
from pdbfixer import PDBFixer

DEFAULT_PDB_ID = "5Y7J"

def main():
    """Main execution function"""
    print("=== Node 9: OpenMM PDBFixer - Clean up protein ===")
    
    # Get PDB ID
    pdb_id = os.environ.get("PDB_ID", DEFAULT_PDB_ID)
    input_pdb_file = f"./{pdb_id}_chain.pdb"
    # Original code output 5Y7J_AB_chains_fixed.pdb, but following workflow image,
    # output 5Y7J_clean.pdb
    output_pdb_file = f"./{pdb_id}_clean.pdb"
    # Also create fixed version for compatibility with original code
    output_fixed_pdb_file = f"./{pdb_id}_AB_chains_fixed.pdb"
    
    print(f"Input file: {input_pdb_file}")
    print(f"Output file: {output_pdb_file}")
    
    if not os.path.exists(input_pdb_file):
        print(f"❌ Error: {input_pdb_file} not found.")
        print("   Please run Node 7 (node_07_extract_chains.py) first.")
        exit(1)
    
    try:
        print(f"Processing {input_pdb_file} with PDBFixer...")
        fixer = PDBFixer(filename=input_pdb_file)
        
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)  # OK if receptor only
        fixer.addMissingHydrogens()
        
        # Write fixed structure
        with open(output_pdb_file, "w") as fout:
            PDBFile.writeFile(fixer.topology, fixer.positions, fout)
        
        # Also create fixed version for compatibility with original code
        import shutil
        shutil.copy(output_pdb_file, output_fixed_pdb_file)
        
        print(f"✅ Cleanup complete: {output_pdb_file}")
        print(f"   Also created {output_fixed_pdb_file} for compatibility.")
    except Exception as e:
        print(f"❌ Error occurred during PDBFixer processing: {e}")
        exit(1)


if __name__ == "__main__":
    main()

