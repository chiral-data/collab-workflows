#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 10: PDB2PQR - Apply AMBER charges - Apply AMBER charges
Input: 5Y7J_clean.pdb
Output: 5Y7J_amber.pdb (PQR file with AMBER charges applied, or PDB file)
"""

import os
import subprocess

DEFAULT_PDB_ID = "5Y7J"

def main():
    """Main execution function"""
    print("=== Node 10: PDB2PQR - Apply AMBER charges ===")
    
    # Get PDB ID
    pdb_id = os.environ.get("PDB_ID", DEFAULT_PDB_ID)
    input_pdb_file = f"./{pdb_id}_clean.pdb"
    output_pqr_file = f"./{pdb_id}_amber.pqr"
    output_pdb_file = f"./{pdb_id}_amber.pdb"  # Also output PDB file following workflow image
    
    print(f"Input file: {input_pdb_file}")
    print(f"Output PQR file: {output_pqr_file}")
    print(f"Output PDB file: {output_pdb_file}")
    
    if not os.path.exists(input_pdb_file):
        print(f"❌ Error: {input_pdb_file} not found.")
        print("   Please run Node 9 (node_09_clean_protein.py) first.")
        exit(1)
    
    try:
        # Run PDB2PQR
        print(f"Applying AMBER charges to {input_pdb_file} with PDB2PQR...")
        subprocess.run(
            ["pdb2pqr", "--ff=AMBER", input_pdb_file, output_pqr_file],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        print(f"✅ Created PQR file with AMBER charges: {output_pqr_file}")
        
        # Generate PDB file from PQR file following workflow image
        # Read PQR file content and output in PDB format
        try:
            from Bio.PDB import PDBIO, PDBParser
            
            parser = PDBParser()
            structure = parser.get_structure("protein", output_pqr_file)
            
            io = PDBIO()
            io.set_structure(structure)
            io.save(output_pdb_file)
            
            print(f"✅ Generated PDB file from PQR file: {output_pdb_file}")
        except Exception as e:
            print(f"Warning: Error occurred during PQR to PDB conversion: {e}")
            print("   Please use PQR file as is.")
        
    except subprocess.CalledProcessError as e:
        print(f"❌ Error occurred during PDB2PQR execution")
        if e.stderr:
            print(f"Error details: {e.stderr.decode()}")
        exit(1)
    except FileNotFoundError:
        print("❌ Error: pdb2pqr not found in system path.")
        exit(1)


if __name__ == "__main__":
    main()

