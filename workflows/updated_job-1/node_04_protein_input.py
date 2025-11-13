#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 4: Protein input - Verify protein input file
Input: pdb.id (environment variable or default value)
Output: 5Y7J.pdb (verify already downloaded file)
"""

import os

DEFAULT_PDB_ID = "5Y7J"

def main():
    """Main execution function"""
    print("=== Node 4: Verify protein input file ===")
    
    # Get PDB ID
    pdb_id = os.environ.get("PDB_ID", DEFAULT_PDB_ID)
    pdb_file = f"./{pdb_id}.pdb"
    
    print(f"PDB ID: {pdb_id}")
    print(f"File to verify: {pdb_file}")
    
    if os.path.exists(pdb_file):
        file_size = os.path.getsize(pdb_file)
        print(f"✅ File exists: {pdb_file} ({file_size / 1024:.2f} KB)")
    else:
        print(f"❌ Error: {pdb_file} not found.")
        print("   Please run Node 2 (node_02_download_pdb.py) first.")
        exit(1)


if __name__ == "__main__":
    main()

