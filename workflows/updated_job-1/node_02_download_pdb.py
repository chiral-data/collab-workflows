#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 2: Download File - Download PDB file
Input: pdb.id (environment variable or default value)
Output: 5Y7J.pdb (based on PDB ID)
"""

import os
import urllib.request

# Default PDB ID
DEFAULT_PDB_ID = "5Y7J"

def main():
    """Main execution function"""
    print("=== Node 2: Download PDB file ===")
    
    # Get PDB ID (from environment variable or default value)
    pdb_id = os.environ.get("PDB_ID", DEFAULT_PDB_ID)
    pdb_file = f"./{pdb_id}.pdb"
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    
    print(f"PDB ID: {pdb_id}")
    print(f"Download URL: {url}")
    print(f"Output file: {pdb_file}")
    
    try:
        urllib.request.urlretrieve(url, pdb_file)
        print(f"✅ Download complete: {pdb_file}")
    except Exception as e:
        print(f"❌ Error occurred during download: {e}")
        exit(1)


if __name__ == "__main__":
    main()

