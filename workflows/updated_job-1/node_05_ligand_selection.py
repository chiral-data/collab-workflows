#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 5: Ligand Selection - Select ligands
Input: ligands.sdf (SDF files in extracted directory)
Output: Ligands_select.sdf (combined SDF file with selected ligands)
"""

import glob
import os
from rdkit import Chem

INPUT_DIR = "./constructed_library"
OUTPUT_FILE = "Ligands_select.sdf"

# For testing: select first 11 ligands (clean_drug108*.sdf)
# For production: change this pattern to select all ligands
LIGAND_PATTERN = "clean_drug108*.sdf"  # For testing
# LIGAND_PATTERN = "clean_drug*.sdf"  # For production

def main():
    """Main execution function"""
    print("=== Node 5: Select ligands ===")
    
    if not os.path.exists(INPUT_DIR):
        print(f"❌ Error: {INPUT_DIR} directory not found.")
        print("   Please run Node 3 (node_03_unpack_ligands.py) first.")
        exit(1)
    
    # Search for ligand files
    ligand_pattern = os.path.join(INPUT_DIR, LIGAND_PATTERN)
    ligand_files = glob.glob(ligand_pattern)
    
    if not ligand_files:
        print(f"❌ Error: No files found matching {ligand_pattern}.")
        exit(1)
    
    print(f"Search pattern: {ligand_pattern}")
    print(f"Number of ligand files found: {len(ligand_files)}")
    
    # Combine selected ligands into a single SDF file
    writer = Chem.SDWriter(OUTPUT_FILE)
    processed_count = 0
    
    for ligand_file in sorted(ligand_files):
        try:
            # Read molecules from SDF file
            supplier = Chem.SDMolSupplier(ligand_file)
            for mol in supplier:
                if mol is not None:
                    writer.write(mol)
                    processed_count += 1
        except Exception as e:
            print(f"Warning: Error occurred while processing {ligand_file}: {e}")
            continue
    
    writer.close()
    
    print(f"✅ Selection complete: {processed_count} ligands saved to {OUTPUT_FILE}.")


if __name__ == "__main__":
    main()

