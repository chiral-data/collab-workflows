#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ligand View - Export ligand information to CSV
Input: selected_compounds/ (directory with individual SDF files)
Output: ligand.csv (CSV file with ligand information)
"""

import os
import csv
import glob
from rdkit import Chem
from rdkit.Chem import Descriptors

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Read prepared ligands from outputs/selected_compounds (from prepare_ligands)
# Output from this node should be in outputs/ directory
INPUT_DIR = os.path.join(SCRIPT_DIR, "outputs", "selected_compounds")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")
OUTPUT_FILE = os.path.join(OUTPUT_DIR, "ligand.csv")

def main():
    """Main execution function"""
    print("=== Ligand View ===")
    
    if not os.path.exists(INPUT_DIR):
        print(f"❌ Error: {INPUT_DIR} directory not found.")
        print("   Please run prepare_ligands first.")
        exit(1)
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Find all SDF files in the selected_compounds directory
    # Include both selected ligands and true_ligand
    sdf_pattern = os.path.join(INPUT_DIR, "*.sdf")
    sdf_files = sorted(glob.glob(sdf_pattern))
    
    if not sdf_files:
        print(f"❌ Error: No SDF files found in {INPUT_DIR}.")
        print("   Please run prepare_ligands first.")
        exit(1)
    
    # Count real_ligand for information
    real_ligand_count = len([f for f in sdf_files if os.path.basename(f) == "real_ligand.sdf"])
    library_ligand_count = len(sdf_files) - real_ligand_count
    
    if real_ligand_count > 0:
        print(f"   Found {real_ligand_count} true ligand(s)")
    if library_ligand_count > 0:
        print(f"   Found {library_ligand_count} selected ligand(s)")
    
    print(f"Input directory: {INPUT_DIR}")
    print(f"Number of SDF files found: {len(sdf_files)}")
    print(f"Output file: {OUTPUT_FILE}")
    
    # Write to CSV file
    with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as csvfile:
        fieldnames = [
            "index",
            "filename",
            "smiles",
            "molecular_weight",
            "logp",
            "num_atoms",
            "num_rings",
            "num_rotatable_bonds",
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        
        count = 0
        for i, sdf_file in enumerate(sdf_files):
            try:
                # Read molecule from SDF file
                supplier = Chem.SDMolSupplier(sdf_file)
                mol = None
                for m in supplier:
                    if m is not None:
                        mol = m
                        break  # Use the first valid molecule
                
                if mol is None:
                    print(f"Warning: No valid molecule found in {os.path.basename(sdf_file)}")
                    continue
                
                # Calculate molecular descriptors
                smiles = Chem.MolToSmiles(mol)
                mw = Descriptors.MolWt(mol)
                logp = Descriptors.MolLogP(mol)
                num_atoms = mol.GetNumAtoms()
                num_rings = Descriptors.RingCount(mol)
                num_rotatable_bonds = Descriptors.NumRotatableBonds(mol)
                
                writer.writerow({
                    "index": i + 1,
                    "filename": os.path.basename(sdf_file),
                    "smiles": smiles,
                    "molecular_weight": f"{mw:.2f}",
                    "logp": f"{logp:.2f}",
                    "num_atoms": num_atoms,
                    "num_rings": num_rings,
                    "num_rotatable_bonds": num_rotatable_bonds,
                })
                count += 1
                
            except Exception as e:
                print(f"Warning: Error occurred while processing {os.path.basename(sdf_file)}: {e}")
                continue
    
    print(f"✅ CSV export complete: {count} ligand information saved to {OUTPUT_FILE}.")


if __name__ == "__main__":
    main()

