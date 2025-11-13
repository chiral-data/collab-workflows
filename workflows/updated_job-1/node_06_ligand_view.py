#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 6: Ligand View - Export ligand information to CSV
Input: ligands_select.sdf (or Ligands_select.sdf)
Output: ligand.csv (CSV file with ligand information)
"""

import os
import csv
from rdkit import Chem
from rdkit.Chem import Descriptors

INPUT_FILE = "Ligands_select.sdf"
OUTPUT_FILE = "ligand.csv"

def main():
    """Main execution function"""
    print("=== Node 6: Export ligand information to CSV ===")
    
    if not os.path.exists(INPUT_FILE):
        print(f"❌ Error: {INPUT_FILE} not found.")
        print("   Please run Node 5 (node_05_ligand_selection.py) first.")
        exit(1)
    
    print(f"Input file: {INPUT_FILE}")
    print(f"Output file: {OUTPUT_FILE}")
    
    # Read molecules from SDF file
    supplier = Chem.SDMolSupplier(INPUT_FILE)
    
    # Write to CSV file
    with open(OUTPUT_FILE, "w", newline="", encoding="utf-8") as csvfile:
        fieldnames = [
            "index",
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
        for i, mol in enumerate(supplier):
            if mol is not None:
                try:
                    # Calculate molecular descriptors
                    smiles = Chem.MolToSmiles(mol)
                    mw = Descriptors.MolWt(mol)
                    logp = Descriptors.MolLogP(mol)
                    num_atoms = mol.GetNumAtoms()
                    num_rings = Descriptors.RingCount(mol)
                    num_rotatable_bonds = Descriptors.NumRotatableBonds(mol)
                    
                    writer.writerow({
                        "index": i + 1,
                        "smiles": smiles,
                        "molecular_weight": f"{mw:.2f}",
                        "logp": f"{logp:.2f}",
                        "num_atoms": num_atoms,
                        "num_rings": num_rings,
                        "num_rotatable_bonds": num_rotatable_bonds,
                    })
                    count += 1
                except Exception as e:
                    print(f"Warning: Error occurred while processing molecule {i+1}: {e}")
                    continue
    
    print(f"✅ CSV export complete: {count} ligand information saved to {OUTPUT_FILE}.")


if __name__ == "__main__":
    main()

