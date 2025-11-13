#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 8: Biopython - Ligand Center Identification - Identify ligand center coordinates
Input: Ligands_select.sdf, ligand.csv, 5Y7J_chain.pdb
Output: config.txt (docking configuration file)
"""

import os
import numpy as np
from Bio.PDB import PDBParser

DEFAULT_PDB_ID = "5Y7J"
LIGAND_NAME = "8OL"
INPUT_SDF = "Ligands_select.sdf"
INPUT_CSV = "ligand.csv"
INPUT_PDB = f"{DEFAULT_PDB_ID}_chain.pdb"
OUTPUT_CONFIG = "config.txt"

def main():
    """Main execution function"""
    print("=== Node 8: Biopython - Ligand center identification ===")
    
    # Get PDB ID
    pdb_id = os.environ.get("PDB_ID", DEFAULT_PDB_ID)
    input_pdb_file = f"./{pdb_id}_chain.pdb"
    
    # Check input files
    if not os.path.exists(input_pdb_file):
        print(f"❌ Error: {input_pdb_file} not found.")
        print("   Please run Node 7 (node_07_extract_chains.py) first.")
        exit(1)
    
    if not os.path.exists(INPUT_SDF):
        print(f"Warning: {INPUT_SDF} not found, but will calculate ligand center from PDB file.")
    
    if not os.path.exists(INPUT_CSV):
        print(f"Warning: {INPUT_CSV} not found, but will calculate ligand center from PDB file.")
    
    print(f"Input PDB file: {input_pdb_file}")
    print(f"Output config file: {OUTPUT_CONFIG}")
    
    # Calculate ligand center coordinates from PDB file
    center, coords_array = calculate_ligand_center(input_pdb_file)
    
    # Calculate grid size
    lig_min = coords_array.min(axis=0)
    lig_max = coords_array.max(axis=0)
    extent = lig_max - lig_min  # Molecular size in x, y, z directions
    padding = 8.0  # Padding around the molecule [Å]
    min_size = 20.0  # Minimum practical grid [Å]
    size_vec = np.maximum(extent + padding, min_size)
    
    print(f"Ligand center coordinates (x, y, z): {center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}")
    print(f"Recommended grid size (Å): size_x={size_vec[0]:.1f}, size_y={size_vec[1]:.1f}, size_z={size_vec[2]:.1f}")
    
    # Generate config file
    config_lines = [
        f"center_x = {center[0]:.3f}",
        f"center_y = {center[1]:.3f}",
        f"center_z = {center[2]:.3f}",
        "size_x   = 15",  # Use fixed value
        "size_y   = 15",
        "size_z   = 15",
        "exhaustiveness = 8",  # Search intensity
        "num_modes = 5",  # Number of output poses
        "energy_range = 4",  # Maximum energy difference between output poses
    ]
    
    with open(OUTPUT_CONFIG, "w", encoding="utf-8") as f:
        f.write("\n".join(config_lines) + "\n")
    
    print(f"✅ Docking configuration file generated: {OUTPUT_CONFIG}")


def calculate_ligand_center(pdb_file):
    """Calculate ligand binding coordinates (center coordinates)"""
    ligand_name = LIGAND_NAME
    extracted_ligand_coords = []
    
    parser = PDBParser()
    target_structure = parser.get_structure("target_protein", pdb_file)
    
    for model in target_structure:
        for chain in model:
            for residue in chain:
                if residue.get_resname() == ligand_name:
                    coords = [atom.get_coord() for atom in residue]
                    if coords:
                        extracted_ligand_coords.append(coords)
    
    # If multiple ligands found, use only the first one
    if extracted_ligand_coords:
        coords_array = np.array(extracted_ligand_coords[0])
        center = np.mean(coords_array, axis=0)
        return center, coords_array
    else:
        print(f"❌ Error: Ligand '{ligand_name}' not found in file '{pdb_file}'.")
        exit(1)


if __name__ == "__main__":
    main()

