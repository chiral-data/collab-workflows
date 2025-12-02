#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ligand Loading - Extract ligand information from PDB
Input: 5Y7J.pdb (full protein structure)
Output: 5Y7J_ligand.txt (ligand information file)
"""

import os
import json
import numpy as np
from Bio.PDB import PDBParser

DEFAULT_PDB_ID = "5Y7J"
# Default ligand name (can be overridden by environment variable or global_params.json)
DEFAULT_LIGAND_NAME = "8OL"
# Default chain ID (can be overridden by environment variable, None means any chain)
DEFAULT_CHAIN_ID = None

def load_global_params():
    """Load parameters from global_params.json"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Try workflow-003/global_params.json first, then workflow-003/global_params.json
    global_params_file = os.path.join(script_dir, "..", "global_params.json")
    if not os.path.exists(global_params_file):
        global_params_file = os.path.join(script_dir, "..", "global_params.json")
    if os.path.exists(global_params_file):
        try:
            with open(global_params_file, "r") as f:
                params = json.load(f)
                return params.get("pdb_id", DEFAULT_PDB_ID), params.get("ligand_name", DEFAULT_LIGAND_NAME)
        except Exception as e:
            print(f"⚠ Warning: Could not load global_params.json: {e}")
    return DEFAULT_PDB_ID, DEFAULT_LIGAND_NAME
# Standard amino acids
STANDARD_AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
}
# Standard nucleotides
STANDARD_NUCLEOTIDES = {"A", "T", "G", "C", "U", "DA", "DT", "DG", "DC", "DU"}
# Water molecules
WATER_NAMES = {"HOH", "WAT", "H2O"}

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Input from other nodes should be in input/ directory
# Output from this node should be in outputs/ directory
INPUT_DIR = os.path.join(SCRIPT_DIR, "input")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")
# PDB file from protein_input/outputs/
PROTEIN_INPUT_DIR = os.path.join(SCRIPT_DIR, "..", "protein_input", "outputs")
# For optional Ligands_select.sdf and ligand.csv (not required)
INPUT_LIGAND_DIR = INPUT_DIR
OUTPUT_LIGAND_LIST = os.path.join(OUTPUT_DIR, "ligand_list.txt")

def main():
    """Main execution function"""
    print("=== Ligand Loading ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load parameters from global_params.json or environment variables
    global_pdb_id, global_ligand_name = load_global_params()
    pdb_id = os.environ.get("PDB_ID") or os.environ.get("PARAM_PDB_ID") or global_pdb_id
    
    # Find PDB file from protein_input/outputs/
    input_pdb_file = os.path.join(PROTEIN_INPUT_DIR, f"{pdb_id}.pdb")
    
    # Fallback to input/ directory
    if not os.path.exists(input_pdb_file):
        input_pdb_file = os.path.join(INPUT_DIR, f"{pdb_id}.pdb")
    
    input_sdf_file = os.path.join(INPUT_LIGAND_DIR, "Ligands_select.sdf")
    input_csv_file = os.path.join(INPUT_LIGAND_DIR, "ligand.csv")
    
    # Check input files
    if not os.path.exists(input_pdb_file):
        print(f"❌ Error: {pdb_id}.pdb not found.")
        print(f"   Searched in: {PROTEIN_INPUT_DIR}")
        print(f"   Searched in: {INPUT_DIR}")
        if os.path.exists(PROTEIN_INPUT_DIR):
            files = [f for f in os.listdir(PROTEIN_INPUT_DIR) if f.endswith('.pdb')]
            if files:
                print(f"   Available PDB files in protein_input/outputs: {files}")
        if os.path.exists(INPUT_DIR):
            files = [f for f in os.listdir(INPUT_DIR) if f.endswith('.pdb')]
            if files:
                print(f"   Available PDB files in input: {files}")
        print("   Please ensure protein_input has completed.")
        exit(1)
    
    print(f"✓ Found PDB file: {input_pdb_file}")
    
    if not os.path.exists(input_sdf_file):
        print(f"Warning: {input_sdf_file} not found, but will calculate ligand center from PDB file.")
    
    if not os.path.exists(input_csv_file):
        print(f"Warning: {input_csv_file} not found, but will calculate ligand center from PDB file.")
    
    print(f"Input PDB file: {input_pdb_file}")
    print(f"Output ligand list file: {OUTPUT_LIGAND_LIST}")
    
    # Detect all ligands in PDB file
    ligands = detect_all_ligands(input_pdb_file)
    
    if not ligands:
        print("❌ Error: No ligands found in PDB file.")
        exit(1)
    
    # Display all ligands grouped by chain
    print(f"\n=== Found {len(ligands)} ligand(s) in PDB file ===")
    for i, (chain_id, ligand_name, center, num_atoms) in enumerate(ligands, 1):
        print(f"{i}. Chain: {chain_id}, Ligand: {ligand_name}")
        print(f"   Center coordinates (x, y, z): {center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}")
        print(f"   Number of atoms: {num_atoms}")
    
    # Save ligand list to file
    save_ligand_list(ligands, OUTPUT_LIGAND_LIST)
    
    # Select ligand for docking config
    # Priority: 1) Environment variables (LIGAND_NAME, PARAM_LIGAND_NAME), 2) global_params.json, 3) DEFAULT values, 4) First ligand found
    selected_ligand_name = os.environ.get("LIGAND_NAME") or os.environ.get("PARAM_LIGAND_NAME") or global_ligand_name
    # If CHAIN_ID is not set, use DEFAULT_CHAIN_ID (None means all chains)
    selected_chain_id = os.environ.get("CHAIN_ID") or DEFAULT_CHAIN_ID
    selected_ligand = None
    
    print(f"\nSearching for ligand '{selected_ligand_name}' from global_params.json...")
    
    # Try to find the specified ligand (matching ligand name, chain is optional)
    matching_ligands = []
    for chain_id, ligand_name, center, num_atoms in ligands:
        ligand_match = ligand_name.upper() == selected_ligand_name.upper()
        chain_match = (selected_chain_id is None) or (chain_id == selected_chain_id)
        
        if ligand_match:
            matching_ligands.append((chain_id, ligand_name, center, num_atoms))
            # If chain also matches, use this one
            if chain_match:
                selected_ligand = (chain_id, ligand_name, center, num_atoms)
                break
    
    # If exact match (ligand + chain) not found, use first matching ligand name
    if selected_ligand is None and matching_ligands:
        selected_ligand = matching_ligands[0]
        print(f"✓ Found ligand '{selected_ligand[1]}' (first match from {len(matching_ligands)} found)")
    
    # If specified ligand not found at all, use the first one
    if selected_ligand is None:
        print(f"\n⚠ Warning: Ligand '{selected_ligand_name}'", end="")
        if selected_chain_id:
            print(f" in chain '{selected_chain_id}'", end="")
        print(f" not found. Using first ligand (Chain: {ligands[0][0]}, Ligand: {ligands[0][1]}).")
        selected_ligand = ligands[0]
    else:
        print(f"\n✓ Using ligand '{selected_ligand[1]}' in chain '{selected_ligand[0]}' for docking configuration.")
    
    chain_id, ligand_name, center, num_atoms = selected_ligand
    
    print(f"\n✓ Selected ligand: Chain {chain_id}, Ligand {ligand_name}")
    print(f"  Center coordinates (x, y, z): {center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}")
    print(f"\n✅ Ligand list saved to: {OUTPUT_LIGAND_LIST}")


def is_ligand(residue_name):
    """Check if a residue name represents a ligand (non-standard amino acid/nucleotide)"""
    residue_name = residue_name.strip().upper()
    if residue_name in STANDARD_AMINO_ACIDS:
        return False
    if residue_name in STANDARD_NUCLEOTIDES:
        return False
    if residue_name in WATER_NAMES:
        return False
    return True


def detect_all_ligands(pdb_file):
    """
    Detect all ligands in PDB file and calculate their center coordinates.
    Ligands are identified by chain ID and ligand name.
    
    Returns:
        list: List of tuples (chain_id, ligand_name, center_coordinates, num_atoms)
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("target_protein", pdb_file)
    
    # Dictionary: {(chain_id, ligand_name): list of atom coordinates}
    ligands_dict = {}
    
    for model in structure:
        for chain in model:
            chain_id = chain.id
            for residue in chain:
                residue_name = residue.get_resname().strip()
                
                # Check if this is a ligand
                if is_ligand(residue_name):
                    # Get all atom coordinates for this residue
                    coords = []
                    for atom in residue:
                        try:
                            coord = atom.get_coord()
                            coords.append(coord)
                        except:
                            continue
                    
                    if coords:
                        # Store coordinates for this ligand in this chain
                        key = (chain_id, residue_name)
                        if key not in ligands_dict:
                            ligands_dict[key] = []
                        ligands_dict[key].extend(coords)
    
    # Calculate center coordinates for each ligand in each chain
    ligands = []
    for (chain_id, ligand_name), coords_list in ligands_dict.items():
        if coords_list:
            coords_array = np.array(coords_list)
            center = np.mean(coords_array, axis=0)
            num_atoms = len(coords_list)
            ligands.append((chain_id, ligand_name, center, num_atoms))
    
    # Sort by chain ID first, then by ligand name for consistent output
    ligands.sort(key=lambda x: (x[0], x[1]))
    
    return ligands


def save_ligand_list(ligands, output_file):
    """Save list of ligands with their center coordinates to a file"""
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("Ligand List from PDB File\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Total number of ligands found: {len(ligands)}\n\n")
        
        # Group by chain for better readability
        chains_dict = {}
        for chain_id, ligand_name, center, num_atoms in ligands:
            if chain_id not in chains_dict:
                chains_dict[chain_id] = []
            chains_dict[chain_id].append((ligand_name, center, num_atoms))
        
        for chain_id in sorted(chains_dict.keys()):
            f.write(f"Chain {chain_id}:\n")
            f.write("-" * 60 + "\n")
            for ligand_name, center, num_atoms in chains_dict[chain_id]:
                f.write(f"  Ligand: {ligand_name}\n")
                f.write(f"    Center coordinates (x, y, z): {center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}\n")
                f.write(f"    Number of atoms: {num_atoms}\n")
                f.write("\n")
            f.write("\n")


if __name__ == "__main__":
    main()

