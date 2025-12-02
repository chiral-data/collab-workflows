#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ligand Center Identification - Generate docking configuration file
Input: 5Y7J_ligand.txt (from ligand_loading), 5Y7J_chain.pdb (from extract_chains)
Output: config.txt (docking configuration file)
"""

import os
import json

DEFAULT_PDB_ID = "5Y7J"
# Default ligand name (can be overridden by environment variable or global_params.json)
DEFAULT_LIGAND_NAME = "8OL"
# Default chain ID (can be overridden by environment variable)
# None means all chains (default behavior)
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
# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Input from other nodes should be in input/ directory
# Output from this node should be in outputs/ directory
INPUT_DIR = os.path.join(SCRIPT_DIR, "input")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")
OUTPUT_CONFIG = os.path.join(OUTPUT_DIR, "config.txt")

def main():
    """Main execution function"""
    print("=== Ligand Center Identification ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load parameters from global_params.json or environment variables
    global_pdb_id, global_ligand_name = load_global_params()
    
    # Get chain selection
    # If CHAIN_ID is not set, use DEFAULT_CHAIN_ID (None means all chains)
    selected_chain_id = os.environ.get("CHAIN_ID") or DEFAULT_CHAIN_ID
    # ligand_list.txt is output from ligand_loading in the same node, so check outputs/ first
    input_ligand_list_file = os.path.join(OUTPUT_DIR, "ligand_list.txt")
    
    print(f"\nInput ligand list file: {input_ligand_list_file}")
    print(f"Output config file: {OUTPUT_CONFIG}")
    
    if not os.path.exists(input_ligand_list_file):
        print(f"\n❌ Error: ligand_list.txt not found.")
        print(f"   Searched in: {OUTPUT_DIR}")
        print("   Please run ligand_loading first.")
        exit(1)
    
    try:
        # Generate config.txt from ligand_list.txt
        generate_config_file(input_ligand_list_file, OUTPUT_CONFIG, selected_chain_id, global_ligand_name)
        print(f"\n✅ Successfully generated docking configuration file: {OUTPUT_CONFIG}")
        
    except Exception as e:
        print(f"\n❌ Error occurred during config file generation: {e}")
        exit(1)


def generate_config_file(ligand_list_file, output_config_file, selected_chain_id, default_ligand_name=None):
    """
    Generate docking configuration file from ligand list file.
    Reads the ligand center coordinates from ligand_list.txt that binds to the selected chain.
    
    Args:
        ligand_list_file: Path to ligand_list.txt file
        output_config_file: Path to output config.txt file
        selected_chain_id: Chain ID selected in protein extraction (should match the chain in the input PDB file)
        default_ligand_name: Default ligand name from global_params.json
    """
    # Get ligand selection from environment variables, global_params.json, or use defaults
    selected_ligand_name = os.environ.get("LIGAND_NAME") or os.environ.get("PARAM_LIGAND_NAME") or default_ligand_name or DEFAULT_LIGAND_NAME
    
    print(f"\nSearching for ligand '{selected_ligand_name}' from global_params.json...")
    
    # Read ligand center coordinates from ligand_list.txt
    center_x = None
    center_y = None
    center_z = None
    current_chain = None
    current_ligand = None
    matching_ligands = []  # Store all matching ligands
    
    # Priority 1: Find all ligands matching the specified ligand name
    with open(ligand_list_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            
            # Check for chain header
            if line.startswith("Chain "):
                current_chain = line.split()[1].rstrip(":")
                continue
            
            # Check for ligand name
            if line.startswith("Ligand: "):
                current_ligand = line.split(":", 1)[1].strip()
                continue
            
            # Check for center coordinates line
            if line.startswith("Center coordinates"):
                if ":" in line and current_ligand:
                    ligand_match = current_ligand.upper() == selected_ligand_name.upper()
                    if ligand_match:
                        coords_str = line.split(":", 1)[1].strip()
                        coords = [float(x.strip()) for x in coords_str.split(",")]
                        if len(coords) == 3:
                            matching_ligands.append((current_chain, current_ligand, coords))
    
    # Use first matching ligand (from global_params.json)
    if matching_ligands:
        if selected_chain_id:
            # Try to find in selected chain first
            for chain, ligand, coords in matching_ligands:
                if chain == selected_chain_id:
                    center_x, center_y, center_z = coords
                    print(f"✓ Found ligand '{ligand}' in chain '{chain}' (matching global_params.json)")
                    break
        
        # If not found in selected chain, use first matching ligand
        if center_x is None:
            chain, ligand, coords = matching_ligands[0]
            center_x, center_y, center_z = coords
            if len(matching_ligands) > 1:
                print(f"✓ Found {len(matching_ligands)} matching ligands. Using first one: '{ligand}' in chain '{chain}'")
            else:
                print(f"✓ Found ligand '{ligand}' in chain '{chain}' (matching global_params.json)")
    
    # Priority 2: If no matching ligand found, use any ligand in the selected chain
    if center_x is None and selected_chain_id:
        print(f"⚠ Ligand '{selected_ligand_name}' not found. Searching for any ligand in chain '{selected_chain_id}'...")
        
        with open(ligand_list_file, "r", encoding="utf-8") as f:
            current_chain = None
            for line in f:
                line = line.strip()
                
                if line.startswith("Chain "):
                    current_chain = line.split()[1].rstrip(":")
                    continue
                
                if line.startswith("Ligand: "):
                    current_ligand = line.split(":", 1)[1].strip()
                    continue
                
                if line.startswith("Center coordinates") and current_chain == selected_chain_id:
                    if ":" in line:
                        coords_str = line.split(":", 1)[1].strip()
                        coords = [float(x.strip()) for x in coords_str.split(",")]
                        if len(coords) == 3:
                            center_x, center_y, center_z = coords
                            print(f"✓ Using ligand '{current_ligand}' from chain '{current_chain}'")
                            break
    
    # Priority 3: If still no ligand found, use first ligand in the file (fallback)
    if center_x is None:
        print(f"⚠ No ligand found. Using first ligand from file as fallback...")
        with open(ligand_list_file, "r", encoding="utf-8") as f:
            current_chain = None
            for line in f:
                line = line.strip()
                
                if line.startswith("Chain "):
                    current_chain = line.split()[1].rstrip(":")
                    continue
                
                if line.startswith("Ligand: "):
                    current_ligand = line.split(":", 1)[1].strip()
                    continue
                
                if line.startswith("Center coordinates"):
                    if ":" in line:
                        coords_str = line.split(":", 1)[1].strip()
                        coords = [float(x.strip()) for x in coords_str.split(",")]
                        if len(coords) == 3:
                            center_x, center_y, center_z = coords
                            print(f"✓ Using first ligand '{current_ligand}' from chain '{current_chain}' (fallback)")
                            break
    
    if center_x is None or center_y is None or center_z is None:
        raise ValueError("Could not parse center coordinates from ligand_list.txt")
    
    print(f"   Center coordinates (x, y, z): {center_x:.3f}, {center_y:.3f}, {center_z:.3f}")
    
    # Generate config file
    config_lines = [
        f"center_x = {center_x:.3f}",
        f"center_y = {center_y:.3f}",
        f"center_z = {center_z:.3f}",
        "size_x   = 15",  # Use fixed value
        "size_y   = 15",
        "size_z   = 15",
        "exhaustiveness = 8",  # Search intensity
        "num_modes = 5",  # Number of output poses
        "energy_range = 4",  # Maximum energy difference between output poses
    ]
    
    with open(output_config_file, "w", encoding="utf-8") as f:
        f.write("\n".join(config_lines) + "\n")


if __name__ == "__main__":
    main()

