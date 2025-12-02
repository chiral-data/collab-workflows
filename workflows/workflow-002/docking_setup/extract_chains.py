#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chain extraction - Extract chains and real ligand
Input: 5Y7J.pdb
Output: 5Y7J_chain.pdb (PDB file with extracted chains), real_ligand.sdf (extracted ligand)

Environment variables:
    PDB_ID: PDB ID (default: "5Y7J")
    LIGAND_NAME: Ligand name to extract (default: "8OL")
    CHAIN_ID: Chain ID(s) to extract, comma-separated for multiple chains (e.g., "A" or "A,B")
              Default is "A". The ligand will be extracted only from the selected chains.
"""

import os
import json
import tempfile
import subprocess
import numpy as np
from Bio.PDB import PDBIO, PDBParser, Select

DEFAULT_PDB_ID = "5Y7J"
LIGAND_NAME = "8OL"
# Default chain ID (can be overridden by environment variable CHAIN_ID)
# Can be a single chain (e.g., "A") or multiple chains separated by commas (e.g., "A,B")
# None means all chains (default behavior)
DEFAULT_CHAIN_ID = None  # Default is all chains

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
                return params.get("pdb_id", DEFAULT_PDB_ID), params.get("ligand_name", LIGAND_NAME)
        except Exception as e:
            print(f"⚠ Warning: Could not load global_params.json: {e}")
    return DEFAULT_PDB_ID, LIGAND_NAME
# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Input from other nodes should be in input/ directory
# Output from this node should be in outputs/ directory
INPUT_DIR = os.path.join(SCRIPT_DIR, "input")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")

def main():
    """Main execution function"""
    print("=== Chain extraction ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load parameters from global_params.json or environment variables
    global_pdb_id, global_ligand_name = load_global_params()
    pdb_id = os.environ.get("PDB_ID") or os.environ.get("PARAM_PDB_ID") or global_pdb_id
    ligand_name = os.environ.get("LIGAND_NAME") or os.environ.get("PARAM_LIGAND_NAME") or global_ligand_name
    
    # Find PDB file from protein_input/outputs/
    input_pdb_file = os.path.join(SCRIPT_DIR, "..", "protein_input", "outputs", f"{pdb_id}.pdb")
    
    # Fallback to input/ directory
    if not os.path.exists(input_pdb_file):
        input_pdb_file = os.path.join(INPUT_DIR, f"{pdb_id}.pdb")
    
    output_pdb_file = os.path.join(OUTPUT_DIR, f"{pdb_id}_chain.pdb")
    output_ligand_sdf = os.path.join(OUTPUT_DIR, "real_ligand.sdf")
    output_config_file = os.path.join(OUTPUT_DIR, "config.txt")
    
    print(f"Input file: {input_pdb_file}")
    print(f"Output file: {output_pdb_file}")
    
    if not os.path.exists(input_pdb_file):
        print(f"❌ Error: {pdb_id}.pdb not found.")
        print(f"   Searched in: {os.path.join(SCRIPT_DIR, '..', 'protein_input', 'outputs')}")
        print(f"   Searched in: {INPUT_DIR}")
        if os.path.exists(os.path.join(SCRIPT_DIR, "..", "protein_input", "outputs")):
            files = [f for f in os.listdir(os.path.join(SCRIPT_DIR, "..", "protein_input", "outputs")) if f.endswith('.pdb')]
            if files:
                print(f"   Available PDB files in protein_input/outputs: {files}")
        if os.path.exists(INPUT_DIR):
            files = [f for f in os.listdir(INPUT_DIR) if f.endswith('.pdb')]
            if files:
                print(f"   Available PDB files in input: {files}")
        print("   Please ensure protein_input has completed.")
        exit(1)
    
    # Read PDB file
    parser = PDBParser()
    structure = parser.get_structure("protein", input_pdb_file)
    
    # Check chains
    chains = [chain.id for model in structure for chain in model]
    unique_chains = sorted(list(set(chains)))
    print(f"\nChains present in PDB file: {', '.join(unique_chains)}")
    
    # Get chain IDs to extract from environment variable
    # If CHAIN_ID is not set, use DEFAULT_CHAIN_ID (None means all chains)
    chain_id_env = os.environ.get("CHAIN_ID") or DEFAULT_CHAIN_ID
    
    # Parse chain IDs (can be comma-separated like "A,B" or single like "A")
    if chain_id_env:
        selected_chains = [c.strip().upper() for c in str(chain_id_env).split(",")]
        selected_chains = [c for c in selected_chains if c]  # Remove empty strings
    else:
        # If not specified, extract all chains (default behavior)
        selected_chains = unique_chains
        print(f"No CHAIN_ID specified. Will extract all chains: {', '.join(selected_chains)}")
    
    # Validate that all selected chains exist in the PDB file
    missing_chains = [c for c in selected_chains if c not in unique_chains]
    if missing_chains:
        print(f"❌ Error: The following chains are not present in the PDB file: {', '.join(missing_chains)}")
        print(f"   Available chains: {', '.join(unique_chains)}")
        print(f"   Please set CHAIN_ID environment variable to one or more of the available chains.")
        print(f"   Example: export CHAIN_ID=A or export CHAIN_ID=A,B")
        exit(1)
    
    print(f"Selected chains to extract: {', '.join(selected_chains)}")
    
    # Class to select specified chains and ligand
    class ChainAndLigandSelect(Select):
        def __init__(self, chain_ids, ligand_name):
            self.chain_ids = set(chain_ids)
            self.ligand_name = ligand_name
        
        def accept_chain(self, chain):
            return chain.id in self.chain_ids
        
        def accept_residue(self, residue):
            return (
                residue.get_parent().id in self.chain_ids
                or residue.get_resname().strip().upper() == self.ligand_name.upper()
            )
    
    # Extract selected chains
    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb_file, ChainAndLigandSelect(selected_chains, ligand_name))
    print(f"✅ Chain extraction complete: {output_pdb_file}")
    
    # Extract real ligand coordinates from selected chains and save as SDF
    ligand_center = extract_ligand_to_sdf(input_pdb_file, ligand_name, output_ligand_sdf, selected_chains)
    
    # Generate config.txt file with ligand center coordinates
    if ligand_center is not None:
        generate_config_file(ligand_center, output_config_file)
        print(f"✅ Config file generated: {output_config_file}")
    else:
        print(f"⚠ Warning: Could not calculate ligand center coordinates. Config file not generated.")


def extract_ligand_to_sdf(pdb_file, ligand_name, output_sdf_file, selected_chains):
    """
    Extract coordinates of a specific ligand from selected chains in PDB file and save as SDF.
    Only extracts coordinates - no bond information or aromaticity processing.
    Also calculates and returns the center coordinates of the ligand.
    
    Args:
        pdb_file: Path to input PDB file
        ligand_name: Name of the ligand residue (from global_params.json)
        output_sdf_file: Path to output SDF file
        selected_chains: List of chain IDs to extract ligand from
    
    Returns:
        tuple: (center_x, center_y, center_z) or None if ligand not found
    """
    try:
        # Parse PDB file
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", pdb_file)
        
        # Calculate ligand center coordinates from PDB structure
        ligand_coords = []
        for model in structure:
            for chain in model:
                if chain.id in selected_chains:
                    for residue in chain:
                        residue_name = residue.get_resname().strip().upper()
                        if residue_name == ligand_name.upper():
                            for atom in residue:
                                try:
                                    coord = atom.get_coord()
                                    ligand_coords.append(coord)
                                except:
                                    continue
        
        # Calculate center coordinates
        ligand_center = None
        if ligand_coords:
            coords_array = np.array(ligand_coords)
            center = np.mean(coords_array, axis=0)
            ligand_center = (center[0], center[1], center[2])
            print(f"✓ Calculated ligand center coordinates: ({ligand_center[0]:.3f}, {ligand_center[1]:.3f}, {ligand_center[2]:.3f})")
        else:
            print(f"⚠ Warning: Ligand '{ligand_name}' not found in selected chains: {', '.join(selected_chains)}")
            return None
        
        # Create a Select class to extract only the specified ligand from selected chains
        class LigandSelect(Select):
            def __init__(self, chain_ids, ligand_name):
                self.chain_ids = set(chain_ids)
                self.ligand_name = ligand_name.upper()
            
            def accept_residue(self, residue):
                residue_name = residue.get_resname().strip().upper()
                chain_id = residue.get_parent().id
                # Only extract ligand if it's in one of the selected chains
                return (residue_name == self.ligand_name and 
                        chain_id in self.chain_ids)
        
        # Save ligand to temporary PDB file (coordinates only)
        temp_fd, temp_pdb_file = tempfile.mkstemp(suffix=".pdb", prefix="ligand_")
        os.close(temp_fd)
        
        io = PDBIO()
        io.set_structure(structure)
        io.save(temp_pdb_file, LigandSelect(selected_chains, ligand_name))
        
        # Check if temporary file has any content (ligand was found)
        if os.path.getsize(temp_pdb_file) == 0:
            print(f"⚠ Warning: Ligand '{ligand_name}' not found in selected chains: {', '.join(selected_chains)}")
            os.remove(temp_pdb_file)
            return None
        
        # Convert PDB to SDF (coordinates only, no processing)
        try:
            abs_temp_pdb = os.path.abspath(temp_pdb_file)
            abs_output_sdf = os.path.abspath(output_sdf_file)
            
            # Convert PDB to SDF directly (coordinates only)
            result = subprocess.run(
                ["obabel", abs_temp_pdb, "-O", abs_output_sdf],
                capture_output=True,
                text=True,
                timeout=60
            )
            
            if result.returncode != 0:
                print(f"⚠ Warning: obabel PDB->SDF conversion failed")
                print(f"   Return code: {result.returncode}")
                print(f"   stderr: {result.stderr}")
                print(f"   stdout: {result.stdout}")
                os.remove(temp_pdb_file)
                return ligand_center
            
            if not os.path.exists(output_sdf_file) or os.path.getsize(output_sdf_file) == 0:
                print(f"⚠ Warning: SDF file was not created or is empty.")
                os.remove(temp_pdb_file)
                return ligand_center
                
        except subprocess.TimeoutExpired:
            print(f"⚠ Warning: obabel conversion timed out.")
            os.remove(temp_pdb_file)
            return ligand_center
        except Exception as e:
            print(f"⚠ Warning: Error running obabel: {e}")
            os.remove(temp_pdb_file)
            return ligand_center
        finally:
            # Clean up temporary file
            if os.path.exists(temp_pdb_file):
                os.remove(temp_pdb_file)
        
        print(f"✅ Extracted ligand '{ligand_name}' coordinates from chain(s) {', '.join(selected_chains)} and saved to: {os.path.basename(output_sdf_file)}")
        
        return ligand_center
        
    except Exception as e:
        print(f"⚠ Warning: Error extracting ligand '{ligand_name}': {e}")
        import traceback
        traceback.print_exc()
        return None


def generate_config_file(ligand_center, output_config_file):
    """
    Generate docking configuration file with ligand center coordinates.
    
    Args:
        ligand_center: Tuple of (center_x, center_y, center_z)
        output_config_file: Path to output config.txt file
    """
    center_x, center_y, center_z = ligand_center
    
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
    
    print(f"   Center coordinates (x, y, z): {center_x:.3f}, {center_y:.3f}, {center_z:.3f}")


if __name__ == "__main__":
    main()

