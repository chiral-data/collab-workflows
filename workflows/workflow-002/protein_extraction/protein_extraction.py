#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OpenMM PDBFixer - Clean protein - Clean up protein structure
Input: protein_input/outputs/*.pdb (PDB file from protein input)
Output: {pdb_id}_clean.pdb (cleaned up PDB file with selected chains)
"""

import os
import json
from openmm.app import PDBFile
from pdbfixer import PDBFixer
from Bio.PDB import PDBIO, PDBParser, Select

DEFAULT_PDB_ID = "5Y7J"
# Default chain ID (can be overridden by environment variable, None means all chains)
DEFAULT_CHAIN_ID = None  # Default is all chains

def load_global_params():
    """Load parameters from global_params.json"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    global_params_file = os.path.join(script_dir, "..", "global_params.json")
    if os.path.exists(global_params_file):
        try:
            with open(global_params_file, "r") as f:
                params = json.load(f)
                return params.get("pdb_id", DEFAULT_PDB_ID)
        except Exception as e:
            print(f"⚠ Warning: Could not load global_params.json: {e}")
    return DEFAULT_PDB_ID

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Input from other nodes should be in input/ directory
# Output from this node should be in outputs/ directory
INPUT_DIR = os.path.join(SCRIPT_DIR, "input")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")
# Get PDB file from protein_input/outputs/
PROTEIN_INPUT_DIR = os.path.join(SCRIPT_DIR, "..", "protein_input", "outputs")

def main():
    """Main execution function"""
    print("=== OpenMM PDBFixer - Clean up protein ===")
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load parameters from global_params.json or environment variables
    global_pdb_id = load_global_params()
    pdb_id = os.environ.get("PDB_ID") or os.environ.get("PARAM_PDB_ID") or global_pdb_id
    # If CHAIN_ID is not set, use DEFAULT_CHAIN_ID (None means all chains)
    selected_chain_id = os.environ.get("CHAIN_ID") or DEFAULT_CHAIN_ID
    
    # Get PDB file from protein_input/outputs/
    input_pdb_file = os.path.join(PROTEIN_INPUT_DIR, f"{pdb_id}.pdb")
    
    # Fallback to input/ directory if not found in protein_input
    if not os.path.exists(input_pdb_file):
        input_pdb_file = os.path.join(INPUT_DIR, f"{pdb_id}.pdb")
    
    if not os.path.exists(input_pdb_file):
        print(f"❌ Error: PDB file not found.")
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
    
    print(f"Input file: {input_pdb_file}")
    
    try:
        # First, check available chains and select if needed
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure("protein", input_pdb_file)
        
        # Get all chains in the structure
        all_chains = []
        for model in structure:
            for chain in model:
                if chain.id not in all_chains:
                    all_chains.append(chain.id)
        
        all_chains = sorted(all_chains)
        print(f"\nAvailable chains in PDB file: {', '.join(all_chains)}")
        
        # Determine which chains to keep
        if selected_chain_id is None:
            chains_to_keep = all_chains
            print(f"Using all chains: {', '.join(chains_to_keep)}")
        else:
            if selected_chain_id not in all_chains:
                print(f"❌ Error: Chain '{selected_chain_id}' not found in PDB file.")
                print(f"   Available chains: {', '.join(all_chains)}")
                exit(1)
            chains_to_keep = [selected_chain_id]
            print(f"Selected chain: {selected_chain_id}")
        
        # Generate output file name based on selected chains
        if len(chains_to_keep) == len(all_chains):
            # All chains selected
            chain_suffix = "clean"
        elif len(chains_to_keep) == 1:
            # Single chain selected
            chain_suffix = f"chain_{chains_to_keep[0]}_clean"
        else:
            # Multiple chains selected
            chain_suffix = f"chain_{''.join(chains_to_keep)}_clean"
        
        output_pdb_file = os.path.join(OUTPUT_DIR, f"{pdb_id}_{chain_suffix}.pdb")
        
        print(f"Output file: {output_pdb_file}")
        
        # Extract selected chains if needed
        temp_pdb_file = input_pdb_file
        if selected_chain_id is not None:
            # Create a temporary file with only selected chains
            import tempfile
            temp_fd, temp_pdb_file = tempfile.mkstemp(suffix=".pdb", prefix="selected_chains_")
            os.close(temp_fd)
            
            class ChainSelect(Select):
                def accept_chain(self, chain):
                    return chain.id in chains_to_keep
            
            io = PDBIO()
            io.set_structure(structure)
            io.save(temp_pdb_file, ChainSelect())
            print(f"Extracted chains {', '.join(chains_to_keep)} to temporary file")
        
        # Process with PDBFixer
        print(f"\nProcessing {temp_pdb_file} with PDBFixer...")
        fixer = PDBFixer(filename=temp_pdb_file)
        
        fixer.findMissingResidues()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.removeHeterogens(keepWater=False)  # OK if receptor only
        fixer.addMissingHydrogens()
        
        # Write fixed structure
        with open(output_pdb_file, "w") as fout:
            PDBFile.writeFile(fixer.topology, fixer.positions, fout)
        
        # Clean up temporary file if created
        if temp_pdb_file != input_pdb_file and os.path.exists(temp_pdb_file):
            os.remove(temp_pdb_file)
        
        print(f"\n✅ Cleanup complete: {output_pdb_file}")
        if selected_chain_id:
            print(f"   Contains chain(s): {', '.join(chains_to_keep)}")
    except Exception as e:
        print(f"❌ Error occurred during PDBFixer processing: {e}")
        import traceback
        traceback.print_exc()
        exit(1)


if __name__ == "__main__":
    main()

