#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
True ligand addition - Download true ligand from PDB
Input: global_params.json (pdb_id, ligand_name)
Output: outputs/true_ligand.sdf (true ligand downloaded from PDB)
Note: This runs before prepare_ligands so that true_ligand is also prepared
"""

import os
import json
import urllib.request
from rdkit import Chem

# Default values
DEFAULT_PDB_ID = "5Y7J"
DEFAULT_LIGAND_NAME = "8OL"

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Input from other nodes should be in input/ directory
# Output from this node should be in outputs/ directory
INPUT_DIR = os.path.join(SCRIPT_DIR, "input")
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")
# true_ligand.sdf will be downloaded from PDB and saved directly to outputs/
REAL_LIGAND_OUTPUT = os.path.join(OUTPUT_DIR, "true_ligand.sdf")

def load_global_params():
    """Load parameters from global_params.json"""
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Try workflow-003/global_params.json
    global_params_file = os.path.join(script_dir, "..", "global_params.json")
    if os.path.exists(global_params_file):
        try:
            with open(global_params_file, "r") as f:
                params = json.load(f)
                return params.get("pdb_id", DEFAULT_PDB_ID), params.get("ligand_name", DEFAULT_LIGAND_NAME)
        except Exception as e:
            print(f"⚠ Warning: Could not load global_params.json: {e}")
    return DEFAULT_PDB_ID, DEFAULT_LIGAND_NAME

def download_ligand_from_pdb(ligand_name, output_file):
    """
    Download ligand SDF file from RCSB PDB
    
    Args:
        ligand_name: Three-letter ligand name (e.g., "RMZ")
        output_file: Path to save the downloaded SDF file
    
    Returns:
        bool: True if download successful, False otherwise
    """
    # RCSB PDB ligand download URL
    # Try ideal SDF first (most common format)
    url = f"https://files.rcsb.org/ligands/download/{ligand_name.upper()}_ideal.sdf"
    
    print(f"Downloading ligand '{ligand_name}' from RCSB PDB...")
    print(f"URL: {url}")
    
    try:
        # Download the file
        urllib.request.urlretrieve(url, output_file)
        
        # Verify the file was downloaded and is not empty
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            # Try to read with RDKit to verify it's a valid SDF
            try:
                suppl = Chem.SDMolSupplier(output_file)
                mol = next(suppl, None)
                if mol is not None:
                    print(f"✅ Successfully downloaded ligand '{ligand_name}' from PDB")
                    print(f"   File saved to: {output_file}")
                    print(f"   File size: {os.path.getsize(output_file)} bytes")
                    return True
                else:
                    print(f"⚠ Warning: Downloaded file exists but could not be parsed as SDF")
                    # Try alternative URL format
                    return download_ligand_alternative(ligand_name, output_file)
            except Exception as e:
                print(f"⚠ Warning: Error reading SDF file: {e}")
                # Try alternative URL format
                return download_ligand_alternative(ligand_name, output_file)
        else:
            print(f"⚠ Warning: Downloaded file is empty, trying alternative URL...")
            return download_ligand_alternative(ligand_name, output_file)
            
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"⚠ Ligand '{ligand_name}' not found at ideal SDF URL, trying alternative format...")
            return download_ligand_alternative(ligand_name, output_file)
        else:
            print(f"❌ HTTP Error {e.code}: {e.reason}")
            return False
    except Exception as e:
        print(f"❌ Error downloading ligand: {e}")
        return False

def download_ligand_alternative(ligand_name, output_file):
    """
    Try alternative URL format for downloading ligand from PDB
    
    Args:
        ligand_name: Three-letter ligand name
        output_file: Path to save the downloaded SDF file
    
    Returns:
        bool: True if download successful, False otherwise
    """
    # Try alternative URL without "_ideal" suffix
    alt_url = f"https://files.rcsb.org/ligands/download/{ligand_name.upper()}.sdf"
    
    print(f"Trying alternative URL: {alt_url}")
    
    try:
        urllib.request.urlretrieve(alt_url, output_file)
        
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            try:
                suppl = Chem.SDMolSupplier(output_file)
                mol = next(suppl, None)
                if mol is not None:
                    print(f"✅ Successfully downloaded ligand '{ligand_name}' from PDB (alternative URL)")
                    print(f"   File saved to: {output_file}")
                    print(f"   File size: {os.path.getsize(output_file)} bytes")
                    return True
            except Exception:
                pass
        
        print(f"❌ Error: Could not download ligand '{ligand_name}' from PDB")
        print(f"   Tried URLs:")
        print(f"   - https://files.rcsb.org/ligands/download/{ligand_name.upper()}_ideal.sdf")
        print(f"   - https://files.rcsb.org/ligands/download/{ligand_name.upper()}.sdf")
        return False
        
    except Exception as e:
        print(f"❌ Error downloading ligand from alternative URL: {e}")
        return False

def main():
    """Main execution function"""
    print("=== True ligand addition ===")
    
    # Ensure output directory exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Load parameters from global_params.json or environment variables
    global_pdb_id, global_ligand_name = load_global_params()
    pdb_id = os.environ.get("PDB_ID") or os.environ.get("PARAM_PDB_ID") or global_pdb_id
    ligand_name = os.environ.get("LIGAND_NAME") or os.environ.get("PARAM_LIGAND_NAME") or global_ligand_name
    
    print(f"PDB ID: {pdb_id}")
    print(f"Ligand name: {ligand_name}")
    
    # Check if true_ligand.sdf already exists in destination
    if os.path.exists(REAL_LIGAND_OUTPUT):
        print(f"⚠ true_ligand.sdf already exists in {OUTPUT_DIR}, overwriting...")
    
    # Download true ligand from PDB
    print(f"\nDownloading true ligand '{ligand_name}' from RCSB PDB...")
    success = download_ligand_from_pdb(ligand_name, REAL_LIGAND_OUTPUT)
    
    if not success:
        print(f"❌ Error: Failed to download ligand '{ligand_name}' from PDB")
        print(f"   Please check that the ligand name '{ligand_name}' is correct in global_params.json")
        exit(1)
    
    # Verify the file was downloaded successfully
    if os.path.exists(REAL_LIGAND_OUTPUT):
        file_size = os.path.getsize(REAL_LIGAND_OUTPUT)
        print(f"\n✅ True ligand addition complete.")
        print(f"   File saved to: {REAL_LIGAND_OUTPUT}")
        print(f"   File size: {file_size} bytes")
        print(f"✓ Verified: true_ligand.sdf is in outputs directory")
    else:
        print(f"❌ Error: true_ligand.sdf was not downloaded successfully")
        exit(1)


if __name__ == "__main__":
    main()

