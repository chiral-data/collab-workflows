#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ligand Selection - Select ligands
Input: ligands.sdf (SDF files in extracted directory)
Output: selected_compounds/ (directory with individual SDF files for selected ligands)
"""

import glob
import os
import random
import shutil

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Input from other nodes should be in input/ directory
# Output from this node should be in outputs/ directory
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs", "selected_compounds")
# Get constructed_library from input/constructed_library (copied by run.sh from ligand_unpack)
INPUT_DIR = os.path.join(SCRIPT_DIR, "input", "constructed_library")
# Fallback to direct path from ligand_unpack/outputs/constructed_library
if not os.path.exists(INPUT_DIR):
    INPUT_DIR = os.path.join(SCRIPT_DIR, "..", "ligand_unpack", "outputs", "constructed_library")
# Fallback to outputs/constructed_library for backward compatibility
if not os.path.exists(INPUT_DIR):
    INPUT_DIR = os.path.join(SCRIPT_DIR, "outputs", "constructed_library")

# Default number of ligands to select randomly
DEFAULT_NUM_LIGANDS = 9

def main():
    """Main execution function"""
    print("=== Select ligands (random selection) ===")
    
    if not os.path.exists(INPUT_DIR):
        print(f"❌ Error: {INPUT_DIR} directory not found.")
        print("   Please run unpack_ligands first.")
        exit(1)
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Get number of ligands to select (from environment variable or use default)
    num_ligands = int(os.environ.get("NUM_LIGANDS", DEFAULT_NUM_LIGANDS))
    
    # Search for all SDF files
    ligand_pattern = os.path.join(INPUT_DIR, "*.sdf")
    all_ligand_files = glob.glob(ligand_pattern)
    
    if not all_ligand_files:
        print(f"❌ Error: No SDF files found in {INPUT_DIR}.")
        print("   Please run unpack_ligands first.")
        exit(1)
    
    print(f"Total number of ligand files found: {len(all_ligand_files)}")
    print(f"Number of ligands to select: {num_ligands}")
    
    # Randomly select ligand files
    if num_ligands >= len(all_ligand_files):
        print(f"⚠ Warning: Requested {num_ligands} ligands, but only {len(all_ligand_files)} available.")
        print(f"   Using all available ligands.")
        selected_files = all_ligand_files
    else:
        selected_files = random.sample(all_ligand_files, num_ligands)
    
    # Sort selected files for consistent output
    selected_files = sorted(selected_files)
    
    print(f"Selected {len(selected_files)} ligand file(s):")
    for i, f in enumerate(selected_files[:5], 1):  # Show first 5
        print(f"  {i}. {os.path.basename(f)}")
    if len(selected_files) > 5:
        print(f"  ... and {len(selected_files) - 5} more")
    
    # Copy selected ligand files to output directory
    copied_count = 0
    
    for ligand_file in selected_files:
        try:
            # Get the filename
            filename = os.path.basename(ligand_file)
            base_name, ext = os.path.splitext(filename)

            if ext.lower() != ".sdf":
                print(f"Warning: Skipping non-SDF file {filename}")
                continue

            # Remove any additional extension fragments such as '.cdx'
            sanitized_base = base_name.split(".")[0] if "." in base_name else base_name

            # Remove leading non-alphanumeric characters (e.g., '--_-CAMPHOR')
            sanitized_base = sanitized_base.lstrip("-_+.")

            if not sanitized_base:
                sanitized_base = base_name.replace(".", "_") or "ligand"

            output_file = os.path.join(OUTPUT_DIR, f"{sanitized_base}.sdf")

            # Avoid overwriting files with same sanitized name
            counter = 1
            while os.path.exists(output_file):
                output_file = os.path.join(OUTPUT_DIR, f"{sanitized_base}_{counter}.sdf")
                counter += 1
            
            # Copy the file
            shutil.copy2(ligand_file, output_file)
            copied_count += 1
        except Exception as e:
            print(f"Warning: Error occurred while copying {os.path.basename(ligand_file)}: {e}")
            continue
    
    print(f"\n✅ Selection complete: {copied_count} ligand file(s) saved to {OUTPUT_DIR}.")


if __name__ == "__main__":
    main()

