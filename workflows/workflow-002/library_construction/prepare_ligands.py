#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Ligand Preparation - Add hydrogens, assign charges, generate 3D structures, and minimize energy
Input: 
  - ligand_selection/outputs/selected_compounds/*.sdf (SDF files from ligand selection)
  - true_ligand_addition/outputs/true_ligand.sdf (true ligand from true ligand addition)
Output: outputs/selected_compounds/*.sdf (prepared SDF files with 3D structures)
"""

import os
import glob
import shutil
import subprocess
import time
import signal
from contextlib import contextmanager
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign

# Timeout settings (in seconds)
MAX_PROCESSING_TIME_PER_LIGAND = 60  # 1 minute per ligand
MAX_STEP_TIME = 30  # 30 seconds per step

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Input from other nodes should be in input/ directory
INPUT_DIR = os.path.join(SCRIPT_DIR, "input")
# Output from this node should be in outputs/ directory
OUTPUT_BASE_DIR = os.path.join(SCRIPT_DIR, "outputs")
# Get ligands from two sources in input/ directory:
# 1. Selected compounds from ligand_selection (copied to input/selected_compounds/ by run.sh)
SELECTED_COMPOUNDS_INPUT_DIR = os.path.join(INPUT_DIR, "selected_compounds")
# Fallback to direct path from ligand_selection/outputs/
if not os.path.exists(SELECTED_COMPOUNDS_INPUT_DIR):
    SELECTED_COMPOUNDS_INPUT_DIR = os.path.join(SCRIPT_DIR, "..", "ligand_selection", "outputs", "selected_compounds")
# 2. True ligand from true_ligand_addition (copied to input/ by run.sh)
TRUE_LIGAND_FILE = os.path.join(INPUT_DIR, "true_ligand.sdf")
# Fallback to direct path from true_ligand_addition/outputs/
if not os.path.exists(TRUE_LIGAND_FILE):
    TRUE_LIGAND_FILE = os.path.join(SCRIPT_DIR, "..", "true_ligand_addition", "outputs", "true_ligand.sdf")
# Output directory for prepared ligands (explicitly in outputs/selected_compounds/)
OUTPUT_DIR = os.path.join(OUTPUT_BASE_DIR, "selected_compounds")


class ProcessingTimeoutError(Exception):
    """Custom timeout exception for ligand processing"""
    pass


@contextmanager
def timeout_context(seconds):
    """Context manager for timeout handling (Unix/Linux/macOS only)"""
    if not hasattr(signal, 'SIGALRM'):
        # Windows doesn't support SIGALRM, just yield without timeout
        yield
        return
    
    def timeout_handler(signum, frame):
        raise ProcessingTimeoutError(f"Operation timed out after {seconds} seconds")
    
    # Set the signal handler
    old_handler = signal.signal(signal.SIGALRM, timeout_handler)
    signal.alarm(seconds)
    
    try:
        yield
    finally:
        # Restore the old handler and cancel the alarm
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)

def main():
    """Main execution function"""
    print("=== Ligand Preparation ===")
    
    # Ensure output base directory exists
    os.makedirs(OUTPUT_BASE_DIR, exist_ok=True)
    # Ensure output directory exists (outputs/selected_compounds/)
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Verify that OUTPUT_DIR is in outputs/ directory (not in root)
    if not OUTPUT_DIR.startswith(OUTPUT_BASE_DIR):
        print(f"❌ Error: OUTPUT_DIR must be in outputs/ directory")
        print(f"   OUTPUT_DIR: {OUTPUT_DIR}")
        print(f"   OUTPUT_BASE_DIR: {OUTPUT_BASE_DIR}")
        exit(1)
    
    # Collect SDF files from two sources and copy to output directory
    sdf_files = []
    
    # Source 1: Selected compounds from ligand_selection (from input/selected_compounds/)
    if os.path.exists(SELECTED_COMPOUNDS_INPUT_DIR):
        selected_pattern = os.path.join(SELECTED_COMPOUNDS_INPUT_DIR, "*.sdf")
        selected_files = sorted(glob.glob(selected_pattern))
        
        # Copy selected compounds to output directory (outputs/selected_compounds/)
        for selected_file in selected_files:
            filename = os.path.basename(selected_file)
            output_file = os.path.join(OUTPUT_DIR, filename)
            try:
                shutil.copy2(selected_file, output_file)
                sdf_files.append(output_file)
            except Exception as e:
                print(f"⚠ Warning: Could not copy {filename}: {e}")
        
        print(f"✓ Copied {len(selected_files)} ligand(s) from {SELECTED_COMPOUNDS_INPUT_DIR} to {OUTPUT_DIR}")
    else:
        print(f"⚠ Warning: {SELECTED_COMPOUNDS_INPUT_DIR} directory not found.")
        print("   Please run ligand_selection first.")
    
    # Source 2: True ligand from true_ligand_addition
    if os.path.exists(TRUE_LIGAND_FILE):
        # Copy true_ligand.sdf to output directory
        output_true_ligand = os.path.join(OUTPUT_DIR, "true_ligand.sdf")
        try:
            shutil.copy2(TRUE_LIGAND_FILE, output_true_ligand)
            sdf_files.append(output_true_ligand)
            print(f"✓ Added true_ligand.sdf from true_ligand_addition")
        except Exception as e:
            print(f"⚠ Warning: Could not copy true_ligand.sdf: {e}")
    else:
        print(f"⚠ Warning: {TRUE_LIGAND_FILE} not found.")
        print("   Please run true_ligand_addition first.")
    
    if not sdf_files:
        print(f"❌ Error: No SDF files found.")
        print("   Please run ligand_selection and true_ligand_addition first.")
        exit(1)
    
    print(f"\nTotal number of SDF files to process: {len(sdf_files)}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"Building ligand library from selected compounds and true ligand...")
    
    # Process each SDF file
    successful_count = 0
    failed_count = 0
    skipped_count = 0
    
    for i, sdf_file in enumerate(sdf_files, 1):
        base_name = os.path.basename(sdf_file)
        print(f"\n[{i}/{len(sdf_files)}] Processing {base_name}...")
        
        # Record start time for this ligand
        start_time = time.time()
        
        try:
            # Check if we've exceeded the time limit before starting
            elapsed = time.time() - start_time
            if elapsed >= MAX_PROCESSING_TIME_PER_LIGAND:
                print(f"  ⏱ Skipping {base_name}: Time limit reached before processing")
                skipped_count += 1
                continue
            
            # Step 1: Add hydrogens using obabel
            print("  Step 1: Adding hydrogens...")
            step_start = time.time()
            if not add_hydrogens_obabel(sdf_file):
                print(f"  ⚠ Warning: Failed to add hydrogens for {base_name}")
                failed_count += 1
                continue
            
            # Check time after step 1
            elapsed = time.time() - start_time
            if elapsed >= MAX_PROCESSING_TIME_PER_LIGAND:
                print(f"  ⏱ Skipping {base_name}: Time limit reached after Step 1 ({elapsed:.1f}s)")
                skipped_count += 1
                continue
            
            # Step 2: Assign partial charges (Gasteiger) using obabel
            print("  Step 2: Assigning Gasteiger charges...")
            if not assign_charges_obabel(sdf_file):
                print(f"  ⚠ Warning: Failed to assign charges for {base_name}")
                failed_count += 1
                continue
            
            # Check time after step 2
            elapsed = time.time() - start_time
            if elapsed >= MAX_PROCESSING_TIME_PER_LIGAND:
                print(f"  ⏱ Skipping {base_name}: Time limit reached after Step 2 ({elapsed:.1f}s)")
                skipped_count += 1
                continue
            
            # Step 3: Generate 3D structure using RDKit
            print("  Step 3: Generating 3D structure...")
            remaining_time = MAX_PROCESSING_TIME_PER_LIGAND - (time.time() - start_time)
            if remaining_time <= 0:
                print(f"  ⏱ Skipping {base_name}: No time remaining for Step 3")
                skipped_count += 1
                continue
            
            if not generate_3d_structure_with_timeout(sdf_file, remaining_time):
                print(f"  ⚠ Warning: Failed to generate 3D structure for {base_name}")
                failed_count += 1
                continue
            
            # Check final time
            elapsed = time.time() - start_time
            if elapsed >= MAX_PROCESSING_TIME_PER_LIGAND:
                print(f"  ⏱ Skipping {base_name}: Time limit reached after Step 3 ({elapsed:.1f}s)")
                skipped_count += 1
                continue
            
            # Step 4: Energy minimization
            print("  Step 4: Energy minimization...")
            remaining_time = MAX_PROCESSING_TIME_PER_LIGAND - (time.time() - start_time)
            if remaining_time <= 0:
                print(f"  ⏱ Skipping {base_name}: No time remaining for Step 4")
                skipped_count += 1
                continue
            
            if not minimize_energy_rdkit(sdf_file, remaining_time):
                print(f"  ⚠ Warning: Failed to minimize energy for {base_name} (continuing anyway)")
                # Continue processing even if minimization fails
            
            elapsed = time.time() - start_time
            print(f"  ✓ Successfully prepared {base_name} (took {elapsed:.1f}s)")
            successful_count += 1
            
        except ProcessingTimeoutError as e:
            elapsed = time.time() - start_time
            print(f"  ⏱ Skipping {base_name}: {e} (elapsed: {elapsed:.1f}s)")
            skipped_count += 1
            continue
        except Exception as e:
            elapsed = time.time() - start_time
            print(f"  ✗ Error processing {base_name}: {e} (elapsed: {elapsed:.1f}s)")
            failed_count += 1
            continue
    
    print(f"\n✅ Ligand preparation complete.")
    print(f"   Successful: {successful_count}/{len(sdf_files)}")
    if failed_count > 0:
        print(f"   Failed: {failed_count}/{len(sdf_files)}")
    if skipped_count > 0:
        print(f"   Skipped (timeout): {skipped_count}/{len(sdf_files)}")


def add_hydrogens_obabel(sdf_file):
    """Add hydrogens to SDF file using obabel with pH consideration"""
    try:
        abs_sdf_file = os.path.abspath(sdf_file)
        # -h: add hydrogens
        # -p 7.4: consider protonation at physiological pH (7.4)
        # This ensures all atoms get hydrogens with proper protonation state
        result = subprocess.run(
            ["obabel", abs_sdf_file, "-O", abs_sdf_file, "-h", "-p", "7.4"],
            capture_output=True,
            text=True,
            timeout=MAX_STEP_TIME
        )
        return result.returncode == 0
    except subprocess.TimeoutExpired:
        print(f"    Error: Timeout after {MAX_STEP_TIME}s")
        return False
    except Exception as e:
        print(f"    Error: {e}")
        return False


def assign_charges_obabel(sdf_file):
    """Assign Gasteiger partial charges using obabel"""
    try:
        abs_sdf_file = os.path.abspath(sdf_file)
        result = subprocess.run(
            ["obabel", abs_sdf_file, "-O", abs_sdf_file, "--partialcharge", "gasteiger"],
            capture_output=True,
            text=True,
            timeout=MAX_STEP_TIME
        )
        return result.returncode == 0
    except subprocess.TimeoutExpired:
        print(f"    Error: Timeout after {MAX_STEP_TIME}s")
        return False
    except Exception as e:
        print(f"    Error: {e}")
        return False


def generate_3d_structure_with_timeout(sdf_file, timeout_seconds):
    """
    Generate 3D structure using RDKit with timeout.
    Uses multiple embedding attempts with different methods for better success rate.
    Note: Energy minimization is performed separately using RDKit UFF, so no optimization here.
    Processes all molecules in the SDF file.
    
    Args:
        sdf_file: Path to SDF file
        timeout_seconds: Maximum time allowed for this operation
    
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Use signal-based timeout (Unix/Linux/macOS only)
        if hasattr(signal, 'SIGALRM'):
            with timeout_context(min(timeout_seconds, MAX_STEP_TIME)):
                return generate_3d_structure(sdf_file, timeout_seconds)
        else:
            # Windows doesn't support SIGALRM, use time-based check
            return generate_3d_structure(sdf_file, timeout_seconds)
    except ProcessingTimeoutError:
        print(f"    Error: Timeout after {timeout_seconds}s")
        return False
    except Exception as e:
        print(f"    Error: {e}")
        return False


def generate_3d_structure(sdf_file, timeout_seconds=None):
    """
    Generate 3D structure using RDKit.
    Uses multiple embedding attempts with different methods for better success rate.
    Note: Energy minimization is performed separately using RDKit UFF, so no optimization here.
    Processes all molecules in the SDF file.
    
    Args:
        sdf_file: Path to SDF file
        timeout_seconds: Optional timeout in seconds (for Windows compatibility)
    """
    start_time = time.time()
    
    try:
        # Read all molecules from SDF
        supplier = Chem.SDMolSupplier(sdf_file)
        mols = [m for m in supplier if m is not None]
        
        if not mols:
            print(f"    Error: Could not read any molecules from {sdf_file}")
            return False
        
        processed_mols = []
        
        # Process each molecule
        # Limit the number of molecules to process to avoid timeout
        max_molecules = 10  # Limit to first 10 molecules to avoid timeout
        molecules_to_process = mols[:max_molecules] if len(mols) > max_molecules else mols
        
        if len(mols) > max_molecules:
            print(f"    Warning: Processing only first {max_molecules} of {len(mols)} molecules to avoid timeout")
        
        for mol_idx, mol in enumerate(molecules_to_process):
            # Check timeout for Windows (or when signal-based timeout is not available)
            if timeout_seconds is not None:
                elapsed = time.time() - start_time
                if elapsed >= timeout_seconds:
                    print(f"    Warning: Timeout reached ({elapsed:.1f}s), stopping processing")
                    break
            if mol is None:
                continue
            
            # Remove existing 3D coordinates if present
            mol = Chem.RemoveHs(mol)
            mol = Chem.AddHs(mol)
            
            # Try to generate 3D structure with multiple attempts
            success = False
            max_attempts = 3  # Reduced from 5 to 3 to save time
            
            # Method 1: Standard ETKDG embedding
            for attempt in range(max_attempts):
                try:
                    # Use ETKDG (Experimental-Torsion-angle preference with Distance Geometry) method
                    params = AllChem.ETKDGv3()
                    params.randomSeed = 42 + attempt  # Different seed for each attempt
                    conf_id = AllChem.EmbedMolecule(mol, params)
                    
                    if conf_id >= 0:
                        # No optimization here - will be done by RDKit UFF later
                        success = True
                        break
                except Exception as e:
                    if attempt < max_attempts - 1:
                        continue
                    else:
                        if mol_idx == 0:  # Only print warning for first molecule
                            print(f"    Warning: ETKDG embedding failed: {e}")
            
            # Method 2: If ETKDG fails, try basic distance geometry
            if not success:
                try:
                    conf_id = AllChem.EmbedMolecule(mol, useRandomCoords=True)
                    if conf_id >= 0:
                        # No optimization here - will be done by RDKit UFF later
                        success = True
                except Exception as e:
                    if mol_idx == 0:  # Only print warning for first molecule
                        print(f"    Warning: Basic embedding failed: {e}")
            
            # Method 3: If still fails, try with 2D coordinates first
            if not success:
                try:
                    # Generate 2D coordinates
                    AllChem.Compute2DCoords(mol)
                    # Then try to embed 3D
                    conf_id = AllChem.EmbedMolecule(mol, useRandomCoords=True)
                    if conf_id >= 0:
                        # No optimization here - will be done by RDKit UFF later
                        success = True
                except Exception as e:
                    if mol_idx == 0:  # Only print warning for first molecule
                        print(f"    Warning: 2D-to-3D conversion failed: {e}")
            
            if success:
                processed_mols.append(mol)
            else:
                if mol_idx == 0:  # Only print error for first molecule
                    print(f"    Warning: Could not generate 3D structure for molecule {mol_idx + 1}")
                # Still add the molecule even if 3D generation failed (it may have existing coords)
                processed_mols.append(mol)
        
        if not processed_mols:
            print(f"    Error: Could not process any molecules")
            return False
        
        # Write all processed molecules back to SDF file
        writer = Chem.SDWriter(sdf_file)
        for mol in processed_mols:
            writer.write(mol)
        writer.close()
        
        return True
        
    except Exception as e:
        print(f"    Error generating 3D structure: {e}")
        return False


def minimize_energy_rdkit(sdf_file, timeout_seconds=None):
    """
    Minimize energy using RDKit UFF (Universal Force Field).
    Simple energy minimization for all molecules in the SDF file.
    
    Args:
        sdf_file: Path to SDF file
        timeout_seconds: Optional timeout in seconds
    
    Returns:
        bool: True if successful, False otherwise
    """
    start_time = time.time()
    
    try:
        # Read all molecules from SDF
        supplier = Chem.SDMolSupplier(sdf_file)
        mols = [m for m in supplier if m is not None]
        
        if not mols:
            print(f"    Error: Could not read any molecules from {sdf_file}")
            return False
        
        processed_mols = []
        
        # Process each molecule
        for mol_idx, mol in enumerate(mols):
            # Check timeout
            if timeout_seconds is not None:
                elapsed = time.time() - start_time
                if elapsed >= timeout_seconds:
                    print(f"    Warning: Timeout reached ({elapsed:.1f}s), stopping minimization")
                    break
            
            if mol is None:
                continue
            
            try:
                # Add hydrogens
                mol = Chem.AddHs(mol)
                
                # Generate 3D structure
                AllChem.EmbedMolecule(mol, useRandomCoords=True)
                
                # Perform UFF energy minimization
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
                
                processed_mols.append(mol)
                    
            except Exception as e:
                if mol_idx == 0:  # Only print warning for first molecule
                    print(f"    Warning: Error processing molecule {mol_idx + 1}: {e}")
                processed_mols.append(mol)
        
        if not processed_mols:
            print(f"    Error: Could not process any molecules")
            return False
        
        # Write all processed molecules back to SDF file
        writer = Chem.SDWriter(sdf_file)
        for mol in processed_mols:
            writer.write(mol)
        writer.close()
        
        return True
        
    except Exception as e:
        print(f"    Error minimizing energy: {e}")
        return False


if __name__ == "__main__":
    main()

