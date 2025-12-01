#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reporting - Generate report
Input: Docking score (docking result files from smina_screening)
Output: Score report (score report)
"""

import glob
import json
import os
import re
import shutil
from pathlib import Path

DEFAULT_PDB_ID = "5Y7J"
DEFAULT_CHAIN_ID = "A"

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
                return params.get("pdb_id", DEFAULT_PDB_ID)
        except Exception as e:
            print(f"⚠ Warning: Could not load global_params.json: {e}")
    return DEFAULT_PDB_ID

# Get script directory and set paths relative to script location
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
# Input from other nodes should be in input/ directory
# Output from this node should be in outputs/results/ directory (as specified in job.toml)
INPUT_DIR = os.path.join(SCRIPT_DIR, "input")
# Base output directory
OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")
# Actual results directory
RESULTS_DIR = os.path.join(OUTPUT_DIR, "results")

# docking_results from docking/outputs/ (files are directly in outputs/, not in a subdirectory)
DOCKING_OUTPUT_DIR = os.path.join(SCRIPT_DIR, "..", "docking", "outputs")
DOCKING_RESULTS_DIR = DOCKING_OUTPUT_DIR
# Fallback to input/ for backward compatibility
if not os.path.exists(DOCKING_RESULTS_DIR):
    DOCKING_RESULTS_DIR = INPUT_DIR
# PDB files from protein_extraction should be in input/ (optional)
PREPARATION_DIR = INPUT_DIR
# Output should be in "outputs/results" directory
RANKING_FILE = os.path.join(RESULTS_DIR, "docking_ranking.txt")

def main():
    """Main execution function"""
    global OUTPUT_DIR, RESULTS_DIR, RANKING_FILE
    
    print("=== Generate report ===")
    
    # Reset directories to ensure correctness
    OUTPUT_DIR = os.path.join(SCRIPT_DIR, "outputs")
    RESULTS_DIR = os.path.join(OUTPUT_DIR, "results")
    RANKING_FILE = os.path.join(RESULTS_DIR, "docking_ranking.txt")
    
    # Find docking results from docking/outputs/ (files are directly in outputs/, not in a subdirectory)
    if not os.path.exists(DOCKING_RESULTS_DIR):
        print(f"❌ Error: docking outputs directory not found.")
        print(f"   Searched in: {DOCKING_OUTPUT_DIR}")
        print(f"   Searched in: {INPUT_DIR}")
        if os.path.exists(os.path.join(SCRIPT_DIR, "..", "docking", "outputs")):
            files = [f for f in os.listdir(os.path.join(SCRIPT_DIR, "..", "docking", "outputs")) if f.endswith(('.sdf', '.txt'))]
            if files:
                print(f"   Available files in docking/outputs: {files[:10]}...")
        if os.path.exists(INPUT_DIR):
            dirs = [d for d in os.listdir(INPUT_DIR) if os.path.isdir(os.path.join(INPUT_DIR, d))]
            if dirs:
                print(f"   Available directories in input: {dirs}")
        print("   Please ensure docking has completed.")
        exit(1)
    
    print(f"✓ Found docking results in: {DOCKING_RESULTS_DIR}")
    
    # Find PDB files in input/ directory (copied by run.sh)
    # PREPARATION_DIR is already set to INPUT_DIR at module level
    if not os.path.exists(PREPARATION_DIR):
        print(f"⚠ Warning: PDB files directory not found in input/.")
    
    # Ensure OUTPUT_DIR is not nested (double-check and fix if needed)
    OUTPUT_DIR = os.path.abspath(OUTPUT_DIR)
    RANKING_FILE = os.path.join(OUTPUT_DIR, "docking_ranking.txt")
    
    # Remove nested outputs directory if it exists
    nested_outputs = os.path.join(OUTPUT_DIR, "outputs")
    if os.path.exists(nested_outputs) and os.path.isdir(nested_outputs):
        print(f"⚠ Warning: Found nested outputs directory: {nested_outputs}")
        print(f"   Removing it to avoid confusion...")
        shutil.rmtree(nested_outputs)
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Verify that OUTPUT_DIR was created
    if not os.path.exists(OUTPUT_DIR):
        print(f"❌ Error: Failed to create output directory: {OUTPUT_DIR}")
        exit(1)
    
    print(f"Input directory (docking_results): {DOCKING_RESULTS_DIR}")
    print(f"Preparation directory (PDB files): {PREPARATION_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"✓ Output directory created: {OUTPUT_DIR}")
    
    # Generate docking result ranking
    try:
        results = generate_docking_ranking(DOCKING_RESULTS_DIR)
    except Exception as e:
        print(f"❌ Error in generate_docking_ranking: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
    
    # Copy top 3 compound files and generate complexes
    try:
        copy_top_compounds(results, DOCKING_RESULTS_DIR, PREPARATION_DIR)
    except Exception as e:
        print(f"❌ Error in copy_top_compounds: {e}")
        import traceback
        traceback.print_exc()
        exit(1)
    
    # Verify that results directory still exists and has files
    if os.path.exists(RESULTS_DIR):
        result_files = os.listdir(RESULTS_DIR)
        print(f"✅ Report generation complete. Results directory contains {len(result_files)} file(s).")
    else:
        print(f"⚠ Warning: Results directory was removed or does not exist: {RESULTS_DIR}")


def parse_smina_log(log_file):
    """
    Extract mode1 affinity value from Smina log file.
    
    Args:
        log_file (str): Path to log file
    
    Returns:
        float or None: Binding energy value (kcal/mol) or None
    """
    try:
        with open(log_file, "r") as f:
            for line in f:
                if line.strip().startswith("1 "):  # Get mode1 line
                    parts = line.split()
                    if len(parts) > 1:
                        try:
                            affinity = float(parts[1])  # Affinity value (kcal/mol)
                            return affinity
                        except ValueError:
                            pass
    except Exception as e:
        print(f"Warning: Error occurred while reading {log_file}: {e}")
    return None  # mode1 data not found


def generate_docking_ranking(input_dir):
    """Generate docking result ranking"""
    print("\n--- Generate docking result ranking ---")
    
    # Get docking log files
    log_files = glob.glob(os.path.join(input_dir, "*_log.txt"))
    
    if not log_files:
        error_msg = f"No log files found in {input_dir}."
        print(f"❌ Error: {error_msg}")
        print(f"   Searched pattern: {os.path.join(input_dir, '*_log.txt')}")
        if os.path.exists(input_dir):
            available_files = os.listdir(input_dir)
            print(f"   Available files: {available_files[:20]}...")  # Show first 20 files
        raise FileNotFoundError(error_msg)
    
    # Extract affinity from each log
    results = []
    true_ligand_result = None
    for log_file in log_files:
        affinity = parse_smina_log(log_file)
        if affinity is not None:
            compound_name = os.path.basename(log_file).replace("_log.txt", "")
            results.append((compound_name, affinity))
            # Check if this is true_ligand
            if compound_name == "true_ligand":
                true_ligand_result = (compound_name, affinity)
    
    if not results:
        error_msg = "No valid docking results found."
        print(f"❌ Error: {error_msg}")
        raise ValueError(error_msg)
    
    # Sort by affinity (binding energy) in ascending order (lower = stronger binding)
    results.sort(key=lambda x: x[1])
    
    # Find true_ligand's rank in the sorted results
    true_ligand_rank = None
    if true_ligand_result:
        for rank, (compound, affinity) in enumerate(results, 1):
            if compound == "true_ligand":
                true_ligand_rank = rank
                break
    
    # Write results to text file
    with open(RANKING_FILE, "w", encoding="utf-8") as f:
        f.write("Docking Result Ranking (sorted by binding strength)\n")
        f.write("=" * 50 + "\n")
        for rank, (compound, affinity) in enumerate(results, 1):
            marker = " [TRUE LIGAND]" if compound == "true_ligand" else ""
            f.write(
                f"Rank {rank}: Compound {compound}, Binding energy: {affinity:.2f} kcal/mol{marker}\n"
            )
        # Add summary for true_ligand if it exists
        if true_ligand_rank:
            f.write("\n" + "=" * 50 + "\n")
            f.write(f"True Ligand (true_ligand) Ranking: {true_ligand_rank} / {len(results)}\n")
            f.write(f"True Ligand Binding Energy: {true_ligand_result[1]:.2f} kcal/mol\n")
    
    # Display results
    print(f"Ranked docking results for {len(results)} compounds")
    if true_ligand_rank:
        print(f"✓ True ligand found: Rank {true_ligand_rank} / {len(results)} (Binding energy: {true_ligand_result[1]:.2f} kcal/mol)")
    print(f"Results saved to {RANKING_FILE}")
    
    # Display top 10 results
    print("\n--- Top 10 compounds ---")
    for rank, (compound, affinity) in enumerate(results[:10], 1):
        marker = " [TRUE LIGAND]" if compound == "true_ligand" else ""
        print(f"Rank {rank}: Compound {compound}, Binding energy: {affinity:.2f} kcal/mol{marker}")
    
    return results


def copy_top_compounds(results, input_dir, preparation_dir):
    """Copy top 3 compound files and receptor PDB file separately"""
    print("\n--- Copy top 3 compound files and receptor ---")
    
    if not results:
        print("Warning: No results available for copying.")
        return
    
    # Load parameters from global_params.json or environment variables
    global_pdb_id = load_global_params()
    pdb_id = os.environ.get("PDB_ID") or os.environ.get("PARAM_PDB_ID") or global_pdb_id
    chain_id = os.environ.get("CHAIN_ID", DEFAULT_CHAIN_ID)
    
    # Find receptor file with multiple fallback options
    receptor_file = None
    
    # Priority 1: Chain-specific cleaned PDB (e.g., 4OHU_chain_A_clean.pdb)
    if chain_id:
        receptor_file = os.path.join(preparation_dir, f"{pdb_id}_chain_{chain_id}_clean.pdb")
        if os.path.exists(receptor_file):
            print(f"Found receptor file: {os.path.basename(receptor_file)}")
    
    # Priority 2: General cleaned PDB (e.g., 4OHU_clean.pdb)
    if not receptor_file or not os.path.exists(receptor_file):
        receptor_file = os.path.join(preparation_dir, f"{pdb_id}_clean.pdb")
        if os.path.exists(receptor_file):
            print(f"Found receptor file: {os.path.basename(receptor_file)}")
    
    # Priority 3: Search for any PDB file with _clean in the name
    if not receptor_file or not os.path.exists(receptor_file):
        clean_files = glob.glob(os.path.join(preparation_dir, f"{pdb_id}*_clean.pdb"))
        if clean_files:
            receptor_file = clean_files[0]
            print(f"Found receptor file: {os.path.basename(receptor_file)}")
    
    # Priority 4: Search for any PDB file starting with pdb_id
    if not receptor_file or not os.path.exists(receptor_file):
        pdb_files = glob.glob(os.path.join(preparation_dir, f"{pdb_id}*.pdb"))
        if pdb_files:
            # Prefer files with "clean" in the name
            clean_files = [f for f in pdb_files if "clean" in os.path.basename(f).lower()]
            if clean_files:
                receptor_file = clean_files[0]
            else:
                receptor_file = pdb_files[0]
            print(f"Found receptor file: {os.path.basename(receptor_file)}")
    
    # Get top 3 compounds
    top_n = min(3, len(results))
    top_compounds = results[:top_n]
    
    # Destination directory
    dst_dir = Path(RESULTS_DIR)
    dst_dir.mkdir(exist_ok=True)
    
    # Copy receptor PDB file (used for docking)
    if receptor_file and os.path.exists(receptor_file):
        receptor_dst = dst_dir / f"receptor_{os.path.basename(receptor_file)}"
        shutil.copy(receptor_file, receptor_dst)
        print(f"✓ Copied receptor PDB: {receptor_dst.name}")
    else:
        print(f"⚠ Warning: Receptor file not found in {preparation_dir}")
        print(f"   Searched for: {pdb_id}_chain_{chain_id}_clean.pdb, {pdb_id}_clean.pdb, {pdb_id}*_clean.pdb")
        if os.path.exists(preparation_dir):
            available_files = [f for f in os.listdir(preparation_dir) if f.endswith('.pdb')]
            print(f"   Available PDB files: {available_files[:10]}...")
    
    copied_count = 0
    
    # Track which compounds have been copied
    copied_compounds = set()
    
    # Copy top 3 compounds
    for rank, (ligand_name, affinity) in enumerate(top_compounds, 1):
        print(f"\nRank {rank}: {ligand_name} (Binding energy: {affinity:.2f} kcal/mol)")
        
        # Copy docked ligand SDF file
        src_sdf = Path(input_dir) / f"{ligand_name}_docked.sdf"
        if src_sdf.exists():
            dst_sdf = dst_dir / f"top{rank}_{ligand_name}_docked.sdf"
            shutil.copy(src_sdf, dst_sdf)
            print(f"  ✓ Copied docked ligand: {dst_sdf.name}")
            copied_count += 1
            copied_compounds.add(ligand_name)
        else:
            print(f"  ⚠ Warning: Docked SDF file not found: {src_sdf.name}")
    
    # Always copy true_ligand if it exists and hasn't been copied yet
    true_ligand_name = "true_ligand"
    true_ligand_src = Path(input_dir) / f"{true_ligand_name}_docked.sdf"
    if true_ligand_src.exists() and true_ligand_name not in copied_compounds:
        # Find true_ligand's rank in the full results
        true_ligand_rank = None
        for rank, (compound, affinity) in enumerate(results, 1):
            if compound == true_ligand_name:
                true_ligand_rank = rank
                break
        
        if true_ligand_rank:
            print(f"\n--- True Ligand (Rank {true_ligand_rank}) ---")
            true_ligand_affinity = next((aff for name, aff in results if name == true_ligand_name), None)
            if true_ligand_affinity:
                print(f"Rank {true_ligand_rank}: {true_ligand_name} (Binding energy: {true_ligand_affinity:.2f} kcal/mol)")
            
            dst_sdf = dst_dir / f"true_ligand_rank{true_ligand_rank}_docked.sdf"
            shutil.copy(true_ligand_src, dst_sdf)
            print(f"  ✓ Copied true ligand: {dst_sdf.name}")
            copied_count += 1
            copied_compounds.add(true_ligand_name)
    
    print(f"\n✅ Copied receptor PDB and {copied_count} ligand file(s) to {RESULTS_DIR}/")
    if true_ligand_name in copied_compounds:
        print(f"   (Including true_ligand)")




if __name__ == "__main__":
    main()

