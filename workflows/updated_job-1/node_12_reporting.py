#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Node 12: Reporting - Generate report
Input: Docking score (docking result files)
Output: Score report (score report)
"""

import glob
import os
import re
import shutil
from pathlib import Path

INPUT_DIR = "docking_results"
OUTPUT_DIR = "results"
RANKING_FILE = f"{OUTPUT_DIR}/docking_ranking.txt"

def main():
    """Main execution function"""
    print("=== Node 12: Generate report ===")
    
    if not os.path.exists(INPUT_DIR):
        print(f"❌ Error: {INPUT_DIR} directory not found.")
        print("   Please run Node 11 (node_11_smina_screening.py) first.")
        exit(1)
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    print(f"Input directory: {INPUT_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    
    # Generate docking result ranking
    generate_docking_ranking()
    
    # Copy top compound file
    copy_top_compound()
    
    print("✅ Report generation complete.")


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


def generate_docking_ranking():
    """Generate docking result ranking"""
    print("\n--- Generate docking result ranking ---")
    
    # Get docking log files
    log_files = glob.glob(f"{INPUT_DIR}/*_log.txt")
    
    if not log_files:
        print(f"❌ Error: No log files found in {INPUT_DIR}.")
        exit(1)
    
    # Extract affinity from each log
    results = []
    for log_file in log_files:
        affinity = parse_smina_log(log_file)
        if affinity is not None:
            compound_name = os.path.basename(log_file).replace("_log.txt", "")
            results.append((compound_name, affinity))
    
    if not results:
        print("❌ Error: No valid docking results found.")
        exit(1)
    
    # Sort by affinity (binding energy) in ascending order (lower = stronger binding)
    results.sort(key=lambda x: x[1])
    
    # Write results to text file
    with open(RANKING_FILE, "w", encoding="utf-8") as f:
        f.write("Docking Result Ranking (sorted by binding strength)\n")
        f.write("=" * 50 + "\n")
        for rank, (compound, affinity) in enumerate(results, 1):
            f.write(
                f"Rank {rank}: Compound {compound}, Binding energy: {affinity:.2f} kcal/mol\n"
            )
    
    # Display results
    print(f"Ranked docking results for {len(results)} compounds")
    print(f"Results saved to {RANKING_FILE}")
    
    # Display top 10 results
    print("\n--- Top 10 compounds ---")
    for rank, (compound, affinity) in enumerate(results[:10], 1):
        print(f"Rank {rank}: Compound {compound}, Binding energy: {affinity:.2f} kcal/mol")


def copy_top_compound():
    """Copy top compound file"""
    print("\n--- Copy top compound file ---")
    
    if not os.path.exists(RANKING_FILE):
        print(f"❌ Error: {RANKING_FILE} not found.")
        return
    
    # Read ranking file
    with open(RANKING_FILE, encoding="utf-8") as f:
        text = f.read()
    
    # Extract top compound name
    m = re.search(r"Rank 1:.*Compound\s+([^\s,，]+)", text)
    if not m:
        print("Warning: Top compound name not found.")
        return
    
    top_ligand = m.group(1)  # e.g., "ligand_1" or "clean_drug1081"
    print(f"Top ligand: {top_ligand}")
    
    # Source SDF file path (try multiple patterns)
    src = Path(f"{INPUT_DIR}/{top_ligand}_docked.sdf")
    if not src.exists():
        # Original code filename pattern (e.g., clean_drug1081_docked.sdf)
        # If _docked is already in filename
        src = Path(f"{INPUT_DIR}/{top_ligand}_docked.sdf")
    
    if not src.exists():
        print(f"Warning: {src} not found.")
        return
    
    # Destination directory (results)
    dst_dir = Path(OUTPUT_DIR)
    dst_dir.mkdir(exist_ok=True)
    
    # Destination path
    dst = dst_dir / src.name
    
    # Copy file
    shutil.copy(src, dst)
    print(f"Copied {src} to {dst}.")


if __name__ == "__main__":
    main()

