#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Docking Results Report Script
Converted from Report/report.ipynb
"""

import glob
import os
import re
import shutil
from pathlib import Path

# Read parameters from environment variables
top_n = int(os.environ.get("PARAM_TOP_N", "10"))


def main():
    """Main execution function"""
    print("Starting docking results report generation...")
    print(f"Top N compounds to display: {top_n}")

    # Check current directory
    print(f"Current directory: {os.getcwd()}")

    # Generate docking ranking
    results = generate_docking_ranking()

    # Copy top compound file (only if we have results)
    if results:
        copy_top_compound(results[0][0])
    else:
        print("\nNo docking results found. Skipping top compound copy.")

    # Copy ranking file
    copy_ranking_file()

    print("Report generation completed.")


def parse_smina_log(log_file):
    """
    Extract mode1 affinity value from Smina log file.

    Args:
        log_file (str): Path to log file

    Returns:
        float or None: Binding energy (kcal/mol) or None
    """
    with open(log_file, "r") as f:
        for line in f:
            if line.strip().startswith("1 "):
                parts = line.split()
                if len(parts) > 1:
                    try:
                        affinity = float(parts[1])
                        return affinity
                    except ValueError:
                        pass
    return None


def generate_docking_ranking():
    """Generate docking results ranking"""
    print("\n=== Generating docking results ranking ===")

    # Get docking log files (silva copies them to ./inputs/ subdirectory)
    log_files = glob.glob("./inputs/*_log.txt")

    # Extract affinity from each log
    results = []
    for log_file in log_files:
        affinity = parse_smina_log(log_file)
        if affinity is not None:
            compound_name = os.path.basename(log_file).replace("_log.txt", "")
            results.append((compound_name, affinity))

    # Sort by affinity (lower = stronger binding)
    results.sort(key=lambda x: x[1])

    # Write results to text file
    os.makedirs("./results", exist_ok=True)
    output_file = "./results/docking_ranking.txt"

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("Docking Results Ranking (sorted by binding strength)\n")
        f.write("----------------------------------------\n")
        if results:
            for rank, (compound, affinity) in enumerate(results, 1):
                f.write(
                    f"Rank {rank}: Compound {compound}, Binding energy: {affinity:.2f} kcal/mol\n"
                )
        else:
            f.write("No docking results found.\n")

    print(f"Ranked {len(results)} compounds")
    print(f"Results saved to {output_file}")

    # Display top N results
    print(f"\n=== Top {top_n} compounds ===")
    if results:
        for rank, (compound, affinity) in enumerate(results[:top_n], 1):
            print(f"Rank {rank}: Compound {compound}, Binding energy: {affinity:.2f} kcal/mol")
    else:
        print("No results to display.")

    return results


def copy_top_compound(top_ligand):
    """Copy top compound file"""
    print("\n=== Copying top compound file ===")
    print(f"Top ligand: {top_ligand}")

    # Source SDF file path (silva copies files to ./inputs/ subdirectory)
    src = Path(f"./inputs/{top_ligand}_docked.sdf")

    # Destination directory
    dst_dir = Path("./results")
    dst_dir.mkdir(exist_ok=True)

    # Destination path
    dst = dst_dir / src.name

    # Copy file
    if src.exists():
        shutil.copy(src, dst)
        print(f"Copied {src} to {dst}.")
    else:
        print(f"Warning: Source file {src} not found.")


def copy_ranking_file():
    """Copy ranking file to results directory"""
    print("\n=== Copying ranking file ===")
    print("Ranking file is already in ./results/.")


if __name__ == "__main__":
    main()
