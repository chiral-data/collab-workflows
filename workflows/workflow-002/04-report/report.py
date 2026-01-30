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
    generate_docking_ranking()

    # Copy top compound file
    copy_top_compound()

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

    # Get docking log files
    log_files = glob.glob("docking_results/*_log.txt")

    # Extract affinity from each log
    results = []
    for log_file in log_files:
        affinity = parse_smina_log(log_file)
        if affinity is not None:
            compound_name = log_file.split("/")[-1].replace("_log.txt", "")
            results.append((compound_name, affinity))

    # Sort by affinity (lower = stronger binding)
    results.sort(key=lambda x: x[1])

    # Write results to text file
    os.makedirs("./results", exist_ok=True)
    output_file = "./results/docking_ranking.txt"

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("Docking Results Ranking (sorted by binding strength)\n")
        f.write("----------------------------------------\n")
        for rank, (compound, affinity) in enumerate(results, 1):
            f.write(
                f"Rank {rank}: Compound {compound}, Binding energy: {affinity:.2f} kcal/mol\n"
            )

    print(f"Ranked {len(results)} compounds")
    print(f"Results saved to {output_file}")

    # Display top N results
    print(f"\n=== Top {top_n} compounds ===")
    for rank, (compound, affinity) in enumerate(results[:top_n], 1):
        print(f"Rank {rank}: Compound {compound}, Binding energy: {affinity:.2f} kcal/mol")


def copy_top_compound():
    """Copy top compound file"""
    print("\n=== Copying top compound file ===")

    # Read ranking file
    with open("./results/docking_ranking.txt", encoding="utf-8") as f:
        text = f.read()

    # Extract top compound name
    m = re.search(r"Rank 1:.*Compound\s+([^\s,]+)", text)
    if not m:
        raise ValueError("Top compound name not found")
    top_ligand = m.group(1)
    print(f"Top ligand: {top_ligand}")

    # Source SDF file path
    src = Path(f"docking_results/{top_ligand}_docked.sdf")

    # Destination directory
    dst_dir = Path("./results")
    dst_dir.mkdir(exist_ok=True)

    # Destination path
    dst = dst_dir / src.name

    # Copy file
    shutil.copy(src, dst)

    print(f"Copied {src} to {dst}.")


def copy_ranking_file():
    """Copy ranking file to results directory"""
    print("\n=== Copying ranking file ===")
    print("Ranking file is already in ./results/.")


if __name__ == "__main__":
    main()
