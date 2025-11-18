#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Docking Results Report Generation Script
Converted from Report/report.ipynb
"""

import glob
import os
import re
import shutil
from pathlib import Path


def main():
    """Main execution function"""
    print("Starting docking results report generation...")

    # Check current directory
    print(f"Current working directory: {os.getcwd()}")

    # Analyze docking results and generate ranking
    generate_docking_ranking()

    # Copy top-ranked compound file
    copy_top_compound()

    # Copy ranking file to Visualization directory
    copy_ranking_file()

    print("Report generation complete.")


def parse_smina_log(log_file):
    """
    Extract mode 1 affinity value from a Smina log file.

    Args:
        log_file (str): Path to the log file

    Returns:
        float or None: Binding energy value (kcal/mol), or None if not found
    """
    with open(log_file, "r") as f:
        for line in f:
            if line.strip().startswith("1 "):  # Get the line for mode 1
                parts = line.split()
                if len(parts) > 1:
                    try:
                        affinity = float(parts[1])  # Affinity value (kcal/mol)
                        return affinity
                    except ValueError:
                        pass
    return None  # No mode 1 data found


def generate_docking_ranking():
    """Generate docking result ranking"""
    print("\n=== Generating docking result ranking ===")

    # Get docking log files
    log_files = glob.glob("docking_results/*_docking.log")

    # Extract affinity values from each log
    results = []
    for log_file in log_files:
        affinity = parse_smina_log(log_file)
        if affinity is not None:
            compound_name = log_file.split("/")[-1].replace("_docking.log", "")
            results.append((compound_name, affinity))

    # Sort by affinity (lower = stronger binding)
    results.sort(key=lambda x: x[1])

    # Output results to a text file
    output_file = "./results/docking_ranking.txt"

    with open(output_file, "w", encoding="utf-8") as f:
        f.write("Docking Result Ranking (Strongest Binding First)\n")
        f.write("----------------------------------------\n")
        for rank, (compound, affinity) in enumerate(results, 1):
            f.write(
                f"Rank {rank}: Compound {compound}, Binding Energy: {affinity:.2f} kcal/mol\n"
            )

    # Display summary
    print(f"Ranked docking results for {len(results)} compounds.")
    print(f"Results saved to {output_file}")

    # Display top 10 results
    print("\n=== Top 10 Compounds ===")
    for rank, (compound, affinity) in enumerate(results[:10], 1):
        print(f"Rank {rank}: Compound {compound}, Binding Energy: {affinity:.2f} kcal/mol")


def copy_top_compound():
    """Copy the top compound’s file"""
    print("\n=== Copying Top Compound File ===")

    # Read ranking file
    with open("./results/docking_ranking.txt", encoding="utf-8") as f:
        text = f.read()

    # Extract top compound name
    m = re.search(r"Rank\s*1:.*Compound\s+([^\s,，]+)", text)
    if not m:
        raise ValueError("Top compound name not found.")
    top_ligand = m.group(1)
    print(f"Top ligand: {top_ligand}")

    # Source SDF file path
    src = Path(f"docking_results/{top_ligand}_docked.sdf")

    # Destination directory (results)
    dst_dir = Path("./results")
    dst_dir.mkdir(exist_ok=True)

    # Destination path
    dst = dst_dir / src.name

    # Copy file
    shutil.copy(src, dst)

    print(f"Copied {src} → {dst}")


def copy_ranking_file():
    """Copy ranking file to Visualization directory"""
    print("\n=== Copying Ranking File ===")

    # Ranking file is already in ./results, so this is just a message
    print("Ranking file already saved in ./results/.")


if __name__ == "__main__":
    main()
