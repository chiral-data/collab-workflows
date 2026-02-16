#!/usr/bin/env python3
"""
Run BoltzGen binder design.

Wraps the `boltzgen` CLI with workflow parameters,
then collects and organizes output files.
"""

import argparse
import glob
import os
import shutil
import subprocess
import sys
from pathlib import Path


def run_boltzgen_design(design_spec, output_dir, protocol, num_designs, budget):
    """Run BoltzGen design and collect results."""
    # Build command
    boltzgen_output = "boltzgen_output"
    cmd = [
        "boltzgen", "run", design_spec,
        "--protocol", protocol,
        "--num_designs", str(num_designs),
        "--budget", str(budget),
        "--devices", "1",
        "--output", boltzgen_output,
    ]

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False, text=True)

    if result.returncode != 0:
        print(f"Error: BoltzGen design failed with return code {result.returncode}")
        sys.exit(1)

    # Collect results from BoltzGen output directory
    input_basename = Path(design_spec).stem
    results_patterns = [
        f"{boltzgen_output}/**/*",
        f"{boltzgen_output}/*",
        f"boltzgen_results_{input_basename}/**/*",
        f"boltzgen_results_{input_basename}/*",
        f"results_{input_basename}/**/*",
        f"results_{input_basename}/*",
    ]

    collected = 0
    for pattern in results_patterns:
        for filepath in glob.glob(pattern, recursive=True):
            if os.path.isfile(filepath):
                dest = os.path.join(output_dir, os.path.basename(filepath))
                # Avoid overwriting files with same name from different subdirs
                if os.path.exists(dest):
                    base, ext = os.path.splitext(os.path.basename(filepath))
                    parent = os.path.basename(os.path.dirname(filepath))
                    dest = os.path.join(output_dir, f"{parent}_{base}{ext}")
                shutil.copy2(filepath, dest)
                print(f"  Collected: {os.path.basename(dest)}")
                collected += 1

    if collected == 0:
        # Try walking any results directory
        for dirpath, dirnames, filenames in os.walk('.'):
            if 'results' in dirpath.lower() or 'output' in dirpath.lower():
                for f in filenames:
                    src = os.path.join(dirpath, f)
                    dest = os.path.join(output_dir, f)
                    if not os.path.exists(dest):
                        shutil.copy2(src, dest)
                        print(f"  Collected: {f}")
                        collected += 1

    print(f"Collected {collected} output files")

    # Report on outputs
    cif_files = glob.glob(os.path.join(output_dir, "*.cif"))
    csv_files = glob.glob(os.path.join(output_dir, "*.csv"))
    print(f"CIF files: {len(cif_files)}, CSV files: {len(csv_files)}")

    if not cif_files:
        print("Warning: No CIF design files found in output")

    # Clean up intermediate directories
    for d in glob.glob("boltzgen_results_*") + glob.glob("results_*"):
        if os.path.isdir(d):
            shutil.rmtree(d)
            print(f"  Cleaned up: {d}")


def main():
    parser = argparse.ArgumentParser(description='Run BoltzGen binder design')
    parser.add_argument('--design-spec', required=True, help='Design specification YAML')
    parser.add_argument('--output-dir', default='.', help='Output directory')
    parser.add_argument('--protocol', default='protein-anything', help='Design protocol')
    parser.add_argument('--num-designs', type=int, default=50, help='Number of designs')
    parser.add_argument('--budget', type=int, default=10, help='Final designs after filtering')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    run_boltzgen_design(
        design_spec=args.design_spec,
        output_dir=args.output_dir,
        protocol=args.protocol,
        num_designs=args.num_designs,
        budget=args.budget,
    )


if __name__ == '__main__':
    main()
