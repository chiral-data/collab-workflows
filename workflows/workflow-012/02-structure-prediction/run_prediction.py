#!/usr/bin/env python3
"""
Run Boltz-2 protein structure prediction.

Wraps the `boltz predict` CLI command with parameters from the workflow,
then collects and organizes output files.
"""

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path


def run_boltz_predict(input_file, output_dir, diffusion_samples, recycling_steps, use_msa_server):
    """Run Boltz-2 prediction and collect results."""
    # Build command
    cmd = [
        "boltz", "predict", input_file,
        "--output_format", "pdb",
        "--diffusion_samples", str(diffusion_samples),
        "--recycling_steps", str(recycling_steps),
        "--devices", "1",
        "--accelerator", "gpu",
        "--no_kernels",
    ]

    if use_msa_server.lower() == "true":
        cmd.append("--use_msa_server")

    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=False, text=True)

    if result.returncode != 0:
        print(f"Error: Boltz-2 prediction failed with return code {result.returncode}")
        sys.exit(1)

    # Collect results from Boltz output directory
    # Boltz writes to boltz_results_<name>/predictions/<name>/
    input_basename = Path(input_file).stem
    results_patterns = [
        f"boltz_results_{input_basename}/predictions/{input_basename}/*",
        f"boltz_results_{input_basename}/*",
    ]

    collected = 0
    for pattern in results_patterns:
        for filepath in glob.glob(pattern):
            if os.path.isfile(filepath):
                dest = os.path.join(output_dir, os.path.basename(filepath))
                shutil.copy2(filepath, dest)
                print(f"  Collected: {os.path.basename(filepath)}")
                collected += 1

    if collected == 0:
        # Try alternative directory structure
        for dirpath, dirnames, filenames in os.walk('.'):
            if 'predictions' in dirpath:
                for f in filenames:
                    src = os.path.join(dirpath, f)
                    dest = os.path.join(output_dir, f)
                    shutil.copy2(src, dest)
                    print(f"  Collected: {f}")
                    collected += 1

    print(f"Collected {collected} output files")

    # Verify expected outputs
    pdb_files = glob.glob(os.path.join(output_dir, "*.pdb"))
    json_files = glob.glob(os.path.join(output_dir, "confidence_*.json"))
    print(f"PDB files: {len(pdb_files)}, Confidence files: {len(json_files)}")

    if not pdb_files:
        print("Warning: No PDB structure files found in output")
    if not json_files:
        print("Warning: No confidence JSON files found in output")

    # Clean up intermediate directories
    for d in glob.glob("boltz_results_*"):
        if os.path.isdir(d):
            shutil.rmtree(d)
            print(f"  Cleaned up: {d}")


def main():
    parser = argparse.ArgumentParser(description='Run Boltz-2 structure prediction')
    parser.add_argument('--input', required=True, help='Input YAML file')
    parser.add_argument('--output-dir', default='.', help='Output directory')
    parser.add_argument('--diffusion-samples', type=int, default=10, help='Number of diffusion samples')
    parser.add_argument('--recycling-steps', type=int, default=5, help='Number of recycling steps')
    parser.add_argument('--use-msa-server', default='true', help='Use MSA server (true/false)')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    run_boltz_predict(
        input_file=args.input,
        output_dir=args.output_dir,
        diffusion_samples=args.diffusion_samples,
        recycling_steps=args.recycling_steps,
        use_msa_server=args.use_msa_server,
    )


if __name__ == '__main__':
    main()
