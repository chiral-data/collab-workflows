#!/usr/bin/env python3
"""
Node 03: Boltz-2 structure prediction.

Wraps the `boltz predict` CLI, collects output PDB files and
confidence metrics (pLDDT, pTM, ipTM, PAE, PDE) into a
boltz_summary.json for use by the visualization node.
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


# ── Prediction ────────────────────────────────────────────────────────────────

def run_boltz_predict(input_file, output_dir, diffusion_samples, recycling_steps, use_msa_server):
    """Run the boltz predict CLI command."""
    cmd = [
        "boltz", "predict", input_file,
        "--output_format", "pdb",
        "--diffusion_samples", str(diffusion_samples),
        "--recycling_steps",   str(recycling_steps),
        "--devices",           "1",
        "--accelerator",       "gpu",
        "--no_kernels",
    ]

    if str(use_msa_server).lower() == "true":
        cmd.append("--use_msa_server")

    print(f"\nRunning: {' '.join(cmd)}\n")
    result = subprocess.run(cmd, capture_output=False, text=True)

    if result.returncode != 0:
        print(f"Error: Boltz-2 failed with return code {result.returncode}")
        sys.exit(1)


# ── Output collection ─────────────────────────────────────────────────────────

def collect_outputs(input_file, output_dir):
    """Collect PDB and confidence JSON files from Boltz output directories."""
    input_basename = Path(input_file).stem

    # Boltz writes to boltz_results_<name>/predictions/<name>/
    search_patterns = [
        f"boltz_results_{input_basename}/predictions/{input_basename}/*",
        f"boltz_results_{input_basename}/*",
    ]

    collected = []
    for pattern in search_patterns:
        for filepath in glob.glob(pattern):
            if os.path.isfile(filepath):
                dest = os.path.join(output_dir, os.path.basename(filepath))
                shutil.copy2(filepath, dest)
                collected.append(os.path.basename(filepath))
                print(f"  Collected: {os.path.basename(filepath)}")

    # Fallback: walk any predictions/ subdirectory
    if not collected:
        for dirpath, _, filenames in os.walk('.'):
            if 'predictions' in dirpath:
                for f in filenames:
                    src = os.path.join(dirpath, f)
                    dest = os.path.join(output_dir, f)
                    shutil.copy2(src, dest)
                    collected.append(f)
                    print(f"  Collected (fallback): {f}")

    # Clean up intermediate Boltz result directories
    for d in glob.glob("boltz_results_*"):
        if os.path.isdir(d):
            shutil.rmtree(d)
            print(f"  Cleaned up: {d}")

    return collected


# ── Confidence parsing ────────────────────────────────────────────────────────

def parse_confidence_files(output_dir):
    """Extract pLDDT, pTM, ipTM, PAE, PDE from Boltz confidence JSONs."""
    confidence_files = glob.glob(os.path.join(output_dir, "confidence_*.json"))
    metrics = []

    for cf in sorted(confidence_files):
        with open(cf, 'r') as f:
            data = json.load(f)

        sample_name = Path(cf).stem  # e.g. confidence_boltz_input_model_0

        entry = {
            'sample':  sample_name,
            'plddt':   data.get('plddt',  None),
            'ptm':     data.get('ptm',    None),
            'iptm':    data.get('iptm',   None),
        }

        # PAE and PDE may be nested or top-level depending on Boltz version
        if 'pae' in data:
            pae = data['pae']
            # If it's a matrix, store the mean
            if isinstance(pae, list):
                flat = [v for row in pae for v in (row if isinstance(row, list) else [row])]
                entry['pae_mean'] = round(sum(flat) / len(flat), 4) if flat else None
            else:
                entry['pae_mean'] = pae

        if 'pde' in data:
            pde = data['pde']
            if isinstance(pde, list):
                flat = [v for row in pde for v in (row if isinstance(row, list) else [row])]
                entry['pde_mean'] = round(sum(flat) / len(flat), 4) if flat else None
            else:
                entry['pde_mean'] = pde

        metrics.append(entry)
        print(f"  {sample_name}: pLDDT={entry['plddt']}, pTM={entry['ptm']}, "
              f"ipTM={entry['iptm']}, PAE={entry.get('pae_mean')}, PDE={entry.get('pde_mean')}")

    return metrics


# ── Summary ───────────────────────────────────────────────────────────────────

def write_summary(output_dir, collected_files, metrics, params):
    """Write boltz_summary.json for the visualization node."""
    pdb_files = [f for f in collected_files if f.endswith('.pdb')]

    summary = {
        'tool':       'boltz2',
        'params':     params,
        'pdb_files':  pdb_files,
        'confidence': metrics,
    }

    summary_path = os.path.join(output_dir, 'boltz_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\n  Summary written to: boltz_summary.json")
    print(f"  PDB files: {len(pdb_files)}")
    print(f"  Confidence entries: {len(metrics)}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Run Boltz-2 structure prediction')
    parser.add_argument('--input',             required=True,       help='Input YAML file from node 02')
    parser.add_argument('--output-dir',        default='.',         help='Output directory')
    parser.add_argument('--diffusion-samples', type=int, default=2, help='Number of diffusion samples')
    parser.add_argument('--recycling-steps',   type=int, default=3, help='Number of recycling steps')
    parser.add_argument('--use-msa-server',    default='true',      help='Use MSA server (true/false)')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Input:             {args.input}")
    print(f"Diffusion samples: {args.diffusion_samples}")
    print(f"Recycling steps:   {args.recycling_steps}")
    print(f"MSA server:        {args.use_msa_server}")

    # 1. Run prediction
    run_boltz_predict(
        input_file=args.input,
        output_dir=args.output_dir,
        diffusion_samples=args.diffusion_samples,
        recycling_steps=args.recycling_steps,
        use_msa_server=args.use_msa_server,
    )

    # 2. Collect output files
    print("\nCollecting outputs:")
    collected = collect_outputs(args.input, args.output_dir)

    # 3. Parse confidence metrics
    print("\nParsing confidence metrics:")
    metrics = parse_confidence_files(args.output_dir)

    # Warn if expected outputs are missing
    pdb_files = [f for f in collected if f.endswith('.pdb')]
    if not pdb_files:
        print("Warning: No PDB files found in output")
    if not metrics:
        print("Warning: No confidence JSON files found")

    # 4. Write summary for node 05
    params = {
        'diffusion_samples': args.diffusion_samples,
        'recycling_steps':   args.recycling_steps,
        'use_msa_server':    args.use_msa_server,
    }
    write_summary(args.output_dir, collected, metrics, params)

    print("\nNode 03 completed ✓")


if __name__ == '__main__':
    main()
