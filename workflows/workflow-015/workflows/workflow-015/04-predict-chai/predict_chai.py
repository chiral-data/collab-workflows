#!/usr/bin/env python3
"""
Node 04: Chai-1 structure prediction.

Wraps the `chai-lab fold` CLI, collects output CIF files and
confidence metrics (pLDDT, pTM, ipTM) into chai_summary.json
for use by the visualization node.

Chai-1 CLI signature:
  chai-lab fold <fasta_file> <output_dir> [OPTIONS]

Chai-1 output files:
  pred.model_idx_N.cif        — structure (CIF format); B-factor column = per-atom pLDDT (0-100)
  scores.model_idx_N.npz     — NumPy archive with pTM, ipTM, aggregate_score
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np


# ── Prediction ────────────────────────────────────────────────────────────────

def run_chai_fold(input_file, output_dir, num_trunk_recycles, num_diffusion_timesteps):
    """Run the chai-lab fold CLI command."""
    cmd = [
        "chai-lab", "fold", input_file, output_dir,
        "--num-trunk-recycles", str(num_trunk_recycles),
        "--num-diffn-timesteps", str(num_diffusion_timesteps),
    ]

    print(f"\nRunning: {' '.join(cmd)}\n")
    result = subprocess.run(cmd, capture_output=False, text=True)

    if result.returncode != 0:
        print(f"Error: Chai-1 failed with return code {result.returncode}")
        sys.exit(1)


# ── Output collection ─────────────────────────────────────────────────────────

def collect_outputs(output_dir):
    """Collect CIF and scores NPZ files from Chai output directory."""
    collected = []
    for pattern in ["pred.*.cif", "scores.*.npz"]:
        for filepath in glob.glob(os.path.join(output_dir, pattern)):
            collected.append(os.path.basename(filepath))
            print(f"  Found: {os.path.basename(filepath)}")
    return collected


# ── Confidence parsing ────────────────────────────────────────────────────────

def parse_plddt_from_cif(cif_path):
    """Extract mean pLDDT from B-factor column in Chai CIF file (0-100 scale).

    Chai writes per-atom pLDDT as _atom_site.B_iso_or_equiv (column index 17).
    """
    b_factors = []
    with open(cif_path) as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                parts = line.split()
                if len(parts) > 17:
                    try:
                        b_factors.append(float(parts[17]))
                    except ValueError:
                        pass
    return round(sum(b_factors) / len(b_factors), 4) if b_factors else None


def parse_confidence_files(output_dir):
    """Extract pLDDT, pTM, ipTM, aggregate_score from Chai output files.

    pLDDT: averaged from B-factors in pred.model_idx_N.cif (0-100 scale)
    pTM, ipTM, aggregate_score: from scores.model_idx_N.npz
    """
    score_files = sorted(glob.glob(os.path.join(output_dir, "scores.*.npz")))
    metrics = []

    for sf in score_files:
        data = np.load(sf, allow_pickle=False)

        # scores.model_idx_0.npz → model_0, pred.model_idx_0.cif
        stem = Path(sf).stem
        parts = stem.split('.')
        idx_part = parts[-1]  # e.g. model_idx_0
        model_label = idx_part.replace('model_idx_', 'model_')
        cif_path = os.path.join(output_dir, f"pred.{idx_part}.cif")

        plddt = parse_plddt_from_cif(cif_path) if os.path.exists(cif_path) else None

        def scalar(key):
            return round(float(data[key].item()), 4) if key in data else None

        entry = {
            'sample': model_label,
            'plddt': plddt,
            'ptm': scalar('ptm'),
            'iptm': scalar('iptm'),
            'aggregate_score': scalar('aggregate_score'),
        }
        metrics.append(entry)
        print(f"  {model_label}: pLDDT={entry['plddt']}, pTM={entry['ptm']}, ipTM={entry['iptm']}, agg={entry['aggregate_score']}")

    return metrics


# ── Summary ───────────────────────────────────────────────────────────────────

def write_summary(output_dir, collected_files, metrics, params):
    """Write chai_summary.json for the visualization node."""
    cif_files = [f for f in collected_files if f.endswith('.cif')]

    summary = {
        'tool': 'chai1',
        'params': params,
        'cif_files': cif_files,
        'confidence': metrics,
    }

    summary_path = os.path.join(output_dir, 'chai_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"\n  Summary written to: chai_summary.json")
    print(f"  CIF files: {len(cif_files)}")
    print(f"  Confidence entries: {len(metrics)}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Run Chai-1 structure prediction')
    parser.add_argument('--input',                    required=True,       help='Input FASTA file from node 02')
    parser.add_argument('--output-dir',               default='./chai_output', help='Output directory (must be empty; chai-lab requirement)')
    parser.add_argument('--num-trunk-recycles',       type=int, default=3, help='Number of trunk recycling steps')
    parser.add_argument('--num-diffusion-timesteps',  type=int, default=200, help='Number of diffusion timesteps')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Input:                  {args.input}")
    print(f"Trunk recycles:         {args.num_trunk_recycles}")
    print(f"Diffusion timesteps:    {args.num_diffusion_timesteps}")

    run_chai_fold(
        input_file=args.input,
        output_dir=args.output_dir,
        num_trunk_recycles=args.num_trunk_recycles,
        num_diffusion_timesteps=args.num_diffusion_timesteps,
    )

    print("\nCollecting outputs:")
    collected = collect_outputs(args.output_dir)

    print("\nParsing confidence metrics:")
    metrics = parse_confidence_files(args.output_dir)

    cif_files = [f for f in collected if f.endswith('.cif')]
    if not cif_files:
        print("Warning: No CIF files found in output")
    if not metrics:
        print("Warning: No scores NPZ files found")

    params = {
        'num_trunk_recycles': args.num_trunk_recycles,
        'num_diffusion_timesteps': args.num_diffusion_timesteps,
    }
    write_summary(args.output_dir, collected, metrics, params)

    # Copy outputs to workspace root so Silva's output glob picks them up
    print("\nCopying outputs to workspace root:")
    for pattern in ["pred.*.cif", "scores.*.npz", "chai_summary.json"]:
        for src in Path(args.output_dir).glob(pattern):
            shutil.copy(src, ".")
            print(f"  Copied: {src.name}")

    print("\nNode 04 completed ✓")


if __name__ == '__main__':
    main()
