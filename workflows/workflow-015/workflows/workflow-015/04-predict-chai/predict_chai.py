#!/usr/bin/env python3
"""
Node 04: Chai-1 structure prediction.

Wraps the `chai-lab fold` CLI, collects output CIF files and
confidence metrics (pLDDT, pTM, ipTM) into chai_summary.json
for use by the visualization node.
"""

import argparse
import glob
import json
import os
import subprocess
import sys
from pathlib import Path


# ── Prediction ────────────────────────────────────────────────────────────────

def run_chai_fold(input_file, output_dir, num_trunk_recycles, num_diffusion_timesteps):
    """Run the chai-lab fold CLI command."""
    cmd = [
        "chai-lab", "fold", input_file,
        "--output-dir", output_dir,
        "--num-trunk-recycles", str(num_trunk_recycles),
        "--num-diffusion-timesteps", str(num_diffusion_timesteps),
    ]

    print(f"\nRunning: {' '.join(cmd)}\n")
    result = subprocess.run(cmd, capture_output=False, text=True)

    if result.returncode != 0:
        print(f"Error: Chai-1 failed with return code {result.returncode}")
        sys.exit(1)


# ── Output collection ─────────────────────────────────────────────────────────

def collect_outputs(output_dir):
    """Collect CIF and scores JSON files from Chai output directory."""
    collected = []
    for pattern in ["*.cif", "scores.*.json"]:
        for filepath in glob.glob(os.path.join(output_dir, pattern)):
            collected.append(os.path.basename(filepath))
            print(f"  Found: {os.path.basename(filepath)}")
    return collected


# ── Confidence parsing ────────────────────────────────────────────────────────

def parse_confidence_files(output_dir):
    """Extract pLDDT, pTM, ipTM from Chai scores JSONs.

    Chai outputs scores.model_idx_N.json with:
      - plddt: list of per-token values (0-100 scale)
      - ptm: float (0-1)
      - iptm: float (0-1)
      - aggregate_score: float
    """
    score_files = sorted(glob.glob(os.path.join(output_dir, "scores.*.json")))
    metrics = []

    for sf in score_files:
        with open(sf, 'r') as f:
            data = json.load(f)

        # Extract model index from filename: scores.model_idx_0.json → model_0
        stem = Path(sf).stem  # scores.model_idx_0
        parts = stem.split('.')
        model_label = parts[-1].replace('model_idx_', 'model_') if len(parts) > 1 else stem

        plddt_raw = data.get('plddt', None)
        if isinstance(plddt_raw, list) and plddt_raw:
            plddt = round(sum(plddt_raw) / len(plddt_raw) / 100.0, 4)
        elif isinstance(plddt_raw, (int, float)):
            plddt = round(float(plddt_raw) / 100.0, 4) if plddt_raw > 1.0 else round(float(plddt_raw), 4)
        else:
            plddt = None

        entry = {
            'sample': model_label,
            'plddt': plddt,
            'ptm': data.get('ptm', None),
            'iptm': data.get('iptm', None),
            'aggregate_score': data.get('aggregate_score', None),
        }
        metrics.append(entry)
        print(f"  {model_label}: pLDDT={entry['plddt']}, pTM={entry['ptm']}, ipTM={entry['iptm']}")

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
    parser.add_argument('--output-dir',               default='.',         help='Output directory')
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
        print("Warning: No scores JSON files found")

    params = {
        'num_trunk_recycles': args.num_trunk_recycles,
        'num_diffusion_timesteps': args.num_diffusion_timesteps,
    }
    write_summary(args.output_dir, collected, metrics, params)

    print("\nNode 04 completed ✓")


if __name__ == '__main__':
    main()
