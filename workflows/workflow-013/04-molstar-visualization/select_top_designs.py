#!/usr/bin/env python3
"""
Select top-ranked BoltzGen designed binder CIF files for Mol* visualization.

Ranks designs by confidence/iPTM score from metrics CSV (or alphabetically
by filename if no CSV available) and copies the top N to output with
standardized naming for the platform Mol* viewer.
"""

import argparse
import csv
import glob
import json
import os
import re
import shutil
from pathlib import Path


def load_design_scores(input_dir):
    """Load and rank designs by score from metrics CSV or CIF files."""
    designs = []

    # Try loading from aggregate metrics CSV
    csv_files = glob.glob(os.path.join(input_dir, 'aggregate_metrics*.csv'))
    if csv_files:
        with open(csv_files[0], 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                design_id = row.get('design_id', row.get('name', ''))
                score = float(row.get('confidence_score', row.get('iptm', 0.0)))

                # Find CIF file matching design_id
                cif_candidates = glob.glob(os.path.join(input_dir, f'*{design_id}*.cif'))
                cif_file = cif_candidates[0] if cif_candidates else None
                if cif_file:
                    designs.append({
                        'design_id': design_id,
                        'score': score,
                        'cif_file': cif_file,
                    })
    else:
        # Fall back to listing CIF files
        for cif_path in sorted(glob.glob(os.path.join(input_dir, '*.cif'))):
            designs.append({
                'design_id': Path(cif_path).stem,
                'score': 0.0,
                'cif_file': cif_path,
            })

    designs.sort(key=lambda x: x['score'], reverse=True)
    return designs


def main():
    parser = argparse.ArgumentParser(description='Select top BoltzGen designs for Mol* visualization')
    parser.add_argument('--input-dir', default='.', help='Directory with design results')
    parser.add_argument('--top-n', type=int, default=5, help='Number of top designs to select')
    parser.add_argument('--output-dir', default='.', help='Output directory')
    args = parser.parse_args()

    designs = load_design_scores(args.input_dir)
    if not designs:
        print("Warning: No designs found")
        return

    top_n = min(args.top_n, len(designs))
    selected = designs[:top_n]

    os.makedirs(args.output_dir, exist_ok=True)

    for i, design in enumerate(selected):
        rank = i + 1
        dest = os.path.join(args.output_dir, f"top_design_{rank}_{design['design_id']}.cif")
        shutil.copy2(design['cif_file'], dest)
        print(f"  Rank {rank}: {design['design_id']} (score={design['score']:.4f}) -> {os.path.basename(dest)}")

    # Write Mol* viewer configuration
    config = {
        'viewer': 'molstar',
        'structures': [
            {
                'rank': i + 1,
                'file': f"top_design_{i+1}_{d['design_id']}.cif",
                'label': f"Design {d['design_id']} (Rank {i+1})",
                'score': d['score'],
                'format': 'cif',
                'style': 'cartoon',
                'color_scheme': 'chain',
            }
            for i, d in enumerate(selected)
        ],
    }

    config_path = os.path.join(args.output_dir, 'molstar_config.json')
    with open(config_path, 'w') as f:
        json.dump(config, f, indent=2)
    print(f"Mol* config written to: {config_path}")


if __name__ == '__main__':
    main()
