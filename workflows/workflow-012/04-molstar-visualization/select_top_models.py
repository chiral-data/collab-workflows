#!/usr/bin/env python3
"""
Select top-ranked Boltz-2 model PDB files for Mol* visualization.

Ranks models by confidence score and copies the top N to output
with standardized naming for the platform Mol* viewer.
"""

import argparse
import glob
import json
import os
import re
import shutil


def load_confidence_scores(input_dir):
    """Load and rank models by confidence score."""
    models = []
    for conf_file in sorted(glob.glob(os.path.join(input_dir, 'confidence_*_model_*.json'))):
        match = re.search(r'confidence_(.+)_model_(\d+)\.json', os.path.basename(conf_file))
        if not match:
            continue

        protein_name = match.group(1)
        model_id = int(match.group(2))

        with open(conf_file, 'r') as f:
            data = json.load(f)

        pdb_file = os.path.join(input_dir, f'{protein_name}_model_{model_id}.pdb')
        if not os.path.exists(pdb_file):
            continue

        models.append({
            'protein_name': protein_name,
            'model_id': model_id,
            'confidence_score': data.get('confidence_score', 0.0),
            'pdb_file': pdb_file,
        })

    models.sort(key=lambda x: x['confidence_score'], reverse=True)
    return models


def main():
    parser = argparse.ArgumentParser(description='Select top Boltz-2 models for Mol* visualization')
    parser.add_argument('--input-dir', default='.', help='Directory with prediction results')
    parser.add_argument('--top-n', type=int, default=5, help='Number of top models to select')
    parser.add_argument('--output-dir', default='.', help='Output directory')
    args = parser.parse_args()

    models = load_confidence_scores(args.input_dir)
    if not models:
        print("Warning: No models found with confidence scores")
        return

    top_n = min(args.top_n, len(models))
    selected = models[:top_n]

    os.makedirs(args.output_dir, exist_ok=True)

    # Copy top models with standardized names
    for i, model in enumerate(selected):
        rank = i + 1
        dest = os.path.join(args.output_dir,
                            f"top_model_{rank}_{model['protein_name']}_model_{model['model_id']}.pdb")
        shutil.copy2(model['pdb_file'], dest)
        print(f"  Rank {rank}: Model {model['model_id']} "
              f"(confidence={model['confidence_score']:.4f}) -> {os.path.basename(dest)}")

    # Write Mol* viewer configuration
    config = {
        'viewer': 'molstar',
        'structures': [
            {
                'rank': i + 1,
                'file': f"top_model_{i+1}_{m['protein_name']}_model_{m['model_id']}.pdb",
                'label': f"Model {m['model_id']} (Rank {i+1})",
                'confidence': m['confidence_score'],
                'format': 'pdb',
                'style': 'cartoon',
                'color_scheme': 'plddt',
            }
            for i, m in enumerate(selected)
        ],
    }

    config_path = os.path.join(args.output_dir, 'molstar_config.json')
    with open(config_path, 'w') as f:
        json.dump(config, f, indent=2)
    print(f"Mol* config written to: {config_path}")


if __name__ == '__main__':
    main()
