#!/usr/bin/env python3
"""
Validate BoltzGen target structure and design specification.

Checks:
- Valid design spec YAML
- Target structure file exists (if referenced)
- Required fields for binder design
"""

import argparse
import json
import os
import shutil
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None


def validate_design_spec(filepath):
    """Validate a BoltzGen design specification YAML file."""
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Design spec not found: {filepath}")

    summary = {
        'design_spec': path.name,
        'entities': [],
        'design_chains': [],
        'valid': True,
    }

    if yaml is not None:
        with open(path, 'r') as f:
            data = yaml.safe_load(f)

        if not isinstance(data, dict):
            raise ValueError("YAML root must be a mapping")

        sequences = data.get('sequences', [])
        if not sequences:
            raise ValueError("No sequences defined in design spec")

        for i, seq_entry in enumerate(sequences):
            if not isinstance(seq_entry, dict):
                continue
            for entity_type, entity_data in seq_entry.items():
                entity_id = entity_data.get('id', f'chain_{i}')
                entity_info = {
                    'type': entity_type,
                    'id': str(entity_id),
                }

                if 'sequence' in entity_data:
                    entity_info['length'] = len(entity_data['sequence'])

                if 'design' in entity_data:
                    design = entity_data['design']
                    entity_info['design_type'] = design.get('type', 'unknown')
                    summary['design_chains'].append(str(entity_id))

                    if design.get('type') == 'de_novo':
                        entity_info['min_length'] = design.get('min_length')
                        entity_info['max_length'] = design.get('max_length')
                    elif design.get('type') == 'scaffold_library':
                        entity_info['scaffold_path'] = design.get('path', '')

                summary['entities'].append(entity_info)
    else:
        with open(path, 'r') as f:
            content = f.read()
        if 'sequences' not in content:
            raise ValueError("Missing 'sequences' field in design spec")
        if 'design' not in content:
            raise ValueError("No design specification found")
        summary['entities'].append({'type': 'unknown', 'note': 'PyYAML not available'})

    return summary


def main():
    parser = argparse.ArgumentParser(description='Validate BoltzGen target and design spec')
    parser.add_argument('--design-spec', required=True, help='Design specification YAML file')
    parser.add_argument('--target-structure', default='', help='Target structure file (PDB/CIF)')
    parser.add_argument('--output-dir', default='.', help='Output directory')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    print(f"Validating design spec: {args.design_spec}")
    summary = validate_design_spec(args.design_spec)

    # Copy design spec to output
    spec_dest = os.path.join(args.output_dir, 'design_spec.yaml')
    shutil.copy2(args.design_spec, spec_dest)
    print(f"Design spec copied to: {spec_dest}")

    # Copy target structure if provided
    if args.target_structure and os.path.exists(args.target_structure):
        target_dest = os.path.join(args.output_dir, os.path.basename(args.target_structure))
        shutil.copy2(args.target_structure, target_dest)
        summary['target_structure'] = os.path.basename(args.target_structure)
        print(f"Target structure copied to: {target_dest}")

    # Copy scaffold directory if referenced
    spec_dir = os.path.dirname(os.path.abspath(args.design_spec))
    for entity in summary.get('entities', []):
        scaffold_path = entity.get('scaffold_path', '')
        if scaffold_path:
            src_scaffold = os.path.join(spec_dir, scaffold_path)
            if os.path.isdir(src_scaffold):
                dest_scaffold = os.path.join(args.output_dir, scaffold_path)
                if not os.path.exists(dest_scaffold):
                    shutil.copytree(src_scaffold, dest_scaffold)
                print(f"Scaffold library copied to: {dest_scaffold}")

    # Write summary
    summary_path = os.path.join(args.output_dir, 'target_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Target summary written to: {summary_path}")

    for entity in summary['entities']:
        print(f"  Entity: {entity.get('id', 'N/A')} ({entity['type']}"
              f"{', design=' + entity.get('design_type', '') if 'design_type' in entity else ''})")

    if summary['design_chains']:
        print(f"Design chains: {', '.join(summary['design_chains'])}")

    print("Validation passed")


if __name__ == '__main__':
    main()
