#!/usr/bin/env python3
"""
Validate Boltz-2 input YAML file.

Checks:
- Valid YAML format
- Required fields (version, sequences)
- Valid entity types (protein, rna, dna, ligand)
- Non-empty sequences
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
    # Fallback: simple YAML-like parsing for basic cases
    yaml = None


def parse_yaml_simple(filepath):
    """Simple YAML parser for Boltz input files when PyYAML is not available."""
    import re
    with open(filepath, 'r') as f:
        content = f.read()

    # Basic validation: check for required keywords
    if 'version' not in content:
        raise ValueError("Missing 'version' field")
    if 'sequences' not in content:
        raise ValueError("Missing 'sequences' field")
    if 'sequence' not in content and 'smiles' not in content:
        raise ValueError("No sequence or smiles data found")

    return content


def validate_boltz_yaml(filepath):
    """Validate a Boltz-2 input YAML file."""
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {filepath}")

    if not path.suffix in ('.yaml', '.yml'):
        raise ValueError(f"Expected YAML file, got: {path.suffix}")

    summary = {
        'filename': path.name,
        'entities': [],
        'valid': True,
    }

    if yaml is not None:
        with open(path, 'r') as f:
            data = yaml.safe_load(f)

        if not isinstance(data, dict):
            raise ValueError("YAML root must be a mapping")

        version = data.get('version', None)
        if version is None:
            raise ValueError("Missing required field: 'version'")

        sequences = data.get('sequences', [])
        if not sequences:
            raise ValueError("No sequences defined in input file")

        valid_types = {'protein', 'rna', 'dna', 'ligand', 'ccd'}

        for i, seq_entry in enumerate(sequences):
            if not isinstance(seq_entry, dict):
                raise ValueError(f"Sequence entry {i} must be a mapping")

            for entity_type, entity_data in seq_entry.items():
                if entity_type not in valid_types:
                    raise ValueError(f"Unknown entity type: {entity_type}")

                entity_info = {
                    'type': entity_type,
                    'id': entity_data.get('id', f'chain_{i}'),
                }

                if entity_type == 'protein':
                    seq = entity_data.get('sequence', '')
                    if not seq:
                        raise ValueError(f"Empty protein sequence for entity {entity_info['id']}")
                    entity_info['length'] = len(seq)
                elif entity_type in ('rna', 'dna'):
                    seq = entity_data.get('sequence', '')
                    if not seq:
                        raise ValueError(f"Empty {entity_type} sequence for entity {entity_info['id']}")
                    entity_info['length'] = len(seq)
                elif entity_type == 'ligand':
                    entity_info['smiles'] = entity_data.get('smiles', '')
                    entity_info['ccd'] = entity_data.get('ccd', '')

                summary['entities'].append(entity_info)
    else:
        parse_yaml_simple(path)
        summary['entities'].append({'type': 'unknown', 'note': 'PyYAML not available for detailed parsing'})

    return summary


def main():
    parser = argparse.ArgumentParser(description='Validate Boltz-2 input YAML')
    parser.add_argument('--input', required=True, help='Input YAML file')
    parser.add_argument('--output', required=True, help='Output validated YAML file')
    args = parser.parse_args()

    print(f"Validating input: {args.input}")
    summary = validate_boltz_yaml(args.input)

    # Copy validated file to output
    shutil.copy2(args.input, args.output)
    print(f"Validated input copied to: {args.output}")

    # Write summary next to the validated output file
    output_dir = str(Path(args.output).parent)
    summary_path = os.path.join(output_dir, 'input_summary.json')
    with open(summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    print(f"Input summary written to: {summary_path}")

    for entity in summary['entities']:
        entity_type = entity['type']
        entity_id = entity.get('id', 'N/A')
        length = entity.get('length', 'N/A')
        print(f"  Entity: {entity_id} ({entity_type}, length={length})")

    print("Validation passed")


if __name__ == '__main__':
    main()
