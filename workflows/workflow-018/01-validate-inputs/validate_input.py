#!/usr/bin/env python3
"""
Validate inputs for Boltz-2 (YAML) and Chai-1 (FASTA).

Boltz checks:
- Valid YAML format
- Required fields (version, sequences)
- Valid entity types (protein, rna, dna, ligand)
- Non-empty sequences

Chai checks:
- Valid FASTA format
- Chai-style headers (>protein|name=<id>)
- Valid entity types (protein, nucleic-acid, ligand, residue)
- Non-empty sequences
"""

import argparse
import re
import shutil
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None


# ── Boltz (YAML) ─────────────────────────────────────────────────────────────

def parse_yaml_simple(filepath):
    """Fallback YAML parser when PyYAML is unavailable."""
    with open(filepath, 'r') as f:
        content = f.read()
    for field in ('version', 'sequences'):
        if field not in content:
            raise ValueError(f"Missing '{field}' field")
    if 'sequence' not in content and 'smiles' not in content:
        raise ValueError("No sequence or smiles data found")
    return content


def validate_boltz_yaml(filepath):
    """Validate a Boltz-2 input YAML file."""
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    if path.suffix not in ('.yaml', '.yml'):
        raise ValueError(f"Expected YAML file, got: {path.suffix}")

    summary = {'filename': path.name, 'format': 'yaml', 'entities': [], 'valid': True}

    if yaml is not None:
        with open(path, 'r') as f:
            data = yaml.safe_load(f)

        if not isinstance(data, dict):
            raise ValueError("YAML root must be a mapping")
        if data.get('version') is None:
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
                    raise ValueError(f"Unknown entity type: '{entity_type}'")
                entity_info = {
                    'type': entity_type,
                    'id': entity_data.get('id', f'chain_{i}'),
                }
                if entity_type in ('protein', 'rna', 'dna'):
                    seq = entity_data.get('sequence', '')
                    if not seq:
                        raise ValueError(f"Empty sequence for entity '{entity_info['id']}'")
                    entity_info['length'] = len(seq)
                elif entity_type == 'ligand':
                    entity_info['smiles'] = entity_data.get('smiles', '')
                    entity_info['ccd'] = entity_data.get('ccd', '')
                summary['entities'].append(entity_info)
    else:
        parse_yaml_simple(path)
        summary['entities'].append({'type': 'unknown', 'note': 'PyYAML not available'})

    return summary


# ── Chai (FASTA) ─────────────────────────────────────────────────────────────

def validate_chai_fasta(filepath):
    """Validate a Chai-1 input FASTA file.

    Chai expects headers like: >protein|name=chainA
    Valid entity types: protein, nucleic-acid, ligand, residue
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"File not found: {filepath}")
    if path.suffix not in ('.fasta', '.fa'):
        raise ValueError(f"Expected FASTA file, got: {path.suffix}")

    summary = {'filename': path.name, 'format': 'fasta', 'entities': [], 'valid': True}

    content = path.read_text().strip()
    if not content:
        raise ValueError("FASTA file is empty")
    if not content.startswith('>'):
        raise ValueError("FASTA file must start with a '>' header line")

    valid_types = {'protein', 'nucleic-acid', 'ligand', 'residue'}
    records = re.split(r'\n(?=>)', content)

    for record in records:
        lines = record.strip().splitlines()
        header = lines[0]
        sequence = ''.join(lines[1:]).strip()

        if not header.startswith('>'):
            raise ValueError(f"Invalid header (missing '>'): {header}")

        # Parse >protein|name=chainA style header
        header_body = header[1:]
        parts = header_body.split('|')
        entity_type = parts[0].strip().lower()

        if entity_type not in valid_types:
            raise ValueError(
                f"Unknown entity type '{entity_type}' in header '{header}'. "
                f"Expected one of: {', '.join(valid_types)}"
            )

        name = None
        for part in parts[1:]:
            if part.startswith('name='):
                name = part.split('=', 1)[1]

        entity_info = {'type': entity_type, 'id': name or header_body}

        if entity_type in ('protein', 'nucleic-acid'):
            if not sequence:
                raise ValueError(f"Empty sequence for entity '{entity_info['id']}'")
            entity_info['length'] = len(sequence)
        elif entity_type == 'ligand':
            entity_info['smiles'] = sequence

        summary['entities'].append(entity_info)

    return summary


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Validate Boltz-2 (YAML) and Chai-1 (FASTA) inputs')
    parser.add_argument('--yaml',  help='Boltz-2 input YAML file')
    parser.add_argument('--fasta', help='Chai-1 input FASTA file')
    args = parser.parse_args()

    if not args.yaml and not args.fasta:
        print("Error: provide at least one of --yaml or --fasta")
        sys.exit(1)

    if args.yaml:
        print(f"\nValidating Boltz YAML: {args.yaml}")
        summary = validate_boltz_yaml(args.yaml)
        shutil.copy2(args.yaml, 'validated_boltz_input.yaml')
        print("  Copied to: validated_boltz_input.yaml")
        for e in summary['entities']:
            print(f"  Entity: {e.get('id')} ({e['type']}, length={e.get('length', 'N/A')})")
        print("  Boltz YAML validation passed ✓")

    if args.fasta:
        print(f"\nValidating Chai FASTA: {args.fasta}")
        summary = validate_chai_fasta(args.fasta)
        shutil.copy2(args.fasta, 'validated_chai_input.fasta')
        print("  Copied to: validated_chai_input.fasta")
        for e in summary['entities']:
            print(f"  Entity: {e.get('id')} ({e['type']}, length={e.get('length', 'N/A')})")
        print("  Chai FASTA validation passed ✓")

    print("\nNode 01 completed ✓")


if __name__ == '__main__':
    main()
