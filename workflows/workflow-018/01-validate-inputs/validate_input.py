#!/usr/bin/env python3
"""
Validate and normalise inputs for Boltz-2 and Chai-1.

Accepted input types
--------------------
- Boltz YAML  (version/sequences structure)
- Chai FASTA  (>protein|name=... headers)
- UniProt / generic FASTA (>sp|..., >tr|..., or plain >id headers)
  → auto-converted to boltz_input.yaml + chai_input.fasta

Outputs
-------
- boltz_input.yaml   – Boltz-2 ready input
- chai_input.fasta   – Chai-1 ready input
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


# ── UniProt / generic FASTA → Boltz YAML + Chai FASTA ───────────────────────

def _is_uniprot_or_generic_fasta(header):
    """Return True for non-Chai FASTA headers (UniProt sp/tr, NCBI, plain)."""
    body = header.lstrip('>')
    first_token = body.split('|')[0].strip().lower()
    chai_types = {'protein', 'nucleic-acid', 'ligand', 'residue', 'rna', 'dna'}
    return first_token not in chai_types


def _parse_fasta_records(filepath):
    """Return list of (header, sequence) tuples from a FASTA file."""
    content = Path(filepath).read_text().strip()
    records = []
    for block in re.split(r'\n(?=>)', content):
        lines = block.strip().splitlines()
        records.append((lines[0], ''.join(lines[1:]).strip()))
    return records


def _uniprot_id_from_header(header):
    """Extract accession from >sp|P69905|... or fall back to first token."""
    body = header.lstrip('>')
    parts = body.split('|')
    if len(parts) >= 2 and parts[0].lower() in ('sp', 'tr'):
        return parts[1]
    return parts[0].split()[0]


def convert_generic_fasta(filepath):
    """Convert UniProt/generic FASTA to boltz_input.yaml and chai_input.fasta."""
    records = _parse_fasta_records(filepath)
    if not records:
        raise ValueError("FASTA file contains no records")

    boltz_sequences = []
    chai_lines = []
    for i, (header, seq) in enumerate(records):
        if not seq:
            raise ValueError(f"Empty sequence in record: {header}")
        chain_id = chr(ord('A') + i)
        acc = _uniprot_id_from_header(header)
        boltz_sequences.append({'protein': {'id': chain_id, 'sequence': seq}})
        chai_lines.append(f'>protein|name={acc}')
        chai_lines.append(seq)

    # Write boltz_input.yaml
    boltz_data = {'version': 1, 'sequences': boltz_sequences}
    if yaml is not None:
        with open('boltz_input.yaml', 'w') as f:
            yaml.dump(boltz_data, f, default_flow_style=False, allow_unicode=True)
    else:
        with open('boltz_input.yaml', 'w') as f:
            f.write('version: 1\nsequences:\n')
            for entry in boltz_sequences:
                p = entry['protein']
                f.write(f"- protein:\n    id: {p['id']}\n    sequence: {p['sequence']}\n")

    # Write chai_input.fasta
    Path('chai_input.fasta').write_text('\n'.join(chai_lines) + '\n')

    print(f"  Converted {len(records)} sequence(s) from generic FASTA")
    for i, (header, seq) in enumerate(records):
        acc = _uniprot_id_from_header(header)
        print(f"  Chain {chr(ord('A') + i)}: {acc} ({len(seq)} aa)")


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
        shutil.copy2(args.yaml, 'boltz_input.yaml')
        print("  Copied to: boltz_input.yaml")
        for e in summary['entities']:
            print(f"  Entity: {e.get('id')} ({e['type']}, length={e.get('length', 'N/A')})")
        print("  Boltz YAML validation passed ✓")

    if args.fasta:
        print(f"\nProcessing FASTA: {args.fasta}")
        # Peek at header to decide format
        first_line = Path(args.fasta).read_text().splitlines()[0]
        if _is_uniprot_or_generic_fasta(first_line):
            print("  Detected UniProt/generic FASTA — converting for Boltz + Chai")
            convert_generic_fasta(args.fasta)
            print("  Wrote boltz_input.yaml + chai_input.fasta ✓")
        else:
            print("  Detected Chai-format FASTA — validating")
            summary = validate_chai_fasta(args.fasta)
            shutil.copy2(args.fasta, 'chai_input.fasta')
            print("  Copied to: chai_input.fasta")
            for e in summary['entities']:
                print(f"  Entity: {e.get('id')} ({e['type']}, length={e.get('length', 'N/A')})")
            print("  Chai FASTA validation passed ✓")

    print("\nNode 01 completed ✓")


if __name__ == '__main__':
    main()
