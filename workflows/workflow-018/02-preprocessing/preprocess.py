#!/usr/bin/env python3
"""
Preprocessing node: convert validated YAML input into:
  - Boltz-2 format (YAML)     → boltz_input.yaml
  - Chai-1 format  (FASTA)    → chai_input.fasta

Boltz YAML structure:
  version: 1
  sequences:
    - protein:
        id: A
        sequence: MKTLL...

Chai FASTA structure:
  >protein|name=A
  MKTLL...
"""

import argparse
import shutil
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None


# ── Readers ──────────────────────────────────────────────────────────────────

def load_boltz_yaml(filepath):
    """Load and return parsed Boltz-2 YAML data."""
    path = Path(filepath)
    if yaml is not None:
        with open(path, 'r') as f:
            return yaml.safe_load(f)
    else:
        raise ImportError("PyYAML is required for preprocessing. Install with: pip install pyyaml")


# ── Writers ──────────────────────────────────────────────────────────────────

def write_boltz_yaml(data, output_path):
    """Write Boltz-2 YAML output (essentially a clean copy of the input)."""
    with open(output_path, 'w') as f:
        yaml.dump(data, f, default_flow_style=False, allow_unicode=True)
    print(f"  Boltz YAML written to: {output_path}")


def write_chai_fasta(data, output_path):
    """Convert Boltz YAML structure to Chai-1 FASTA format.

    Chai FASTA header format: >entity_type|name=<id>
    Supported mappings:
      protein → >protein|name=<id>
      rna     → >nucleic-acid|name=<id>
      dna     → >nucleic-acid|name=<id>
      ligand  → >ligand|name=<id>  (SMILES as sequence)
    """
    type_map = {
        'protein': 'protein',
        'rna':     'nucleic-acid',
        'dna':     'nucleic-acid',
        'ligand':  'ligand',
        'ccd':     'ligand',
    }

    lines = []
    sequences = data.get('sequences', [])

    for i, seq_entry in enumerate(sequences):
        for entity_type, entity_data in seq_entry.items():
            chai_type = type_map.get(entity_type)
            if chai_type is None:
                print(f"  Warning: skipping unknown entity type '{entity_type}'")
                continue

            entity_id = entity_data.get('id', f'chain_{i}')
            header = f">{chai_type}|name={entity_id}"

            if entity_type in ('protein', 'rna', 'dna'):
                sequence = entity_data.get('sequence', '')
                if not sequence:
                    raise ValueError(f"Empty sequence for entity '{entity_id}'")
            elif entity_type in ('ligand', 'ccd'):
                # Chai uses SMILES string as the sequence for ligands
                sequence = entity_data.get('smiles', '')
                if not sequence:
                    raise ValueError(f"No SMILES string for ligand '{entity_id}'")

            lines.append(header)
            lines.append(sequence)

    with open(output_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')

    print(f"  Chai FASTA written to: {output_path}")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Preprocess input for Boltz-2 and Chai-1')
    parser.add_argument('--input',        required=True, help='Validated input YAML file')
    parser.add_argument('--boltz-output', required=True, help='Output path for Boltz-2 YAML')
    parser.add_argument('--chai-output',  required=True, help='Output path for Chai-1 FASTA')
    args = parser.parse_args()

    print(f"\nLoading: {args.input}")
    data = load_boltz_yaml(args.input)

    print("\nWriting outputs:")
    write_boltz_yaml(data, args.boltz_output)
    write_chai_fasta(data, args.chai_output)

    print("\nNode 02 completed ✓")


if __name__ == '__main__':
    main()
