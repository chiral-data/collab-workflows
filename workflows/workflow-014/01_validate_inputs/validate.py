#!/usr/bin/env python3
"""Validate input CSV and standardize SMILES strings."""

import csv
import json
import os
import sys

try:
    from rdkit import Chem
except ImportError:
    Chem = None


def load_params():
    global_params = {}
    for path in ["./inputs/global_params.json", "../global_params.json"]:
        if os.path.exists(path):
            with open(path) as f:
                global_params = json.load(f)
            break
    return global_params


def validate_smiles(smiles):
    """Validate and canonicalize a SMILES string. Returns canonical SMILES or None."""
    if not smiles or not smiles.strip():
        return None
    smiles = smiles.strip()
    if Chem is not None:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol)
    # Fallback: basic string validation if RDKit unavailable
    if len(smiles) < 1 or len(smiles) > 5000:
        return None
    return smiles


def main():
    params = load_params()
    input_file = params.get("input_file", "drugbank_approved.csv")
    smiles_column = params.get("smiles_column", "smiles")

    # Find input file
    input_path = None
    for candidate in [
        f"./inputs/{input_file}",
        f"./inputs/input_files/{input_file}",
        f"../input_files/{input_file}",
    ]:
        if os.path.exists(candidate):
            input_path = candidate
            break

    if input_path is None:
        print(f"ERROR: Input file '{input_file}' not found", flush=True)
        sys.exit(1)

    print(f"Reading input file: {input_path}", flush=True)

    valid_rows = []
    invalid_count = 0
    duplicate_count = 0
    seen_smiles = set()

    with open(input_path, newline="") as f:
        reader = csv.DictReader(f)

        if smiles_column not in reader.fieldnames:
            print(f"ERROR: Column '{smiles_column}' not found. Available: {reader.fieldnames}", flush=True)
            sys.exit(1)

        fieldnames = reader.fieldnames

        for row in reader:
            raw_smiles = row.get(smiles_column, "")
            canonical = validate_smiles(raw_smiles)

            if canonical is None:
                invalid_count += 1
                continue

            if canonical in seen_smiles:
                duplicate_count += 1
                continue

            seen_smiles.add(canonical)
            row[smiles_column] = canonical
            valid_rows.append(row)

    total_input = len(valid_rows) + invalid_count + duplicate_count
    print(f"Total input molecules: {total_input}", flush=True)
    print(f"Valid unique molecules: {len(valid_rows)}", flush=True)
    print(f"Invalid SMILES removed: {invalid_count}", flush=True)
    print(f"Duplicates removed: {duplicate_count}", flush=True)

    # Write validated CSV
    os.makedirs("outputs", exist_ok=True)
    with open("outputs/validated_molecules.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(valid_rows)

    # Write validation report
    report = {
        "total_input": total_input,
        "valid_unique": len(valid_rows),
        "invalid_removed": invalid_count,
        "duplicates_removed": duplicate_count,
        "smiles_column": smiles_column,
        "input_file": input_file,
    }
    with open("outputs/validation_report.json", "w") as f:
        json.dump(report, f, indent=2)

    print("Validation complete.", flush=True)


if __name__ == "__main__":
    main()
