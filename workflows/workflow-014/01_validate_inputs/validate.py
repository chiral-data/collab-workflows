#!/usr/bin/env python3
"""Validate input CSV and standardize SMILES strings using RDKit + pandas."""

import argparse
import json
import os
import sys

import pandas as pd

try:
    from rdkit import Chem
    from rdkit.Chem.MolStandardize import rdMolStandardize
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False
    print("WARNING: RDKit not available, falling back to basic validation", flush=True)


def standardize_smiles(smiles):
    """Validate, standardize, and canonicalize a SMILES string.

    Returns canonical SMILES or None if invalid.
    """
    if not smiles or not str(smiles).strip():
        return None

    smiles = str(smiles).strip()

    if HAS_RDKIT:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        # Standardize: remove fragments, neutralize, canonicalize
        try:
            uncharger = rdMolStandardize.Uncharger()
            mol = uncharger.uncharge(mol)
            mol = rdMolStandardize.Cleanup(mol)
        except Exception:
            pass  # If standardization fails, use the parsed mol as-is
        return Chem.MolToSmiles(mol)
    else:
        # Basic fallback: check length and printable characters
        if len(smiles) < 1 or len(smiles) > 5000:
            return None
        return smiles


def main():
    parser = argparse.ArgumentParser(description="Validate and standardize SMILES input")
    parser.add_argument("--input", required=True, help="Input CSV file path")
    parser.add_argument("--smiles-column", default="smiles", help="Name of the SMILES column")
    args = parser.parse_args()

    input_path = args.input
    smiles_column = args.smiles_column

    if not os.path.exists(input_path):
        print(f"ERROR: Input file '{input_path}' not found", flush=True)
        sys.exit(1)

    print(f"Reading input file: {input_path}", flush=True)
    print(f"SMILES column: {smiles_column}", flush=True)
    print(f"RDKit available: {HAS_RDKIT}", flush=True)

    df = pd.read_csv(input_path)
    total_input = len(df)
    print(f"Total rows in input: {total_input}", flush=True)

    if smiles_column not in df.columns:
        print(f"ERROR: Column '{smiles_column}' not found. Available: {list(df.columns)}", flush=True)
        sys.exit(1)

    # Validate and standardize SMILES
    df["_canonical_smiles"] = df[smiles_column].apply(standardize_smiles)

    invalid_mask = df["_canonical_smiles"].isna()
    invalid_count = int(invalid_mask.sum())
    df = df[~invalid_mask].copy()

    # Remove duplicates by canonical SMILES
    before_dedup = len(df)
    df = df.drop_duplicates(subset="_canonical_smiles", keep="first")
    duplicate_count = before_dedup - len(df)

    # Replace original SMILES with canonical version
    df[smiles_column] = df["_canonical_smiles"]
    df = df.drop(columns=["_canonical_smiles"])

    print(f"Valid unique molecules: {len(df)}", flush=True)
    print(f"Invalid SMILES removed: {invalid_count}", flush=True)
    print(f"Duplicates removed: {duplicate_count}", flush=True)

    # Write standardized CSV
    os.makedirs("outputs", exist_ok=True)
    df.to_csv("outputs/standardized_molecules.csv", index=False)

    # Write validation report
    report = {
        "total_input": total_input,
        "valid_unique": len(df),
        "invalid_removed": invalid_count,
        "duplicates_removed": duplicate_count,
        "smiles_column": smiles_column,
        "input_file": os.path.basename(input_path),
        "rdkit_available": HAS_RDKIT,
    }
    with open("outputs/validation_report.json", "w") as f:
        json.dump(report, f, indent=2)

    print("Validation complete.", flush=True)


if __name__ == "__main__":
    main()
