#!/usr/bin/env python3
"""Run ADMET-AI predictions on validated molecules."""

import json
import os
import subprocess
import sys


def load_params():
    global_params = {}
    for path in ["./inputs/global_params.json", "../global_params.json"]:
        if os.path.exists(path):
            with open(path) as f:
                global_params = json.load(f)
            break
    return global_params


def main():
    params = load_params()
    smiles_column = params.get("smiles_column", "smiles")

    input_path = "./inputs/standardized_molecules.csv"
    if not os.path.exists(input_path):
        print("ERROR: standardized_molecules.csv not found in inputs/", flush=True)
        sys.exit(1)

    os.makedirs("outputs", exist_ok=True)
    output_path = "outputs/admet_predictions.csv"

    print(f"Running ADMET-AI predictions (smiles_column={smiles_column})...", flush=True)

    result = subprocess.run(
        [
            "admet_predict",
            "--data_path", input_path,
            "--save_path", output_path,
            "--smiles_column", smiles_column,
        ],
        capture_output=True,
        text=True,
    )

    if result.stdout:
        print(result.stdout, flush=True)
    if result.stderr:
        print(result.stderr, file=sys.stderr, flush=True)

    if result.returncode != 0:
        print(f"ERROR: admet_predict exited with code {result.returncode}", flush=True)
        sys.exit(result.returncode)

    if not os.path.exists(output_path):
        print("ERROR: admet_predict did not produce output file", flush=True)
        sys.exit(1)

    # Count results
    with open(output_path) as f:
        line_count = sum(1 for _ in f) - 1  # subtract header
    print(f"ADMET predictions computed for {line_count} molecules.", flush=True)


if __name__ == "__main__":
    main()
