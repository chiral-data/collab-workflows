#!/usr/bin/env python3

import sys
import json
import pandas as pd
from pathlib import Path


def choose_pocket(workdir: Path, output_file: Path):
    """
    Select the highest-probability pocket from P2Rank predictions
    and write its center coordinates to JSON.
    """

    # P2Rank outputs predictions.csv in the working directory
    csv_path = workdir / "5kir.pdb_predictions.csv"

    if not csv_path.exists():
        print("[Node6] ERROR: 5kir.pdb_predictions.csv not found in working directory.")
        sys.exit(1)

    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip()

    required_cols = {"rank", "center_x", "center_y", "center_z"}
    if not required_cols.issubset(df.columns):
        print("[Node6] ERROR: predictions.csv missing required columns.")
        sys.exit(1)

    # Select pocket with rank = 1 (highest probability by P2Rank definition)
    row = df[df["rank"] == 1].iloc[0]

    center = {
        "center_x": float(row["center_x"]),
        "center_y": float(row["center_y"]),
        "center_z": float(row["center_z"]),
    }

    output_file.write_text(json.dumps(center, indent=2))

    print("[Node6] Selected pocket rank = 1")
    print(f"[Node6] Pocket center saved to {output_file}")


if __name__ == "__main__":
    # Silva uses the working directory "./"
    workdir = Path("./")

    # Default output file
    output_file = Path("selected_pocket.json")

    # Optional CLI override (does NOT break Silva)
    if len(sys.argv) >= 2:
        output_file = Path(sys.argv[1])

    choose_pocket(workdir, output_file)
