#!/usr/bin/env python3

import sys
import pandas as pd
from pathlib import Path

def choose_pocket(p2rank_dir: Path, output_file: Path):
    csv_path = p2rank_dir / "predictions.csv"
    df = pd.read_csv(csv_path)

    print(df[["rank", "center_x", "center_y", "center_z", "score"]])

    pocket = input("Enter pocket rank to use (default = 1): ") or "1"
    row = df[df["rank"] == int(pocket)].iloc[0]

    center = {
        "center_x": float(row["center_x"]),
        "center_y": float(row["center_y"]),
        "center_z": float(row["center_z"]),
    }

    import json
    output_file.write_text(json.dumps(center, indent=2))

    print(f"[Node6] Saved chosen pocket to {output_file}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python choose_pocket.py <P2RANK_DIR> <OUTPUT_JSON>")
        sys.exit(2)

    choose_pocket(Path(sys.argv[1]), Path(sys.argv[2]))
