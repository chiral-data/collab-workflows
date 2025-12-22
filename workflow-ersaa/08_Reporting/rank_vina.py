#!/usr/bin/env python3

import sys
from pathlib import Path
import pandas as pd


def extract_affinities(log_dir: Path):
    data = []

    for file in log_dir.glob("*.log"):
        with file.open() as f:
            for line in f:
                parts = line.split()

                # Vina pose lines start with an integer pose index
                if len(parts) >= 2 and parts[0].isdigit():
                    try:
                        affinity = float(parts[1])
                    except ValueError:
                        continue

                    data.append({
                        "Ligand": file.stem,
                        "Affinity (kcal/mol)": affinity
                    })
                    break # only best pose

    if not data:
        print("[Node8] No affinity values found.")
        sys.exit(1)

    df = pd.DataFrame(data).sort_values("Affinity (kcal/mol)")
    out = log_dir / "binding_affinities.xlsx"

    df.to_excel(out, index=False)

    print(f"[Node8] Saved ranking to {out}")
    print(df.head())


if __name__ == "__main__":
    # Default to current working directory (Silva workspace)
    log_dir = Path(sys.argv[1]) if len(sys.argv) > 1 else Path.cwd()

    if not log_dir.exists():
        print(f"[Node8] Directory not found: {log_dir}")
        sys.exit(2)
    extract_affinities(log_dir)
