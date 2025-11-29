#!/usr/bin/env python3

import os
import pandas as pd
from pathlib import Path

def extract_affinities(log_dir: Path):
    data = []

    for file in log_dir.iterdir():
        if file.suffix == ".log":
            with open(file) as f:
                for line in f:
                    if line.strip().startswith("1 "):
                        parts = line.split()
                        if len(parts) > 1:
                            data.append({
                                "Ligand": file.stem,
                                "Affinity (kcal/mol)": float(parts[1])
                            })
                        break

    if data:
        df = pd.DataFrame(data).sort_values("Affinity (kcal/mol)")
        out = log_dir / "binding_affinities.xlsx"
        df.to_excel(out, index=False)
        print(f"[Node8] Saved ranking to {out}")
        print(df.head())
    else:
        print("[Node8] No affinity values found.")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python rank_vina.py <LOG_DIR>")
        sys.exit(2)

    extract_affinities(Path(sys.argv[1]))
