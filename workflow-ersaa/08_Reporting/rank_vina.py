#!/usr/bin/env python3

import os
import pandas as pd
from pathlib import Path

if __name__ == "__main__":
    results_dir = Path("/workspace/input/vina_results")
    outfile = Path("/workspace/out/results.xlsx")

    rows = []

    for log in results_dir.glob("*.log"):
        with open(log) as f:
            for line in f:
                if line.strip().startswith("1 "):  # best mode
                    parts = line.split()
                    affinity = float(parts[1])
                    rows.append({"Ligand": log.stem, "Affinity": affinity})
                    break

    df = pd.DataFrame(rows).sort_values("Affinity")
    df.to_excel(outfile, index=False)

    print("[Node8] Ranking complete")
