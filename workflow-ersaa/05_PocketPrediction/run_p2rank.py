#!/usr/bin/env python3

import subprocess
from pathlib import Path
import sys

def run_p2rank():
    workdir = Path("./")

    pdb_files = list(workdir.glob("*.pdb"))
    if not pdb_files:
        print("[Node5] ERROR: No .pdb file found in the working directory.")
        sys.exit(1)

    protein = pdb_files[0]
    protein_abs = protein.resolve()

    print(f"[Node5] Found protein file: {protein.name}")
    print(f"[Node5] Absolute protein path: {protein_abs}")

    prank_path = Path("/opt/p2rank/p2rank_2.5.1/prank")

    if not prank_path.exists():
        print(f"[Node5] ERROR: prank executable not found at: {prank_path}")
        sys.exit(1)

    # IMPORTANT: Use -o instead of --output_dir
    cmd = [
        str(prank_path),
        "predict",
        "-f", str(protein_abs),
        "-o", str(workdir.resolve())
    ]

    print("[Node5] Running:", " ".join(cmd))

    subprocess.run(cmd, check=True)

    print("[Node5] P2Rank prediction completed successfully.")

if __name__ == "__main__":
    run_p2rank()
