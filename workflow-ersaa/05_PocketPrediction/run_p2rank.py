#!/usr/bin/env python3

import os
import sys
import subprocess
from pathlib import Path

def run_p2rank(protein_pdb: Path, output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)

    print("[Node5] Running P2Rank...")
    subprocess.run(["prank", "predict", str(protein_pdb), "--output_dir", str(output_dir)], check=True)
    print("[Node5] P2Rank completed")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python run_p2rank.py <PROTEIN_PDB> <OUTPUT_DIR>")
        sys.exit(2)

    run_p2rank(Path(sys.argv[1]), Path(sys.argv[2]))
