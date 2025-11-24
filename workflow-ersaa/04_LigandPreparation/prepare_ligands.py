#!/usr/bin/env python3

import os
import sys
import subprocess
from pathlib import Path

def convert_ligands(input_dir: Path, output_dir: Path):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ligands = [f for f in input_dir.iterdir() if f.suffix.lower() == ".pdb"]
    print(f"[Node4] Found {len(ligands)} ligands")

    for lig in ligands:
        out = output_dir / lig.name.replace(".pdb", ".pdbqt")

        print(f"[Node4] Converting {lig.name} -> {out.name}")
        cmd = [
            "python2", "/opt/mgltools/bin/prepare_ligand4.py",
            "-l", str(lig),
            "-o", str(out),
        ]
        subprocess.run(cmd, check=True)

    print(f"[Node4] Finished ligand preparation")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python prepare_ligands.py <INPUT_DIR> <OUTPUT_DIR>")
        sys.exit(2)
    convert_ligands(Path(sys.argv[1]), Path(sys.argv[2]))
