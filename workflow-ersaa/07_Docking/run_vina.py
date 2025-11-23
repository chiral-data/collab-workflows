#!/usr/bin/env python3

import json
import subprocess
import sys
import os
from pathlib import Path

def run_vina(receptor, ligand, center, box_size, outdir):
    outdir.mkdir(parents=True, exist_ok=True)

    x, y, z = center["center_x"], center["center_y"], center["center_z"]

    out = outdir / (ligand.stem + "_out.pdbqt")
    log = out.with_suffix(".log")

    cmd = [
        "vina",
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--center_x", str(x),
        "--center_y", str(y),
        "--center_z", str(z),
        "--size_x", str(box_size),
        "--size_y", str(box_size),
        "--size_z", str(box_size),
        "--out", str(out),
        "--log", str(log),
    ]

    subprocess.run(cmd, check=True)
    print(f"[Node7] Docked: {ligand.name}")


if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python run_vina.py <RECEPTOR_PDBQT> <LIGAND_PDBQT> <POCKET_JSON> <BOX_SIZE> <OUTPUT_DIR>")
        sys.exit(2)

    receptor = Path(sys.argv[1])
    ligand = Path(sys.argv[2])
    center = json.loads(Path(sys.argv[3]).read_text())
    box = int(sys.argv[4])
    outdir = Path(sys.argv[5])

    run_vina(receptor, ligand, center, box, outdir)

