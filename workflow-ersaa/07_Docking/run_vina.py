#!/usr/bin/env python3

import os
import json
import subprocess
from pathlib import Path


def run_vina(receptor, ligand, box, out_prefix, exhaust, modes, energy):
    cmd = [
        "vina",
        "--receptor", str(receptor),
        "--ligand", str(ligand),
        "--center_x", str(box["center_x"]),
        "--center_y", str(box["center_y"]),
        "--center_z", str(box["center_z"]),
        "--size_x", str(box["size_x"]),
        "--size_y", str(box["size_y"]),
        "--size_z", str(box["size_z"]),
        "--exhaustiveness", str(exhaust),
        "--num_modes", str(modes),
        "--energy_range", str(energy),
        "--out", f"{out_prefix}.pdbqt",
        "--log", f"{out_prefix}.log"
    ]

    subprocess.run(cmd, check=True)
    print(f"[✔] Docked {ligand.name} → {out_prefix}.pdbqt")


if __name__ == "__main__":
    receptor = Path("/workspace/input/receptor.pdbqt")
    ligands_dir = Path("/workspace/input/ligands_prepared")
    grids_file = Path("/workspace/input/grids.json")
    out_dir = Path("/workspace/out/vina_results")

    out_dir.mkdir(parents=True, exist_ok=True)

    # -------------------------------
    # Read job-specific parameters
    # -------------------------------
    exhaust = os.getenv("JOB_PARAM_EXHAUSTIVENESS", "8")
    modes = os.getenv("JOB_PARAM_NUM_MODES", "9")        # <-- updated
    energy = os.getenv("JOB_PARAM_ENERGY_RANGE", "4")

    # Load boxes
    with open(grids_file) as f:
        boxes = json.load(f)

    ligands = sorted(ligands_dir.glob("*.pdbqt"))

    print(f"[Node7] Docking {len(ligands)} ligands using {len(boxes)} grid boxes")
    print(f"[Node7] Parameters → exhaust={exhaust}, num_modes={modes}, energy_range={energy}")

    for lig in ligands:
        for i, box in enumerate(boxes):
            out_prefix = out_dir / f"{lig.stem}_p{i}"
            run_vina(receptor, lig, box, out_prefix, exhaust, modes, energy)

    print("[Node7] Docking complete")


