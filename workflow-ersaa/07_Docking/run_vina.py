#!/usr/bin/env python3

import json
import subprocess
from pathlib import Path

WORKDIR = Path(".").resolve()
RECEPTOR = (WORKDIR / "receptor.pdbqt").resolve()
POCKET_FILE = WORKDIR / "selected_pocket.json"
PARAMS_FILE = WORKDIR / "params.json"

# ---------------------------
# Load pocket center
# ---------------------------
with open(POCKET_FILE) as f:
    pocket = json.load(f)

center_x = float(pocket["center_x"])
center_y = float(pocket["center_y"])
center_z = float(pocket["center_z"])

# ---------------------------
# Fixed box size
# ---------------------------
size_x = size_y = size_z = 80

# ---------------------------
# Load params
# ---------------------------
with open(PARAMS_FILE) as f:
    params = json.load(f)

exhaustiveness = int(str(params["exhaustiveness"]).strip())
num_modes = int(str(params["num_modes"]).strip())
energy_range = float(str(params["energy_range"]).strip())

print(f"[DEBUG] exhaustiveness={exhaustiveness}, num_modes={num_modes}, energy_range={energy_range}")

# ---------------------------
# Collect ligands (no renaming)
# ---------------------------
ligands = [
    f for f in WORKDIR.glob("*.pdbqt")
    if f.name != RECEPTOR.name
]

if not ligands:
    raise RuntimeError("No ligand PDBQT files found.")

# ---------------------------
# Run Vina
# ---------------------------
for ligand in ligands:
    out_file = WORKDIR / f"{ligand.stem}_vina_out.pdbqt"

    cmd = [
        "vina",
        "--receptor", str(RECEPTOR),
        "--ligand", str(ligand),
        "--center_x", str(center_x),
        "--center_y", str(center_y),
        "--center_z", str(center_z),
        "--size_x", str(size_x),
        "--size_y", str(size_y),
        "--size_z", str(size_z),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", str(num_modes),
        "--energy_range", str(energy_range),
        "--out", str(out_file)
    ]

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )

    # write log
    log_file = WORKDIR / f"{ligand.stem}.log"
    with open(log_file, "w") as f:
        f.write(result.stdout)

    print(result.stdout)
    print(result.stderr)

    if result.returncode != 0:
        print(f"❌ Docking failed for {ligand.name}")
        break

print("\nNode 7 finished successfully ✅")
