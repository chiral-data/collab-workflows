#!/usr/bin/env python3

import os

# Force OpenBabel 2.4.1 (avoid using conda's obabel 3.x)
os.environ["PATH"] = "/usr/local/bin:" + os.environ.get("PATH", "")
os.environ["LD_LIBRARY_PATH"] = "/usr/local/lib:" + os.environ.get("LD_LIBRARY_PATH", "")

import sys
import subprocess
from pathlib import Path

def convert_ligands(input_dir: Path, output_dir: Path):
    input_dir = Path("./")
    output_dir = Path("./")

    ligands = [f for f in input_dir.iterdir() if f.suffix.lower() == ".sdf"]
    print(f"[Node4] Found {len(ligands)} ligands (.sdf)")

    if not ligands:
        print("[Node4] ERROR: No SDF ligand files found.")
        sys.exit(1)

    for lig in ligands:
        lig_name = lig.stem
        pdb_temp = output_dir / f"{lig_name}.pdb"
        pdbqt_out = output_dir / f"{lig_name}.pdbqt"

        print(f"[Node4] Step 1: Converting {lig.name} → {pdb_temp.name} (OpenBabel)...")
        result = subprocess.run(
            ["obabel", str(lig), "-O", str(pdb_temp)],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        if result.returncode != 0:
            print(f"[Node4] ERROR: OpenBabel failed for {lig.name}")
            print(result.stderr.decode())
            sys.exit(1)

        print(f"[Node4] Step 2: Converting {pdb_temp.name} → {pdbqt_out.name} (prepare_ligand4.py)...")
        cmd = [
            "/opt/mgltools_install/bin/pythonsh",
            "/opt/mgltools_install/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_ligand4.py",
            "-l", str(pdb_temp),
            "-o", str(pdbqt_out),
        ]
        subprocess.run(cmd, check=True)

if __name__ == "__main__":
    try:
        input_dir = Path(sys.argv[1])
        output_dir = Path(sys.argv[2])
    except Exception:
        input_dir = Path("./")
        output_dir = Path("./")

    convert_ligands(input_dir, output_dir)
