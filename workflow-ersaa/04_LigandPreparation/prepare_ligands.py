#!/usr/bin/env python3

import os
from pathlib import Path

if __name__ == "__main__":
    lig_in = Path("/workspace/input/ligands")
    lig_out = Path("/workspace/out/ligands_prepared")
    lig_out.mkdir(parents=True, exist_ok=True)

    for sdf in lig_in.glob("*.sdf"):
        base = sdf.stem
        pdbqt = lig_out / f"{base}.pdbqt"

        os.system(
            f"/opt/mgltools/bin/prepare_ligand4.py "
            f"-l {sdf} -o {pdbqt}"
        )

    print("[Node4] All ligands prepared")
