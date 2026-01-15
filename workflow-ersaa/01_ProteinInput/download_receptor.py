#!/usr/bin/env python3

import os
from pathlib import Path
from Bio.PDB import PDBList

def download_receptor(pdb_id: str, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    pdb_id = pdb_id.strip().lower()

    print(f"[Node1] Downloading receptor: {pdb_id}")
    pdbl = PDBList()

    fn = pdbl.retrieve_pdb_file(pdb_id, pdir=str(outdir), file_format="pdb")
    saved = Path(fn)
    target = outdir / f"{pdb_id}.pdb"

    try:
        saved.rename(target)
    except Exception:
        import shutil
        shutil.copy(str(saved), str(target))

    print(f"[Node1] Saved receptor to {target}")
    return target


if __name__ == "__main__":
    pdb_id = os.getenv("PARAM_PROTEIN_ID")
    if not pdb_id:
        raise ValueError("Missing global parameter: PARAM_PROTEIN_ID")

    output_dir = Path("./")
    download_receptor(pdb_id, output_dir)
