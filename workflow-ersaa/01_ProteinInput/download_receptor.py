#!/usr/bin/env python3

import sys
import os
from pathlib import Path
from Bio.PDB import PDBList

def download_receptor(pdb_id: str, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    pdb_id = pdb_id.lower().strip()

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
    if len(sys.argv) < 3:
        print("Usage: python download_receptor.py <PDB_ID> <OUTPUT_DIR>")
        sys.exit(2)

    pdb_id = sys.argv[1]
    outdir = Path(sys.argv[2])
    download_receptor(pdb_id, outdir)
