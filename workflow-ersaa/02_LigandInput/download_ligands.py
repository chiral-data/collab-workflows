#!/usr/bin/env python3

import os
import json
import urllib.request
from pathlib import Path

def download_pubchem_cid(cid, outdir: Path):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{cid}/SDF?record_type=3d"
    outfile = outdir / f"{cid}.sdf"

    try:
        urllib.request.urlretrieve(url, outfile)
        print(f"[Node2] Downloaded CID {cid}")
    except Exception as e:
        print(f"[Node2] Failed to download CID {cid}: {e}")


if __name__ == "__main__":
    cid_string = os.getenv("PARAM_LIGAND_CIDS")
    if not cid_string:
        raise ValueError("Missing global parameter: PARAM_LIGAND_CIDS")

    cids = [c.strip() for c in cid_string.split(",")]

    outdir = Path("/workspace/out/ligands")
    outdir.mkdir(parents=True, exist_ok=True)

    for cid in cids:
        download_pubchem_cid(cid, outdir)
