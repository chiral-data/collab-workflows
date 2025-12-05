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
        print(f"[Node2] Downloaded CID {cid} → {outfile}")
    except Exception as e:
        print(f"[Node2] Failed to download CID {cid}: {e}")


if __name__ == "__main__":
    cid_string = os.getenv("PARAM_LIGAND_IDS")
    if not cid_string:
        raise ValueError("Missing global parameter: PARAM_LIGAND_IDS")

    try:
        cids = json.loads(cid_string)
    except Exception:
        raise ValueError(f"PARAM_LIGAND_IDS is not valid JSON: {cid_string}")

    outdir = Path("./")
    outdir.mkdir(parents=True, exist_ok=True)

    for cid in cids:
        download_pubchem_cid(cid, outdir)
