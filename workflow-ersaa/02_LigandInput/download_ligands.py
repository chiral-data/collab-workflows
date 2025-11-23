#!/usr/bin/env python3

import sys
import os
import requests
from pathlib import Path

REQUEST_TIMEOUT = 30

def download_ligand(cid: str, outdir: Path, record_type="3d"):
    outdir.mkdir(parents=True, exist_ok=True)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF/?record_type={record_type}"

    outpath = outdir / f"{cid}_{record_type}.sdf"
    print(f"[Node2] Downloading ligand CID={cid} (record_type={record_type})")

    try:
        r = requests.get(url, timeout=REQUEST_TIMEOUT)
        r.raise_for_status()
        text = r.text

        outpath.write_text(text)
        print(f"[Node2] Saved {outpath}")
        return outpath

    except Exception as e:
        print(f"[Node2] ERROR for CID {cid}: {e}")
        return None


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python download_ligands.py <CID1,CID2,...> <OUTPUT_DIR> <RECORD_TYPE>")
        sys.exit(2)

    cids = sys.argv[1].split(",")
    outdir = Path(sys.argv[2])
    record = sys.argv[3]

    for cid in cids:
        download_ligand(cid.strip(), outdir, record)
