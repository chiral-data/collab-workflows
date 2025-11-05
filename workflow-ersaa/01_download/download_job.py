#!/usr/bin/env python


import os
import sys
import time
from pathlib import Path

import requests

# Biopython import
try:
    from Bio.PDB import PDBList
except Exception:
    PDBList = None
    # We'll raise later if needed

OUTDIR = Path(os.environ.get("OUTDIR", "/workspace/out"))
OUTDIR.mkdir(parents=True, exist_ok=True)

PDB_ID = os.environ.get("PDB_ID", "").strip()
CIDS_ENV = os.environ.get("CIDS", "").strip()
PUBCHEM_RECORD_TYPE = os.environ.get("PUBCHEM_RECORD_TYPE", "3d").strip().lower()

# Basic settings
REQUEST_TIMEOUT = 30  # seconds


def download_pdb(pdb_id: str, outdir: Path) -> Path:
    if PDBList is None:
        raise RuntimeError(
            "Biopython is required but not installed. Install with: pip install biopython"
        )
    pdbl = PDBList()
    print(f"[protein] downloading PDB '{pdb_id}' ...")
    # Save into outdir; PDBList will create a file like pdbXXXX.ent, so we rename to pdb_id.pdb
    fn = pdbl.retrieve_pdb_file(pdb_id, pdir=str(outdir), file_format="pdb")
    saved = Path(fn)
    target = outdir / f"{pdb_id.lower()}.pdb"
    try:
        saved.rename(target)
    except Exception:
        import shutil

        shutil.copy(str(saved), str(target))
    print(f"[protein] saved to {target}")
    return target


def download_ligand_sdf(cid: str, outdir: Path, record_type="3d") -> Path:
    cid = cid.strip()
    if not cid:
        raise ValueError("Empty CID")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/record/SDF/?record_type={record_type}"
    outpath = outdir / f"{cid}_{record_type}.sdf"
    print(f"[ligand] fetching CID={cid} -> record_type={record_type}")
    try:
        r = requests.get(url, timeout=REQUEST_TIMEOUT)
        r.raise_for_status()
        text = r.text
        if len(text) < 200:
            # short response likely means PubChem has no 3D data; warn but still save
            print(f"[ligand] warning: small response for CID {cid} (len={len(text)}).")
        outpath.write_text(text)
        print(f"[ligand] saved {outpath}")
        return outpath
    except Exception as e:
        print(f"[ligand] ERROR downloading CID {cid}: {e}", file=sys.stderr)
        return None


def main():
    start = time.time()
    did_something = False

    PDB_ID = "1ALK"
    if PDB_ID:
        did_something = True
        try:
            download_pdb(PDB_ID, OUTDIR)
        except Exception as e:
            print(f"[protein] ERROR: {e}", file=sys.stderr)

    CIDS_ENV = "5381226"
    PUBCHEM_RECORD_TYPE = "2d"
    if CIDS_ENV:
        did_something = True
        cids = [c.strip() for c in CIDS_ENV.split(",") if c.strip()]
        print(
            f"[ligand] will download {len(cids)} CIDs (record_type={PUBCHEM_RECORD_TYPE})"
        )
        for cid in cids:
            download_ligand_sdf(cid, OUTDIR, record_type=PUBCHEM_RECORD_TYPE)

    if not did_something:
        print("No PDB_ID or CIDS provided. Nothing to do.", file=sys.stderr)
        sys.exit(2)

    elapsed = time.time() - start
    print(f"Done. elapsed: {elapsed:.1f}s. outputs in: {OUTDIR}")
    sys.exit(0)


if __name__ == "__main__":
    main()
