#!/usr/bin/env python3
"""Node 01 — Acquire Target.

Resolve the design target: prefer a structure uploaded into ./inputs/ (.cif, .pdb or
.cif.zst), otherwise fetch PARAM_PDB_ID from RCSB. Normalise to mmCIF as target.cif and
emit target_inventory.json.

The inventory must carry BOTH numberings per chain. RFD4 contigs and sel() selectors use
the original mmCIF label_seq_id, while RCSB files are author-numbered, so a hotspot copied
from a paper ("A45") is usually author numbering and needs translating in node 02.
"""
import glob
import os
import sys

# TODO(scaffold): implement download + mmCIF normalisation + inventory extraction.


def main():
    os.makedirs("./outputs", exist_ok=True)

    pdb_id = os.environ.get("PARAM_PDB_ID", "")
    target_chain = os.environ.get("PARAM_TARGET_CHAIN", "A")

    uploads = sorted(
        glob.glob("./inputs/*.cif")
        + glob.glob("./inputs/*.pdb")
        + glob.glob("./inputs/*.cif.zst")
    )

    if not uploads and not pdb_id:
        print(
            "ERROR: No target structure. Upload a .cif/.pdb/.cif.zst file or set pdb_id.",
            flush=True,
        )
        sys.exit(1)

    source = uploads[0] if uploads else f"RCSB:{pdb_id}"
    print(f"Target source: {source}", flush=True)
    print(f"Target chain: {target_chain}", flush=True)

    print("ERROR: Node 01 is not implemented yet (scaffold only).", flush=True)
    sys.exit(1)


if __name__ == "__main__":
    main()
