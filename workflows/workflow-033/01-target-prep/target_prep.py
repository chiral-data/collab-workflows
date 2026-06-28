"""Node 01 — Target Prep.

Fetches target PDB from RCSB, extracts chain, remaps hotspot residues,
calls MSA-search NIM for target A3M.
"""
from __future__ import annotations
import json
import os
import pathlib
import urllib.request

PARAMS = json.loads(pathlib.Path("inputs/global_params.json").read_text())

PDB_ID     = PARAMS["target_pdb_id"].upper()
CHAIN      = PARAMS["target_chain"]
HOTSPOTS   = [h.strip() for h in PARAMS["hotspot_residues"].split(",")]
NIM_MODE   = PARAMS.get("nim_mode", "hosted")

OUT = pathlib.Path("outputs")

# TODO: implement fetch_pdb, extract_chain, remap_hotspots, call_msa_search_nim

def main() -> None:
    raise NotImplementedError("target_prep.py is not yet implemented")

if __name__ == "__main__":
    main()
