"""Node 04 — Co-fold and Score.

Shortlists top candidates by ProteinMPNN NLL, co-folds each binder–target
complex with Boltz2 NIM, extracts ipTM + binder pLDDT from mmCIF,
computes self-consistency RMSD, writes run manifest.
"""
from __future__ import annotations
import json
import pathlib
import sys
import time

PARAMS         = json.loads(pathlib.Path("inputs/global_params.json").read_text())
SEQ_MANIFEST   = json.loads(pathlib.Path("inputs/sequence_manifest.json").read_text())
TARGET_SEQ     = pathlib.Path("inputs/chain_seq.txt").read_text().strip()
TARGET_A3M     = pathlib.Path("inputs/target_a3m.txt").read_text()

SHORTLIST_N          = PARAMS["cofold_shortlist_n"]
RECYCLING_STEPS      = PARAMS.get("boltz2_recycling_steps", 3)
SAMPLING_STEPS       = PARAMS.get("boltz2_sampling_steps", 50)
NIM_MODE             = PARAMS.get("nim_mode", "hosted")

OUT = pathlib.Path("outputs")

# TODO: implement nim_post with retry/backoff, parse_plddt_from_cif,
#       ca_rmsd (uses pdb_utils + metrics), cofold_score loop, manifest write

def main() -> None:
    raise NotImplementedError("cofold_score.py is not yet implemented")

if __name__ == "__main__":
    main()
