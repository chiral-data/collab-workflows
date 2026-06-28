"""Node 03 — Design Sequences.

Calls ProteinMPNN NIM per backbone to design binder sequences.
Redesigns binder chain only; drops native/WT row from mfasta output.
"""
from __future__ import annotations
import json
import pathlib

PARAMS    = json.loads(pathlib.Path("inputs/global_params.json").read_text())
BB_LIST   = json.loads(pathlib.Path("inputs/backbone_list.json").read_text())

SEQS_PER_BACKBONE = PARAMS["seqs_per_backbone"]
SAMPLING_TEMP     = PARAMS.get("proteinmpnn_sampling_temp", 0.1)
NIM_MODE          = PARAMS.get("nim_mode", "hosted")

OUT = pathlib.Path("outputs")

# TODO: implement nim_post, parse_mfasta, drop_native_row, design_sequences loop

def main() -> None:
    raise NotImplementedError("design_sequences.py is not yet implemented")

if __name__ == "__main__":
    main()
