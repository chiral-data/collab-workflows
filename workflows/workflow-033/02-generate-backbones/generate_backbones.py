"""Node 02 — Generate Backbones.

Calls RFdiffusion NIM N times with distinct seeds to generate de novo
binder backbone PDBs conditioned on the target hotspots.
"""
from __future__ import annotations
import json
import pathlib

PARAMS   = json.loads(pathlib.Path("inputs/global_params.json").read_text())
HOTSPOTS = json.loads(pathlib.Path("inputs/hotspots.json").read_text())
TARGET   = pathlib.Path("inputs/target.pdb").read_text()

N_BACKBONES      = PARAMS["n_backbones"]
BINDER_LEN_MIN   = PARAMS["binder_length_min"]
BINDER_LEN_MAX   = PARAMS["binder_length_max"]
DIFFUSION_STEPS  = PARAMS["diffusion_steps"]
NIM_MODE         = PARAMS.get("nim_mode", "hosted")

OUT = pathlib.Path("outputs")

# TODO: implement nim_post, build_contig, generate_backbone loop

def main() -> None:
    raise NotImplementedError("generate_backbones.py is not yet implemented")

if __name__ == "__main__":
    main()
