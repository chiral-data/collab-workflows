#!/usr/bin/env python3
"""Node 02 — Build Condition Spec.

Render the task-specific RFD4 condition spec (input / contig / conditions, with
C_CRDGRP, C_SEQ, C_HOT, C_CRD and sel() selectors) from the workflow parameters, then
validate it on CPU by converting it to an annotated CIF before any GPU time is spent.

This node exists because sel() and contig errors are the most common RFD4 failure mode.
Upstream's scripts/write_condition_spec_cifs.py renders a spec to an annotated CIF without
running the model, which turns a 20-minute GPU failure into a 10-second CPU failure — and
spec_preview.cif is itself useful, since the user can open it and see exactly which atoms
are conditioned.
"""
import os
import sys

TASKS = ("binder", "motif_scaffold", "monomer", "binder_redesign")

# TODO(scaffold): add templates/<task>.yaml.j2, auth->label residue translation, and the
# write_condition_spec_cifs.py validation pass.


def main():
    os.makedirs("./outputs", exist_ok=True)

    task = os.environ.get("PARAM_TASK", "binder")
    target_chain = os.environ.get("PARAM_TARGET_CHAIN", "A")
    binder_length = os.environ.get("PARAM_BINDER_LENGTH", "50-120")
    contig = os.environ.get("PARAM_CONTIG", "")
    hotspot_residues = os.environ.get("PARAM_HOTSPOT_RESIDUES", "")
    motif_residues = os.environ.get("PARAM_MOTIF_RESIDUES", "")
    redesign_chain = os.environ.get("PARAM_REDESIGN_CHAIN", "B")
    residue_numbering = os.environ.get("PARAM_RESIDUE_NUMBERING", "auth")

    if task not in TASKS:
        print(f"ERROR: Unknown task '{task}'. Expected one of: {', '.join(TASKS)}", flush=True)
        sys.exit(1)

    if residue_numbering not in ("auth", "label"):
        print(
            f"ERROR: residue_numbering must be 'auth' or 'label', got '{residue_numbering}'",
            flush=True,
        )
        sys.exit(1)

    target_cif = "./inputs/target.cif"
    inventory = "./inputs/target_inventory.json"
    for path in (target_cif, inventory):
        if not os.path.exists(path):
            print(f"ERROR: Required input not found: {path}", flush=True)
            sys.exit(1)

    print(f"Task: {task} (numbering: {residue_numbering})", flush=True)
    print(f"Target chain: {target_chain}", flush=True)
    print(f"Contig: {contig or f'templated from binder_length={binder_length}'}", flush=True)
    print(f"Hotspots: {hotspot_residues or 'none'}", flush=True)
    print(f"Motif: {motif_residues or 'none'}", flush=True)
    print(f"Redesign chain: {redesign_chain}", flush=True)

    print("ERROR: Node 02 is not implemented yet (scaffold only).", flush=True)
    sys.exit(1)


if __name__ == "__main__":
    main()
