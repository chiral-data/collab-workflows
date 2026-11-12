#!/usr/bin/env python3
"""Node 03 — Generate Designs (GPU).

Invoke rfd4-generate on the validated condition spec, with the checkpoint revision pinned
and the baked-in weights used offline, then flatten inference_samples/ into ./outputs/.

Two deliberate choices:
  * argv is a Python list passed to subprocess.run WITHOUT shell=True. The `~metrics`
    override is a shell-expansion hazard (upstream quotes it specifically to stop `~`
    becoming $HOME); bypassing the shell removes that class of bug entirely.
  * HF_HUB_OFFLINE=1, so an image with baked-in weights never silently hits the network.

Supports --dry-run: print the resolved argv and exit 0 without touching the GPU.
"""
import os
import subprocess
import sys

# TODO(scaffold): build argv, run generation, flatten inference_samples/, write manifest.
# TODO(build): the package import name is unresolved upstream — the GitHub README shows
# `import rfd4_proteina`, the HF model card shows `from rfproteina.proteina import Proteina`.
# Resolve empirically inside the built image before any Python here imports the package.


def main():
    dry_run = "--dry-run" in sys.argv

    os.makedirs("./outputs", exist_ok=True)

    num_replicates = int(os.environ.get("PARAM_NUM_REPLICATES", "4"))
    seed = int(os.environ.get("PARAM_SEED", "42"))
    disable_compile = os.environ.get("PARAM_DISABLE_COMPILE", "true").lower() == "true"
    hf_ckpt_revision = os.environ.get(
        "PARAM_HF_CKPT_REVISION", "a06e2509c7033b010fc43614e9fd840303588502"
    )

    spec_path = "./inputs/design_spec.yaml"
    if not os.path.exists(spec_path):
        print(f"ERROR: Condition spec not found: {spec_path}", flush=True)
        sys.exit(1)

    argv = [
        "rfd4-generate",
        f"--config={spec_path}",
        f"inference.num_replicates={num_replicates}",
        f"inference.seed={seed}",
        "~metrics",
    ]
    if disable_compile:
        argv.append("compile_lengths=null")

    print(f"Checkpoint revision: {hf_ckpt_revision}", flush=True)
    print(f"Resolved argv: {argv}", flush=True)

    if dry_run:
        print("Dry run requested — not invoking the GPU.", flush=True)
        return

    print("ERROR: Node 03 is not implemented yet (scaffold only).", flush=True)
    sys.exit(1)


if __name__ == "__main__":
    main()
