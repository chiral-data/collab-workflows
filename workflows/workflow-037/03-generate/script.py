#!/usr/bin/env python3
"""Node 03 — Generate Designs (GPU).

Invoke rfd4-generate on the validated condition spec, with the checkpoint revision pinned
and the baked-in weights used offline, then flatten inference_samples/ into ./outputs/.

Three deliberate choices:
  * argv is a Python list passed to subprocess.run WITHOUT shell=True. The `~metrics`
    override is a shell-expansion hazard (upstream quotes it specifically to stop `~`
    becoming $HOME); bypassing the shell removes that class of bug entirely.
  * HF_HUB_OFFLINE=1, so an image with baked-in weights never silently hits the network.
  * generation runs in a scratch directory holding the spec and target.cif side by side,
    because the spec's `input:` is a bare relative filename.

Supports --dry-run: print the resolved argv and exit 0 without touching the GPU.

The override names are Hydra keys confirmed against RFD4-Proteina 0.0.3's
configs/inference.yaml and the upstream README: `inputs` (the spec path, not --config),
top-level `num_replicates`, `seed`, `out_dir`, `compile_lengths` and `hf_ckpt_revision`,
with the `metrics` config group removed by `~metrics`.
"""
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time

SPEC_NAME = "design_spec.yaml"
TARGET_NAME = "target.cif"


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def check_gpu():
    try:
        result = subprocess.run(["nvidia-smi"], capture_output=True)
        if result.returncode == 0:
            return
    except FileNotFoundError:
        pass
    try:
        import torch

        if torch.cuda.is_available():
            return
    except ImportError:
        pass
    fail("No GPU available (nvidia-smi not found, torch fallback failed)")


def gpu_info():
    try:
        result = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=name,driver_version,memory.total",
                "--format=csv,noheader",
            ],
            capture_output=True,
            text=True,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except FileNotFoundError:
        pass
    return "unknown"


def package_version():
    try:
        from importlib.metadata import version

        return version("rfd4-proteina")
    except Exception:
        return "unknown"


def build_argv(spec_name, out_dir, num_replicates, seed, disable_compile, revision):
    # Upstream's documented generation call:
    #   rfd4-generate inputs=binder.yaml num_replicates=2 out_dir=./runs/first '~metrics'
    # The shell quoting upstream notes around ~metrics is moot here: argv goes to
    # subprocess.run as a list, so nothing expands ~ to $HOME.
    argv = [
        "rfd4-generate",
        f"inputs={spec_name}",
        f"out_dir={out_dir}",
        f"num_replicates={num_replicates}",
        f"seed={seed}",
        f"hf_ckpt_revision={revision}",
        "~metrics",
    ]
    if disable_compile:
        argv.append("compile_lengths=null")
    return argv


def collect_designs(run_dir, dest_dir):
    """Flatten the generated CIFs out of inference_samples/ into ./outputs/.

    The real layout, confirmed by running generation and listing the tree, is

        <out_dir>/inference_samples/design_0/design_0_generated.cif
        <out_dir>/inference_metrics/per_example_metrics.csv

    so each design sits in its own subdirectory rather than directly under
    inference_samples/, and the glob has to descend. (inference_metrics/ is written even
    with ~metrics, holding no scores; it contains no CIF, so it is simply not collected.)
    """
    found = sorted(
        glob.glob(
            os.path.join(run_dir, "**", "inference_samples", "**", "*.cif"), recursive=True
        )
    )
    if not found:
        # Fall back to any CIF produced, so a layout change surfaces as a clear error
        # rather than an empty success.
        found = sorted(
            path
            for path in glob.glob(os.path.join(run_dir, "**", "*.cif"), recursive=True)
            if os.path.basename(path) != TARGET_NAME
        )
        if found:
            print(
                "WARNING: no inference_samples/ directory; collecting CIFs from the "
                "output tree instead. The upstream output layout may have changed.",
                flush=True,
            )
    collected = []
    for index, path in enumerate(found):
        base = os.path.basename(path)
        if not base.endswith("_generated.cif"):
            base = f"design_{index:03d}_generated.cif"
        dest = os.path.join(dest_dir, base)
        if os.path.exists(dest):
            # Flattening a tree can collide; keep both rather than overwrite one silently.
            parent = os.path.basename(os.path.dirname(path)) or f"dup{index}"
            base = f"{parent}_{base}"
            dest = os.path.join(dest_dir, base)
        shutil.copyfile(path, dest)
        collected.append(os.path.basename(dest))
    return collected


def main():
    dry_run = "--dry-run" in sys.argv

    os.makedirs("./outputs", exist_ok=True)

    num_replicates = int(os.environ.get("PARAM_NUM_REPLICATES", "4"))
    seed = int(os.environ.get("PARAM_SEED", "42"))
    disable_compile = os.environ.get("PARAM_DISABLE_COMPILE", "true").lower() == "true"
    hf_ckpt_revision = os.environ.get(
        "PARAM_HF_CKPT_REVISION", "a06e2509c7033b010fc43614e9fd840303588502"
    )

    if num_replicates < 1:
        fail(f"num_replicates must be at least 1, got {num_replicates}")

    spec_path = f"./inputs/{SPEC_NAME}"
    if not os.path.exists(spec_path):
        fail(f"Condition spec not found: {spec_path}")
    target_path = f"./inputs/{TARGET_NAME}"

    argv = build_argv(SPEC_NAME, "out", num_replicates, seed, disable_compile, hf_ckpt_revision)

    print(f"rfd4-proteina version: {package_version()}", flush=True)
    print(f"Checkpoint revision: {hf_ckpt_revision}", flush=True)
    print(f"Resolved argv: {argv}", flush=True)

    if dry_run:
        print("Dry run requested — not invoking the GPU.", flush=True)
        return

    check_gpu()
    print(f"GPU: {gpu_info()}", flush=True)

    run_dir = tempfile.mkdtemp(prefix="rfd4_")
    shutil.copyfile(spec_path, os.path.join(run_dir, SPEC_NAME))
    if os.path.exists(target_path):
        # The spec's `input:` is a bare filename, resolved relative to the run directory.
        shutil.copyfile(target_path, os.path.join(run_dir, TARGET_NAME))

    env = dict(os.environ)
    env["HF_HUB_OFFLINE"] = "1"  # weights are baked in; never reach for the network

    print(f"Generating {num_replicates} design(s) in {run_dir}", flush=True)
    started = time.time()
    result = subprocess.run(argv, cwd=run_dir, env=env)
    elapsed = round(time.time() - started, 1)

    if result.returncode != 0:
        # The scratch directory is left in place deliberately: it holds the spec and the
        # partial output tree, which is what a post-mortem needs.
        fail(
            f"rfd4-generate exited {result.returncode} after {elapsed}s (output above). "
            f"Run directory kept for inspection: {run_dir}"
        )

    designs = collect_designs(run_dir, "./outputs")
    if not designs:
        fail(
            f"rfd4-generate exited 0 but produced no structures. Run directory kept for "
            f"inspection: {run_dir}"
        )
    if len(designs) != num_replicates:
        print(
            f"WARNING: expected {num_replicates} design(s), collected {len(designs)}.",
            flush=True,
        )

    manifest = {
        "argv": argv,
        "num_replicates": num_replicates,
        "designs_collected": len(designs),
        "designs": designs,
        "seed": seed,
        "disable_compile": disable_compile,
        "hf_ckpt_revision": hf_ckpt_revision,
        "rfd4_proteina_version": package_version(),
        "gpu": gpu_info(),
        "elapsed_seconds": elapsed,
        "metrics_disabled": True,
        "securedna_screened": False,
        "biosafety_note": (
            "The RFD4 CLI/Python path performs no SecureDNA screening. These designs are "
            "unscreened and require biosecurity review before synthesis."
        ),
    }
    with open("./outputs/generation_manifest.json", "w") as handle:
        json.dump(manifest, handle, indent=2)

    shutil.rmtree(run_dir, ignore_errors=True)
    print(f"Collected {len(designs)} design(s) into ./outputs in {elapsed}s", flush=True)


if __name__ == "__main__":
    main()
