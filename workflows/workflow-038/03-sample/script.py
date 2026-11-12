#!/usr/bin/env python3
"""Node 03 — Sample Conformations (GPU).

Draw `num_samples` structures from the prepared system, check their physicality, and write
the ensemble as samples.xtc + topology.cif + physicality.json.

This node uses the **Python API, not the `bioemu2` CLI**, for one concrete reason: the CLI
exposes no flag for supplying a pre-computed alignment (`_common_options` in bioemu2/cli.py
has --system/--pdb-id/--checkpoint/--out/--batch-size/--seed/--sampler/... and nothing for an
MSA). Node 02 fetched the alignment precisely so that this node needs no network, and only
`EmuSession(..., msa=...)` can accept it. The Python path also lets this node write
`save_frames()` output straight into ./outputs/ with the exact declared filenames.

HF_HUB_OFFLINE=1 is set before any import: the checkpoint is baked into the image at a
pinned revision, and the alignment arrives from node 02, so nothing here should ever touch
the network. If the image were built wrong, this turns a silent multi-gigabyte download into
an immediate, legible error.

Supports --dry-run: resolve and print every setting, then exit 0 without loading the model
or touching the GPU.

Mode semantics (upstream): `structure` applies heuristic sharpening at step_scale 1.5 for one
accurate folded structure; `equilibrium` uses step_scale 1.0 unsharpened for a range of
conformations. Returned frames are independent samples, NOT a time series.
"""
import glob
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
import time

# Must precede the bioemu2 / huggingface_hub import to take effect.
os.environ.setdefault("HF_HUB_OFFLINE", "1")

MODES = {"equilibrium", "structure"}
SAMPLERS = {"edm", "edm_gd", "edm_fk", "edm_fkc"}
PHYSICALITY_MODES = {
    "high_quality",
    "high_quality_per_residue",
    "no_violations",
    "no_filter",
}
PRECISIONS = {"bf16-mixed", "32-true"}
MATMUL_PRECISIONS = {"highest", "high", "medium"}


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
    fail(
        "No GPU available (nvidia-smi not found, torch fallback failed). BioEmu2 requires "
        "an NVIDIA GPU with compute capability 7.0 or newer."
    )


def gpu_info():
    try:
        result = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=name,compute_cap,driver_version,memory.total",
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

        return version("bioemu2")
    except Exception:
        return "unknown"


def read_json(path, label):
    if not os.path.exists(path):
        return {}
    try:
        with open(path) as handle:
            return json.load(handle)
    except Exception as exc:
        print(f"WARNING: could not read {label}: {exc}", flush=True)
        return {}


def choice(name, value, allowed):
    if value not in allowed:
        fail(f"{name} must be one of {sorted(allowed)}, got '{value}'")
    return value


def load_alignment(system):
    """Load node 02's A3M files. Returns an MSA, or None when there is nothing to load."""
    a3m_files = sorted(glob.glob("./inputs/msa/*.a3m"))
    if not a3m_files:
        print(
            "No alignment files from node 02 — sampling without an explicit MSA. For a "
            "protein system this is unexpected; check node 02's msa_summary.json.",
            flush=True,
        )
        return None

    from bioemu2 import get_msa

    source = {os.path.splitext(os.path.basename(p))[0]: p for p in a3m_files}
    print(f"Loading {len(source)} alignment(s) from disk: {sorted(source)}", flush=True)
    try:
        # A dict maps one A3M per protein chain; a lone chain can be passed as a bare path.
        if len(source) == 1:
            return get_msa(system, source=next(iter(source.values())))
        return get_msa(system, source=source)
    except Exception as exc:
        fail(
            f"Could not load the alignment from node 02: {exc}. The first record of each "
            "A3M must be the query and must match the system's sequence after gaps and "
            "insertions are removed."
        )


def main():
    dry_run = "--dry-run" in sys.argv

    os.makedirs("./outputs", exist_ok=True)

    mode = choice("mode", os.environ.get("PARAM_MODE", "equilibrium"), MODES)
    num_samples = int(os.environ.get("PARAM_NUM_SAMPLES", "8"))
    batch_size = int(os.environ.get("PARAM_BATCH_SIZE", "0"))
    seed = int(os.environ.get("PARAM_SEED", "42"))
    sampler = choice("sampler", os.environ.get("PARAM_SAMPLER", "edm"), SAMPLERS)
    physicality_filter = choice(
        "physicality_filter",
        os.environ.get("PARAM_PHYSICALITY_FILTER", "high_quality"),
        PHYSICALITY_MODES,
    )
    precision = choice("precision", os.environ.get("PARAM_PRECISION", "bf16-mixed"), PRECISIONS)
    matmul_precision = choice(
        "matmul_precision",
        os.environ.get("PARAM_MATMUL_PRECISION", "high"),
        MATMUL_PRECISIONS,
    )
    compile_denoiser = os.environ.get("PARAM_COMPILE_DENOISER", "false").lower() == "true"
    hf_ckpt_revision = os.environ.get(
        "PARAM_HF_CKPT_REVISION", "dd4b0c1b69984b68edcb6200a544213be4d75242"
    )

    if num_samples < 1:
        fail(f"num_samples must be at least 1, got {num_samples}")
    if batch_size < 0:
        fail(f"batch_size must be 0 (BioEmu2 default) or positive, got {batch_size}")
    if mode == "structure" and num_samples > 1:
        print(
            f"NOTE: mode=structure with num_samples={num_samples}. Structure mode is for one "
            "accurate folded structure; for a range of conformations use mode=equilibrium.",
            flush=True,
        )

    # 0 and -1 are the Silva-side encodings of "let BioEmu2 decide".
    resolved_batch_size = batch_size if batch_size > 0 else None
    resolved_seed = seed if seed >= 0 else None

    system_path = "./inputs/system.yaml"
    if not os.path.exists(system_path):
        fail(f"System YAML not found: {system_path}")
    shutil.copyfile(system_path, "./outputs/system.yaml")

    # Pass the alignment through so node 04 can propagate without re-querying the MSA
    # server. Node 04 depends on this node, not on node 02.
    os.makedirs("./outputs/msa", exist_ok=True)
    for a3m in sorted(glob.glob("./inputs/msa/*.a3m")):
        shutil.copyfile(a3m, os.path.join("./outputs/msa", os.path.basename(a3m)))

    msa_summary = read_json("./inputs/msa_summary.json", "MSA summary")

    settings = {
        "mode": mode,
        "num_samples": num_samples,
        "batch_size": resolved_batch_size,
        "seed": resolved_seed,
        "sampler": sampler,
        "physicality_filter": physicality_filter,
        "precision": precision,
        "matmul_precision": matmul_precision,
        "compile_denoiser": compile_denoiser,
        "checkpoint": "bioemu2-preview",
        "hf_ckpt_revision": hf_ckpt_revision,
    }
    print(f"bioemu2 version: {package_version()}", flush=True)
    print(f"Checkpoint revision (pinned at image build): {hf_ckpt_revision}", flush=True)
    print(f"Resolved settings: {json.dumps(settings, indent=2)}", flush=True)

    if dry_run:
        print("Dry run requested — not loading the model or invoking the GPU.", flush=True)
        return

    check_gpu()
    print(f"GPU: {gpu_info()}", flush=True)

    from bioemu2 import EmuSession, System, check_physicality, load_model, save_frames

    try:
        system = System.from_yaml(Path(system_path))
    except Exception as exc:
        fail(f"Could not load system YAML '{system_path}': {exc}")

    msa = load_alignment(system)

    print("Loading the model", flush=True)
    try:
        model = load_model(
            device="cuda",
            precision=precision,
            compile_denoiser=compile_denoiser,
            matmul_precision=matmul_precision,
        )
    except Exception as exc:
        fail(
            f"Could not load the model: {exc}. The checkpoint is baked into the image under "
            "$HF_HOME and HF_HUB_OFFLINE=1 is set, so this is an image-build problem rather "
            "than a network one."
        )

    session = EmuSession(system, model, mode, msa=msa)

    print(f"Sampling {num_samples} structure(s) in {mode} mode", flush=True)
    started = time.time()
    try:
        samples = session.sample(
            num_samples,
            batch_size=resolved_batch_size,
            sampler=sampler,
            seed=resolved_seed,
        )
    except Exception as exc:
        message = str(exc)
        if "out of memory" in message.lower():
            fail(
                f"Out of GPU memory: {exc}. Lower batch_size (currently "
                f"{resolved_batch_size or num_samples}); if it already fails at 1, the system "
                "is too large for one GPU and context-parallel sampling, which needs four, "
                "is out of scope for this workflow."
            )
        fail(f"Sampling failed: {exc}")
    elapsed = round(time.time() - started, 1)

    if samples.n_frames != num_samples:
        print(
            f"WARNING: requested {num_samples} sample(s) but got {samples.n_frames}.",
            flush=True,
        )

    print(f"Checking physicality (mode={physicality_filter})", flush=True)
    try:
        physicality = check_physicality(samples, system=system, mode=physicality_filter)
    except Exception as exc:
        fail(f"Physicality check failed: {exc}")

    # save_frames writes exactly samples.xtc, topology.cif and physicality.json into the
    # directory given — which is why ./outputs is passed directly and nothing is copied.
    try:
        save_frames(samples, "./outputs", physicality=physicality)
    except Exception as exc:
        fail(f"Could not save the ensemble: {exc}")

    passed = getattr(physicality, "passed", None)
    checked = bool(getattr(physicality, "checked", False))
    n_passed = int(sum(bool(x) for x in passed)) if passed is not None else None

    manifest = {
        "settings": settings,
        "bioemu2_version": package_version(),
        "n_frames": int(samples.n_frames),
        "n_atoms": int(samples.n_atoms),
        "physicality_checked": checked,
        "n_passed": n_passed,
        "pass_rate": (
            round(n_passed / samples.n_frames, 3)
            if n_passed is not None and samples.n_frames
            else None
        ),
        "msa_source": msa_summary.get("source"),
        "msa_min_depth": msa_summary.get("min_depth"),
        "gpu": gpu_info(),
        "elapsed_seconds": elapsed,
        "frames_are_independent_samples": True,
        "caveat": (
            "Sampling frequency is not a calibrated Boltzmann probability, and these frames "
            "are independent samples rather than a trajectory. Frame counts are not free "
            "energies."
        ),
    }
    with open("./outputs/sampling_manifest.json", "w") as handle:
        json.dump(manifest, handle, indent=2)

    for name in ("samples.xtc", "topology.cif", "physicality.json"):
        path = os.path.join("./outputs", name)
        if not os.path.exists(path):
            fail(f"Declared output '{name}' was not written by save_frames().")

    print(
        f"Wrote {samples.n_frames} frame(s), {samples.n_atoms} atoms in {elapsed}s",
        flush=True,
    )
    if n_passed is not None and checked:
        print(f"Physicality: {n_passed}/{samples.n_frames} frame(s) passed", flush=True)
        if n_passed == 0:
            print(
                "WARNING: no frame passed the physicality checks. Rejected frames are kept "
                "by design, but consider a steering sampler (edm_gd or edm_fk).",
                flush=True,
            )
    else:
        print("Physicality: not checked (no eligible protein atoms)", flush=True)


if __name__ == "__main__":
    main()
