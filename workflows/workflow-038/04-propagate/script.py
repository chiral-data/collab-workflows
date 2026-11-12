#!/usr/bin/env python3
"""Node 04 — Propagate Ensemble (GPU).

Draw one time-lagged endpoint per conformation from node 03, using BioEmu2's implicit
transfer operator model, and write propagated_samples.xtc + propagated_topology.cif.

**This node always writes its declared outputs.** Propagation doubles GPU cost and most runs
will not want it, but Silva has no conditionally-skipped node, and a declared output that is
never written is a silent failure (SILVA_MIGRATION_GUIDE.md §9, pitfall 5). So when
`propagate` is false this node copies node 03's ensemble through unchanged and records
`{"enabled": false, "lag_time_ns": 0.0}` in the manifest — an identity propagation. Nodes 05
and 06 read `enabled` and must not present identity frames as time-lagged.

The starting ensemble is bound to the system's own topology (`system.topology`) rather than
parsed from topology.cif: the XTC holds coordinates only, and mdtraj reading PDBx/mmCIF is
the weak link in that chain. The System is the authoritative topology on both sides, so the
atom order is guaranteed to match by construction.

Propagation requires equilibrium embeddings — `EmuSession.propagate()` warns on other modes —
so the session is always created with mode="equilibrium" here regardless of node 03's mode.

Supports --dry-run.
"""
import glob
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
import time

os.environ.setdefault("HF_HUB_OFFLINE", "1")

SAMPLERS = {"edm", "edm_gd", "edm_fk", "edm_fkc"}
PRECISIONS = {"bf16-mixed", "32-true"}
MATMUL_PRECISIONS = {"highest", "high", "medium"}

# Upstream's recommended range for --lag-time-ns.
LAG_MIN_NS = 1.0
LAG_MAX_NS = 1000.0


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def choice(name, value, allowed):
    if value not in allowed:
        fail(f"{name} must be one of {sorted(allowed)}, got '{value}'")
    return value


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
            ["nvidia-smi", "--query-gpu=name,compute_cap,memory.total", "--format=csv,noheader"],
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


def write_manifest(payload):
    with open("./outputs/propagation_manifest.json", "w") as handle:
        json.dump(payload, handle, indent=2)


def passthrough(reason):
    """Identity propagation: copy node 03's ensemble into the declared output names."""
    shutil.copyfile("./inputs/samples.xtc", "./outputs/propagated_samples.xtc")
    shutil.copyfile("./inputs/topology.cif", "./outputs/propagated_topology.cif")
    write_manifest(
        {
            "enabled": False,
            "lag_time_ns": 0.0,
            "reason": reason,
            "note": (
                "Identity pass-through: these frames are node 03's equilibrium samples, "
                "copied unchanged. They are NOT time-lagged endpoints."
            ),
        }
    )
    print(f"Propagation disabled ({reason}) — ensemble passed through unchanged.", flush=True)


def main():
    dry_run = "--dry-run" in sys.argv

    os.makedirs("./outputs", exist_ok=True)

    enabled = os.environ.get("PARAM_PROPAGATE", "false").lower() == "true"
    lag_time_ns = float(os.environ.get("PARAM_LAG_TIME_NS", "1.0"))
    batch_size = int(os.environ.get("PARAM_BATCH_SIZE", "0"))
    seed = int(os.environ.get("PARAM_SEED", "42"))
    sampler = choice("sampler", os.environ.get("PARAM_SAMPLER", "edm"), SAMPLERS)
    precision = choice("precision", os.environ.get("PARAM_PRECISION", "bf16-mixed"), PRECISIONS)
    matmul_precision = choice(
        "matmul_precision",
        os.environ.get("PARAM_MATMUL_PRECISION", "high"),
        MATMUL_PRECISIONS,
    )

    for name in ("samples.xtc", "topology.cif", "system.yaml"):
        path = f"./inputs/{name}"
        if not os.path.exists(path):
            fail(f"Required input not found: {path}")

    if batch_size < 0:
        fail(f"batch_size must be 0 (BioEmu2 default) or positive, got {batch_size}")
    resolved_batch_size = batch_size if batch_size > 0 else None
    resolved_seed = seed if seed >= 0 else None

    if not enabled:
        if dry_run:
            print("Dry run: propagation is disabled; node would pass the ensemble through.",
                  flush=True)
            return
        passthrough("propagate parameter is false")
        return

    if lag_time_ns <= 0:
        fail(f"lag_time_ns must be positive, got {lag_time_ns}")
    if not (LAG_MIN_NS <= lag_time_ns <= LAG_MAX_NS):
        print(
            f"WARNING: lag_time_ns={lag_time_ns} is outside the range upstream recommends "
            f"({LAG_MIN_NS:g}-{LAG_MAX_NS:g} ns), which is the range the model was "
            "trained for.",
            flush=True,
        )

    settings = {
        "enabled": True,
        "lag_time_ns": lag_time_ns,
        "batch_size": resolved_batch_size,
        "seed": resolved_seed,
        "sampler": sampler,
        "precision": precision,
        "matmul_precision": matmul_precision,
        "session_mode": "equilibrium",
    }
    print(f"bioemu2 version: {package_version()}", flush=True)
    print(f"Resolved settings: {json.dumps(settings, indent=2)}", flush=True)

    if dry_run:
        print("Dry run requested — not loading the model or invoking the GPU.", flush=True)
        return

    check_gpu()
    print(f"GPU: {gpu_info()}", flush=True)

    import mdtraj

    from bioemu2 import EmuSession, System, get_msa, load_model, save_frames

    try:
        system = System.from_yaml(Path("./inputs/system.yaml"))
    except Exception as exc:
        fail(f"Could not load system YAML: {exc}")

    # Bind the XTC to the System's own topology: guaranteed atom-order match, no CIF parsing.
    try:
        x0 = mdtraj.load_xtc("./inputs/samples.xtc", top=system.topology)
    except Exception as exc:
        fail(f"Could not load the starting ensemble: {exc}")
    if x0.n_frames == 0:
        fail("The starting ensemble contains no frames.")
    print(f"Starting from {x0.n_frames} conformation(s)", flush=True)

    a3m_files = sorted(glob.glob("./inputs/msa/*.a3m"))
    msa = None
    if a3m_files:
        source = {os.path.splitext(os.path.basename(p))[0]: p for p in a3m_files}
        try:
            msa = (
                get_msa(system, source=next(iter(source.values())))
                if len(source) == 1
                else get_msa(system, source=source)
            )
        except Exception as exc:
            print(f"WARNING: could not reuse node 02's alignment: {exc}", flush=True)

    try:
        model = load_model(
            device="cuda",
            precision=precision,
            matmul_precision=matmul_precision,
        )
    except Exception as exc:
        fail(f"Could not load the model: {exc}")

    # Propagation needs equilibrium embeddings whatever node 03 sampled in.
    session = EmuSession(system, model, "equilibrium", msa=msa)

    print(f"Propagating by {lag_time_ns} ns", flush=True)
    started = time.time()
    try:
        endpoints = session.propagate(
            x0,
            lag_time_ns=lag_time_ns,
            batch_size=resolved_batch_size,
            sampler=sampler,
            seed=resolved_seed,
        )
    except Exception as exc:
        message = str(exc)
        if "out of memory" in message.lower():
            fail(
                f"Out of GPU memory: {exc}. Lower batch_size (currently "
                f"{resolved_batch_size or x0.n_frames})."
            )
        fail(
            f"Propagation failed: {exc}. If the error mentions checkpoint support, the "
            "pinned bioemu2-preview checkpoint may not support propagation."
        )
    elapsed = round(time.time() - started, 1)

    # save_frames writes samples.xtc / topology.cif, so write to a private directory and
    # rename into the declared propagated_* names.
    staging = "./_propagated"
    shutil.rmtree(staging, ignore_errors=True)
    try:
        save_frames(endpoints, staging)
    except Exception as exc:
        fail(f"Could not save the propagated ensemble: {exc}")

    shutil.copyfile(f"{staging}/samples.xtc", "./outputs/propagated_samples.xtc")
    shutil.copyfile(f"{staging}/topology.cif", "./outputs/propagated_topology.cif")
    shutil.rmtree(staging, ignore_errors=True)

    settings.update(
        {
            "bioemu2_version": package_version(),
            "n_start_frames": int(x0.n_frames),
            "n_endpoint_frames": int(endpoints.n_frames),
            "gpu": gpu_info(),
            "elapsed_seconds": elapsed,
            "note": (
                "One time-lagged endpoint per starting frame, in input order. This is a "
                "lagged endpoint per start point, not a continuous trajectory."
            ),
        }
    )
    write_manifest(settings)

    if endpoints.n_frames != x0.n_frames:
        print(
            f"WARNING: {x0.n_frames} starting frame(s) produced {endpoints.n_frames} "
            "endpoint(s); propagation should be one-to-one and in input order.",
            flush=True,
        )

    print(
        f"Wrote {endpoints.n_frames} endpoint(s) after {lag_time_ns} ns in {elapsed}s",
        flush=True,
    )


if __name__ == "__main__":
    main()
