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

# Mirrors node 02's naming. See node 03 for why the two sets must not be mixed.
PAIRED_SUFFIX = ".paired.a3m"
# Samplers whose presets enable resampling, which propagation cannot use. `edm_fk` is
# GuidedFKSampler() with resampling_interval=3 and `edm_fkc` is FKCSampler() with
# ess_threshold=0.9, while `edm_gd` is GuidedFKSampler(resampling_interval=None) and so is
# fine. propagate_with_embeddings() requires every block of resampling_group_size consecutive
# x0 frames to be byte-identical, because resampling copies particles within a group and each
# group must share one conditioning frame. Node 03's samples are all different by
# construction, so these two can never work here.
RESAMPLING_SAMPLERS = {"edm_fk", "edm_fkc"}
PHYSICALITY_MODES = {
    "high_quality",
    "high_quality_per_residue",
    "no_violations",
    "no_filter",
}
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


def n_frames_in(xtc_path):
    """Frame count of an XTC, for a per-frame record written without loading the model."""
    try:
        import mdtraj

        with mdtraj.open(xtc_path) as handle:
            return len(handle)
    except Exception as exc:
        print(f"WARNING: could not count frames in {xtc_path}: {exc}", flush=True)
        return 0


def read_json(path, label):
    if not os.path.exists(path):
        print(f"WARNING: {label} not found at {path}", flush=True)
        return {}
    try:
        with open(path) as handle:
            return json.load(handle)
    except Exception as exc:
        print(f"WARNING: could not read {label}: {exc}", flush=True)
        return {}


def write_manifest(payload):
    with open("./outputs/propagation_manifest.json", "w") as handle:
        json.dump(payload, handle, indent=2)


def write_unavailable_physicality(n_frames, reason):
    """Record one `unavailable` verdict per frame.

    Every declared output has to exist, and an absent physicality file would read downstream
    as "not checked" — which is a real BioEmu2 state meaning the thresholds were disabled, not
    this one. `unavailable` is deliberately not in upstream's vocabulary of passed / rejected /
    not_checked, so it cannot be mistaken for any of them.
    """
    records = [
        {"frame": frame, "status": "unavailable", "reason": reason} for frame in range(n_frames)
    ]
    with open("./outputs/propagated_physicality.json", "w") as handle:
        json.dump(records, handle, indent=2)


def passthrough(reason):
    """Identity propagation: copy node 03's ensemble into the declared output names."""
    shutil.copyfile("./inputs/samples.xtc", "./outputs/propagated_samples.xtc")
    shutil.copyfile("./inputs/topology.cif", "./outputs/propagated_topology.cif")
    # These frames *are* node 03's frames, so its verdicts describe them exactly. Copying the
    # report forward keeps the declared output present and keeps the physicality column
    # populated for a pass-through run, rather than blanking it.
    if os.path.exists("./inputs/physicality.json"):
        shutil.copyfile("./inputs/physicality.json", "./outputs/propagated_physicality.json")
    else:
        write_unavailable_physicality(
            n_frames_in("./inputs/samples.xtc"),
            "node 03's physicality.json was not staged into node 04",
        )
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
    num_steps = int(os.environ.get("PARAM_NUM_STEPS", "0"))
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

    for name in ("samples.xtc", "topology.cif", "system.yaml"):
        path = f"./inputs/{name}"
        if not os.path.exists(path):
            fail(f"Required input not found: {path}")

    sampling_manifest = read_json("./inputs/sampling_manifest.json", "sampling manifest")

    if batch_size < 0:
        fail(f"batch_size must be 0 (BioEmu2 default) or positive, got {batch_size}")
    if num_steps < 0:
        fail(f"num_steps must be 0 (the sampler preset's default) or positive, got {num_steps}")
    resolved_batch_size = batch_size if batch_size > 0 else None
    resolved_seed = seed if seed >= 0 else None
    resolved_num_steps = num_steps if num_steps > 0 else None

    if not enabled:
        if dry_run:
            print("Dry run: propagation is disabled; node would pass the ensemble through.",
                  flush=True)
            return
        passthrough("propagate parameter is false")
        return

    if sampler in RESAMPLING_SAMPLERS:
        fail(
            f"sampler='{sampler}' cannot be used with propagation. Its preset enables "
            "resampling, and resampling copies particles within a group of consecutive "
            "samples, so bioemu2 requires every frame in a group to be identical. The "
            "starting frames here are node 03's samples, which are all different, so this "
            "would fail inside propagate_with_embeddings with a message about frames "
            "differing. Use sampler='edm', or 'edm_gd', whose preset sets "
            "resampling_interval=None and is supported."
        )

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
        "num_steps": resolved_num_steps,
        "physicality_filter": physicality_filter,
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

    from bioemu2 import (
        MSA,
        EmuSession,
        System,
        check_physicality,
        get_msa,
        load_model,
        save_frames,
    )

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

    # Whether an alignment is due is node 03's recorded finding, not a guess from whatever
    # happens to be on disk here. A warning is the wrong response to a missing one: handing
    # EmuSession msa=None makes it fetch an alignment itself from the public ColabFold
    # server, so the node would propagate against a different alignment from the one node 03
    # sampled with, over the network, in a node that sets HF_HUB_OFFLINE and claims to need
    # none. Nothing downstream could tell.
    a3m_files = sorted(
        Path(p)
        for p in glob.glob("./inputs/msa/*.a3m")
        if not p.endswith(PAIRED_SUFFIX)
    )
    paired_paths = {
        os.path.basename(p)[: -len(PAIRED_SUFFIX)]: Path(p)
        for p in sorted(glob.glob("./inputs/msa/*" + PAIRED_SUFFIX))
    }
    msa_due = bool(sampling_manifest.get("msa_loaded"))
    paired_due = bool(sampling_manifest.get("msa_paired"))
    msa = None
    msa_record = {"loaded": False, "how": "no alignment was due", "files": [], "paired": False}

    if msa_due and not a3m_files:
        fail(
            "Node 03 sampled with an alignment "
            f"({sampling_manifest.get('msa_files') or 'files not recorded'}) but none was "
            "staged into ./inputs/msa here. Propagating would have BioEmu2 fetch its own "
            "alignment over the network, giving endpoints that do not match the ensemble "
            "they started from. Check that node 03 declares `msa` as an output and this node "
            "declares it as an input, by directory name: Silva matches an input entry by "
            "basename, so 'msa/*.a3m' stages nothing."
        )
    if paired_due and set(paired_paths) != {path.stem for path in a3m_files}:
        fail(
            "Node 03 sampled with a paired alignment but the paired files staged here cover "
            f"{sorted(paired_paths)} against chains "
            f"{sorted(path.stem for path in a3m_files)}. Propagating from the unpaired "
            "alignment alone would condition the endpoints differently from the ensemble they "
            "start from, which nothing downstream could detect."
        )

    if a3m_files:
        # Values are Path: get_msa calls Path methods on its source and raises on a str.
        source = {path.stem: path for path in a3m_files}
        try:
            msa = (
                get_msa(system, source=next(iter(source.values())))
                if len(source) == 1
                else get_msa(system, source=source)
            )
            if paired_paths:
                # load_a3m_msa returns an unpaired MSA whatever is on disk, so the paired half
                # is attached here exactly as node 03 does it. The model treats a missing
                # paired feature as zeros, so skipping this silently changes the run.
                msa = MSA(
                    unpaired_a3m=dict(msa.unpaired_a3m),
                    paired_a3m={
                        chain: path.read_text()
                        for chain, path in sorted(paired_paths.items())
                    },
                )
        except Exception as exc:
            fail(
                f"Could not reuse node 03's alignment: {exc}. Continuing without it would "
                "silently propagate against a freshly fetched alignment instead."
            )
        msa_record = {
            "loaded": True,
            "how": (
                "reused node 03's A3M files, paired and unpaired"
                if paired_paths
                else "reused node 03's A3M files (unpaired only; none was due)"
            ),
            "files": (
                [f"msa/{path.name}" for path in a3m_files]
                + [f"msa/{path.name}" for path in sorted(paired_paths.values())]
            ),
            "paired": bool(paired_paths),
            "chains": sorted(source),
        }
    elif not msa_due:
        print(
            "No alignment to reuse, and node 03 recorded none, so BioEmu2 needs none.",
            flush=True,
        )

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

    # As in node 03: the string preset is the public path and num_steps lives on the sampler
    # object, so a non-default step count is the only reason to build one here. Propagation
    # always uses the equilibrium step scale of 1.0, which is what bioemu2 would pick anyway
    # and what it warns about if overridden.
    sampler_arg = sampler
    if resolved_num_steps is not None:
        try:
            from bioemu2.sampling import sampler_from_preset

            sampler_arg = sampler_from_preset(
                sampler, mode="equilibrium", num_steps=resolved_num_steps
            )
        except Exception as exc:
            fail(
                f"num_steps={resolved_num_steps} was requested but the sampler could not be "
                f"configured: {exc}. Set num_steps=0 for the preset's default of 200 steps."
            )
        print(f"Sampler configured with num_steps={resolved_num_steps}", flush=True)

    print(f"Propagating by {lag_time_ns} ns", flush=True)
    started = time.time()
    try:
        endpoints = session.propagate(
            x0,
            lag_time_ns=lag_time_ns,
            batch_size=resolved_batch_size,
            sampler=sampler_arg,
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

    # A propagated endpoint can be unphysical exactly as a sampled frame can, and upstream
    # checks every frame at the end of a multi-step rollout for that reason. Without this,
    # node 03's frames carried a verdict and node 04's carried none, so the report showed a
    # blank column for the half of the ensemble that had been through an extra model pass.
    # Frames are reported, never removed, which is the same contract as node 03.
    physicality = None
    try:
        physicality = check_physicality(endpoints, system=system, mode=physicality_filter)
        physicality.save("./outputs/propagated_physicality.json")
    except Exception as exc:
        print(
            f"WARNING: could not check the propagated endpoints' physicality: {exc}. The "
            "endpoints are still written; node 05 reports the verdict as unavailable rather "
            "than as a pass.",
            flush=True,
        )
        write_unavailable_physicality(int(endpoints.n_frames), f"check failed: {exc}")

    settings.update(
        {
            "bioemu2_version": package_version(),
            "n_start_frames": int(x0.n_frames),
            "n_endpoint_frames": int(endpoints.n_frames),
            "gpu": gpu_info(),
            "elapsed_seconds": elapsed,
            "msa_loaded": msa_record["loaded"],
            "msa_how": msa_record["how"],
            "msa_files": msa_record["files"],
            "msa_paired": msa_record["paired"],
            "physicality_checked": bool(
                physicality is not None and getattr(physicality, "checked", False)
            ),
            "n_passed": (
                int(sum(bool(x) for x in physicality.passed))
                if physicality is not None and getattr(physicality, "checked", False)
                else None
            ),
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
