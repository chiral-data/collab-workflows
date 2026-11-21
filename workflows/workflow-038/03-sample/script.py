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

# Node 02 writes <chain>.a3m for the independent alignment and <chain>.paired.a3m for the
# positionally paired one a heteromer also gets. The two must never be mixed in one mapping.
PAIRED_SUFFIX = ".paired.a3m"


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


def advise_batch_size(n_tokens, num_samples, batch_size):
    """Print upstream's own batch-size advice without acting on it.

    bioemu2's `warn_batch_size` logs this estimate of the largest batch that fits:

        max(1, int(0.8 * total_memory / (0.3e6 * num_tokens) - num_tokens / 75))

    It only warns, and that warning goes to a logger this node does not configure, so it is
    recomputed and printed here. It is deliberately **not** applied: batch_size changes which
    samples are drawn, so deriving it from the GPU would make a run unreproducible on a
    different card while reporting the same seed. The number a run used is the one in the
    manifest, and the advice is advice.
    """
    if not isinstance(n_tokens, int) or n_tokens <= 0:
        print(
            "Batch-size advice unavailable: node 01's token count did not reach this node.",
            flush=True,
        )
        return
    effective = batch_size if batch_size is not None else num_samples
    try:
        import torch

        if not torch.cuda.is_available():
            return
        total_memory = torch.cuda.get_device_properties(0).total_memory
    except Exception as exc:
        print(f"Batch-size advice unavailable: {exc}", flush=True)
        return

    estimate = max(1, int(0.8 * total_memory / (0.3e6 * n_tokens) - n_tokens / 75))
    print(
        f"Batch size in use: {effective} "
        f"({'explicit' if batch_size is not None else 'the BioEmu2 default, num_samples'}). "
        f"Upstream's estimate of the largest batch that fits {n_tokens} token(s) in "
        f"{round(total_memory / 1024 ** 3, 1)} GiB is {estimate}.",
        flush=True,
    )
    if min(effective, num_samples) > estimate:
        print(
            f"WARNING: batch size {effective} is above that estimate, so this run may run out "
            f"of GPU memory. Set batch_size={estimate} or lower if it does.",
            flush=True,
        )


def unpaired_a3m_paths(directory):
    """The per-chain unpaired A3M paths, excluding node 02's paired files.

    `<chain>.paired.a3m` has the stem `<chain>.paired`, which is not a chain id, and
    bioemu2's `load_a3m_msa` requires the mapping's keys to equal the system's protein chain
    ids exactly. A plain `*.a3m` glob would therefore make it raise.
    """
    return sorted(
        Path(path)
        for path in glob.glob(os.path.join(directory, "*.a3m"))
        if not path.endswith(PAIRED_SUFFIX)
    )


def paired_a3m_paths(directory):
    """Node 02's positionally paired A3M files, as {chain_id: Path}.

    The chain id is the basename with `.paired.a3m` removed. `Path.stem` is not enough:
    for `A.paired.a3m` it gives `A.paired`, and bioemu2 keys an MSA by chain id.
    """
    return {
        os.path.basename(path)[: -len(PAIRED_SUFFIX)]: Path(path)
        for path in sorted(glob.glob(os.path.join(directory, "*" + PAIRED_SUFFIX)))
    }


def staged_alignment(msa_summary):
    """Return the A3M paths Silva staged, refusing to continue if any are missing.

    Sampling without an alignment must never be a silent fallback. `EmuSession` builds its
    own MSA when handed `msa=None`, which means querying the public ColabFold server from
    inside a node that is supposed to be offline, and the run then looks exactly like a
    successful one: frames are produced, exit status is 0, and the manifest records the
    alignment node 02 fetched rather than the one the model actually used.

    Two independent bugs produced that outcome before this check existed. The declarations
    said `msa/*.a3m`, and Silva matches an `inputs` entry by **basename**, so the pattern
    could never match and the directory was never staged. `get_msa` was then handed `str`
    paths where it requires `Path`, so even a staged alignment would have been rejected.
    Whether an alignment is due is node 02's finding, not this node's guess: it is read from
    `protein_chains` in msa_summary.json (see CLAUDE.md, "don't re-derive upstream facts").
    """
    expected = list(msa_summary.get("protein_chains") or [])
    a3m_files = unpaired_a3m_paths("./inputs/msa")

    if not a3m_files and expected:
        fail(
            f"Node 02 reported protein chain(s) {expected} and so fetched an alignment, but "
            "no .a3m file was staged into ./inputs/msa. Sampling would continue by having "
            "BioEmu2 fetch its own alignment over the network, so this stops here instead. "
            "Check that node 02 declares `msa` as an output and this node declares it as an "
            "input: Silva matches an input entry by basename, so a pattern such as "
            "'msa/*.a3m' stages nothing and fails silently."
        )
    if not a3m_files:
        print(
            "No alignment to load, and node 02 reported no protein chain, so BioEmu2 needs "
            "none. Sampling without an MSA.",
            flush=True,
        )
    return a3m_files, expected


def load_alignment(system, a3m_files, msa_summary):
    """Load node 02's A3M files. Returns (msa, record); record states the path taken."""
    if not a3m_files:
        return None, {
            "loaded": False,
            "how": "no alignment was due (node 02 reported no protein chain)",
            "files": [],
            "paired": False,
        }

    from bioemu2 import MSA, get_msa

    # Keys are chain ids, since node 02 writes one <chain>.a3m per protein chain. Values are
    # Path, not str: get_msa calls Path methods on them and raises on a plain string.
    source = {path.stem: path for path in a3m_files}
    print(f"Loading {len(source)} alignment(s) from disk: {sorted(source)}", flush=True)
    try:
        # A dict maps one A3M per protein chain; a lone chain can be passed as a bare path.
        # Either way this goes through load_a3m_msa, which checks every A3M's query row
        # against the system's own sequence — the validation is the reason to come in through
        # get_msa rather than constructing an MSA by hand.
        if len(source) == 1:
            msa = get_msa(system, source=next(iter(source.values())))
        else:
            msa = get_msa(system, source=source)
    except Exception as exc:
        fail(
            f"Could not load the alignment from node 02: {exc}. The first record of each "
            "A3M must be the query and must match the system's sequence after gaps and "
            "insertions are removed."
        )

    files = [f"msa/{path.name}" for path in a3m_files]
    paired_paths = paired_a3m_paths("./inputs/msa")
    paired_expected = bool(msa_summary.get("paired"))

    if paired_expected and set(paired_paths) != set(source):
        fail(
            f"Node 02 recorded a paired alignment for chains {sorted(source)} but the staged "
            f"paired files cover {sorted(paired_paths)}. A heteromer's paired alignment is "
            "positionally matched across chains, so a partial set cannot be used and "
            "continuing would quietly drop the paired feature the model expects."
        )

    if paired_paths:
        # load_a3m_msa only ever returns an unpaired MSA, so the paired half has to be
        # attached here. MSA is a frozen dataclass whose __post_init__ requires the paired and
        # unpaired chain id sets to match, which is the check that this is well formed.
        try:
            msa = MSA(
                unpaired_a3m=dict(msa.unpaired_a3m),
                paired_a3m={
                    chain: path.read_text() for chain, path in sorted(paired_paths.items())
                },
            )
        except Exception as exc:
            fail(f"Could not attach node 02's paired alignment: {exc}")
        files += [f"msa/{path.name}" for path in sorted(paired_paths.values())]
        print(
            f"Attached the paired alignment for {len(paired_paths)} chain(s). BioEmu2's "
            "pairmixer uses it as a feature and substitutes zeros without it, so a heteromer "
            "sampled from the unpaired alignment alone is a different and worse run.",
            flush=True,
        )
    elif not paired_expected:
        print(
            "No paired alignment, and node 02 recorded none: a system with one unique protein "
            "sequence needs only the unpaired alignment.",
            flush=True,
        )

    return msa, {
        "loaded": True,
        "how": (
            "read from node 02's A3M files, paired and unpaired"
            if paired_paths
            else "read from node 02's A3M files (unpaired only; none was due)"
        ),
        "files": files,
        "paired": bool(paired_paths),
        "chains": sorted(source),
    }


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
    num_steps = int(os.environ.get("PARAM_NUM_STEPS", "0"))
    hf_ckpt_revision = os.environ.get(
        "PARAM_HF_CKPT_REVISION", "dd4b0c1b69984b68edcb6200a544213be4d75242"
    )

    if num_samples < 1:
        fail(f"num_samples must be at least 1, got {num_samples}")
    if batch_size < 0:
        fail(f"batch_size must be 0 (BioEmu2 default) or positive, got {batch_size}")
    if num_steps < 0:
        fail(f"num_steps must be 0 (the sampler preset's default) or positive, got {num_steps}")
    if mode == "structure" and num_samples > 1:
        print(
            f"NOTE: mode=structure with num_samples={num_samples}. Structure mode is for one "
            "accurate folded structure; for a range of conformations use mode=equilibrium.",
            flush=True,
        )

    # 0 and -1 are the Silva-side encodings of "let BioEmu2 decide".
    resolved_batch_size = batch_size if batch_size > 0 else None
    resolved_seed = seed if seed >= 0 else None
    resolved_num_steps = num_steps if num_steps > 0 else None

    system_path = "./inputs/system.yaml"
    if not os.path.exists(system_path):
        fail(f"System YAML not found: {system_path}")
    shutil.copyfile(system_path, "./outputs/system.yaml")

    # Pass the alignment through so node 04 can propagate without re-querying the MSA
    # server. Node 04 depends on this node, not on node 02.
    # Every A3M, paired files included, so node 04 propagates against the same alignment.
    os.makedirs("./outputs/msa", exist_ok=True)
    for a3m in sorted(glob.glob("./inputs/msa/*.a3m")):
        shutil.copyfile(a3m, os.path.join("./outputs/msa", os.path.basename(a3m)))

    msa_summary = read_json("./inputs/msa_summary.json", "MSA summary")
    inventory = read_json("./inputs/system_inventory.json", "system inventory")
    if os.path.exists("./inputs/system_inventory.json"):
        shutil.copyfile(
            "./inputs/system_inventory.json", "./outputs/system_inventory.json"
        )
    else:
        # Keep the declared output present rather than letting the node fail on a missing
        # file it only forwards. The empty marker says which it is.
        with open("./outputs/system_inventory.json", "w") as handle:
            json.dump({"unavailable": "node 01's inventory was not staged into node 03"}, handle)

    # Checked before the dry-run exit, so `--dry-run` is a real test of the staging: the
    # defect this guards against costs nothing to reproduce and needs no GPU.
    a3m_files, expected_chains = staged_alignment(msa_summary)
    print(
        f"Alignment staged: {[p.name for p in a3m_files] or 'none'} "
        f"(node 02 reported protein chain(s) {expected_chains or 'none'})",
        flush=True,
    )

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
        # None means the preset's own default, which is 200 for every sampler bioemu2 ships.
        "num_steps": resolved_num_steps,
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
    advise_batch_size(inventory.get("n_tokens"), num_samples, resolved_batch_size)

    from bioemu2 import EmuSession, System, check_physicality, load_model, save_frames

    try:
        system = System.from_yaml(Path(system_path))
    except Exception as exc:
        fail(f"Could not load system YAML '{system_path}': {exc}")

    msa, msa_record = load_alignment(system, a3m_files, msa_summary)

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

    # A string preset is the public API and is what a default run passes. num_steps lives on
    # the sampler object rather than on sample(), so a non-default step count means building
    # the sampler here with sampler_from_preset, which applies the mode's step scale (1.5 for
    # structure, 1.0 for equilibrium) exactly as the string path would.
    sampler_arg = sampler
    if resolved_num_steps is not None:
        try:
            from bioemu2.sampling import sampler_from_preset

            sampler_arg = sampler_from_preset(
                sampler, mode=mode, num_steps=resolved_num_steps
            )
        except Exception as exc:
            fail(
                f"num_steps={resolved_num_steps} was requested but the sampler could not be "
                f"configured: {exc}. Set num_steps=0 to use the preset's default of 200 "
                "steps, which goes through bioemu2's public string-preset path."
            )
        print(f"Sampler configured with num_steps={resolved_num_steps}", flush=True)

    print(f"Sampling {num_samples} structure(s) in {mode} mode", flush=True)
    started = time.time()
    try:
        samples = session.sample(
            num_samples,
            batch_size=resolved_batch_size,
            sampler=sampler_arg,
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

    # `checked` is False in two unrelated situations, and the distinction has to survive into
    # the manifest. bioemu2's check_physicality returns passed=all-True, checked=False both
    # when mode is no_filter (the thresholds were switched off) and when the system has no
    # eligible protein atoms (there was nothing to check). Counting that all-True array as
    # passes reported pass_rate 1.0 for a run where nothing was ever examined.
    passed = getattr(physicality, "passed", None)
    checked = bool(getattr(physicality, "checked", False))
    if physicality_filter == "no_filter":
        not_checked_because = "physicality_filter=no_filter, so no thresholds were applied"
    elif not checked:
        not_checked_because = (
            "no eligible protein atoms: bioemu2 checks only multi-residue chains containing "
            "a standard amino acid"
        )
    else:
        not_checked_because = None

    n_passed = int(sum(bool(x) for x in passed)) if (checked and passed is not None) else None
    reasons = list(getattr(physicality, "reasons", None) or [])
    reason_counts = {}
    if checked:
        for reason in reasons:
            if reason:
                reason_counts[str(reason)] = reason_counts.get(str(reason), 0) + 1

    manifest = {
        "settings": settings,
        "bioemu2_version": package_version(),
        "n_frames": int(samples.n_frames),
        "n_atoms": int(samples.n_atoms),
        "physicality_checked": checked,
        "physicality_not_checked_because": not_checked_because,
        "n_passed": n_passed,
        "pass_rate": (
            round(n_passed / samples.n_frames, 3)
            if n_passed is not None and samples.n_frames
            else None
        ),
        # First failed check per rejected frame, counted. These are bioemu2's threshold names:
        # bad_bonds_backbone, bad_bonds_sidechain, bad_bonds_sulfur, within_residue_clashes,
        # valence_exceeded, isolated_atoms, same_chain_residue_clashes,
        # cross_chain_residue_clashes, duplicate_atoms, atoms_mismatch, bonds_mismatch, or
        # non_finite_atoms.
        "rejection_reasons": reason_counts,
        "msa_source": msa_summary.get("source"),
        "msa_min_depth": msa_summary.get("min_depth"),
        # What the model was actually given, as distinct from what node 02 fetched. These
        # differed silently while the alignment was failing to stage.
        "msa_loaded": msa_record["loaded"],
        "msa_how": msa_record["how"],
        "msa_files": msa_record["files"],
        "msa_paired": msa_record["paired"],
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
        if reason_counts:
            print(f"  first failed check per rejected frame: {reason_counts}", flush=True)
        if n_passed == 0:
            print(
                "WARNING: no frame passed the physicality checks. Rejected frames are kept "
                "by design, but consider a steering sampler (edm_gd or edm_fk).",
                flush=True,
            )
        elif n_passed * 2 < samples.n_frames:
            print(
                f"WARNING: more than half the frames ({samples.n_frames - n_passed}/"
                f"{samples.n_frames}) failed the physicality checks. Upstream treats this as "
                "a sign that the checkpoint or the system needs review.",
                flush=True,
            )
    else:
        print(f"Physicality: not checked — {not_checked_because}", flush=True)


if __name__ == "__main__":
    main()
