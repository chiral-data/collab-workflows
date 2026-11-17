#!/usr/bin/env python3
"""Node 03 — Generate Designs (GPU).

Invoke rfd4-generate on the validated condition spec, with the checkpoint revision pinned
and the baked-in weights used offline, then flatten inference_samples/ into ./outputs/.

Three deliberate choices:
  * argv is a Python list passed to subprocess.Popen WITHOUT shell=True. The `~metrics`
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
import signal
import subprocess
import sys
import tempfile
import threading
import time
from collections import deque

SPEC_NAME = "design_spec.yaml"
TARGET_NAME = "target.cif"

# The allocator logs these and carries on. None of them stops the process or changes its
# exit status, so they are only ever visible by reading the output.
OOM_MARKERS = (
    "CUDACachingAllocator",
    "out of memory",
    "OutOfMemoryError",
    "CUDA out of memory",
)
TAIL_LINES = 40


def preserved_chain_ids(spec_path):
    """Target chain ids of the designed system, read off the contig node 02 wrote.

    Not the chain ids in target.cif, which answer a different question. That file is the
    whole normalised entry: every polymer chain the contig did not preserve, plus one
    label_asym_id per ligand and solvent group. On 4ZQK it carries A, B, C, D and E where
    the designed system holds one chain -- and the generator draws the designed chain's
    id from the same alphabet, so the design comes back as B, which is already in that
    list. Node 04 is then left with no candidate at all, falls back to the first chain in
    the file, and the manifest written to spare it that guess has made it no wiser.

    The contig is where node 02 said which chains the design is built against, so that is
    the fact to carry forward rather than re-derive. Returns [] for a spec with no contig,
    which is the monomer passthrough: there is no target, and the design is the only
    chain in the output.
    """
    import re

    try:
        import yaml

        with open(spec_path) as handle:
            spec = yaml.safe_load(handle) or {}
        contig = (spec.get("design") or {}).get("contig")
    except Exception as exc:
        print(f"WARNING: could not read the contig from '{spec_path}': {exc}", flush=True)
        return None
    if not contig:
        return []
    chains = []
    # Segments are "<length>", "<length>-<length>", "/0" and "<chain><start>-<end>".
    # Only the last names a chain, and it is the only one that starts with a letter.
    for token in str(contig).split(","):
        match = re.match(r"^([A-Za-z]+)-?\d", token.strip())
        if match and match.group(1) not in chains:
            chains.append(match.group(1))
    if not chains:
        print(
            f"WARNING: no preserved chain found in contig '{contig}'; falling back to "
            "the chain ids in target.cif.",
            flush=True,
        )
        return None
    return chains


def target_chain_ids(path):
    """Every chain id in the target file -- the fallback when the contig cannot be read.

    Superseded by preserved_chain_ids(), which reads node 02's decision instead of
    re-deriving it. Kept because a spec this node cannot parse should still leave node 04
    something rather than nothing, and because for a single-chain target the two agree.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    try:
        cif = MMCIF2Dict(path)
    except Exception as exc:
        print(f"WARNING: could not read chain ids from '{path}': {exc}", flush=True)
        return None
    asym = cif.get("_atom_site.label_asym_id")
    if not asym:
        print(f"WARNING: '{path}' has no _atom_site records; no target chains recorded.",
              flush=True)
        return None
    seen = []
    for chain in asym:
        if chain not in seen:
            seen.append(chain)
    return seen


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def int_param(env_name, display_name, default):
    """Read a whole-number parameter, failing the way every other bad parameter does.

    Casting inline at the point of use turns a typo into a raw ValueError traceback, which
    reads as a crash in the workflow rather than as the user error it is.
    """
    raw = os.environ.get(env_name)
    if raw is None or raw.strip() == "":
        return default
    try:
        return int(raw)
    except ValueError:
        fail(f"{display_name} must be a whole number, got '{raw}'.")


def bool_param(env_name, display_name, default):
    """Read a true/false parameter, accepting what a person would reasonably type.

    Comparing against the literal "true" quietly treats "1", "yes" and "True " as false,
    which reverses the default rather than reporting anything. Silva emits "true"/"false"
    for a boolean parameter, so this mostly matters when the variable is set by hand.
    """
    raw = os.environ.get(env_name)
    if raw is None or raw.strip() == "":
        return default
    value = raw.strip().lower()
    if value in ("true", "1", "yes", "on"):
        return True
    if value in ("false", "0", "no", "off"):
        return False
    fail(f"{display_name} must be true or false, got '{raw}'.")


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


def build_argv(
    spec_name, out_dir, num_replicates, seed, disable_compile, revision, batch_size
):
    # Upstream's documented generation call:
    #   rfd4-generate inputs=binder.yaml num_replicates=2 out_dir=./runs/first '~metrics'
    # The shell quoting upstream notes around ~metrics is moot here: argv goes to
    # Popen as a list, so nothing expands ~ to $HOME.
    argv = [
        "rfd4-generate",
        f"inputs={spec_name}",
        f"out_dir={out_dir}",
        f"num_replicates={num_replicates}",
        f"seed={seed}",
        f"hf_ckpt_revision={revision}",
        f"batch_size={batch_size}",
        "~metrics",
    ]
    if disable_compile:
        argv.append("compile_lengths=null")
    return argv


def kill_tree(proc):
    """End the whole generation process group, not just the process we started."""
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(os.getpgid(proc.pid), sig)
        except (ProcessLookupError, PermissionError, OSError):
            proc.kill() if sig == signal.SIGKILL else proc.terminate()
        try:
            proc.wait(timeout=20)
            return
        except subprocess.TimeoutExpired:
            continue


def run_generation(argv, cwd, env, timeout_s):
    """Run rfd4-generate under a wall-clock budget, streaming and retaining its output.

    Two reasons this is not subprocess.run(argv, timeout=...):

    The budget has to be enforceable. subprocess.run's timeout signals the direct child
    only, while generation spawns workers that go on holding the GPU. The child is given
    its own process group here so the whole tree can be ended.

    The output has to be read. When the GPU runs out of memory the CUDA caching allocator
    prints a warning and continues: the process does not raise, does not exit, and does
    not change its status. It simply stops making progress while holding the device.
    Without keeping the output there is no way to tell that apart from a slow job, and
    the failure that reaches the user is a timeout with no cause attached.

    Returns (returncode, timed_out, saw_oom, tail).
    """
    proc = subprocess.Popen(
        argv,
        cwd=cwd,
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        start_new_session=True,
    )
    tail = deque(maxlen=TAIL_LINES)
    # The reader is a daemon and the join below is bounded, so the snapshot of the tail
    # can in principle run while a line is still being appended: a grandchild that
    # survives the kill keeps the pipe open and the thread alive. Copying a deque that
    # is being mutated raises, and the traceback would land exactly where the timeout
    # explanation belongs. Not reproducible in 200k tries -- copying 40 lines fits
    # inside one GIL slice -- but cheap to make impossible rather than merely unlikely.
    tail_lock = threading.Lock()
    state = {"oom": False}

    def pump():
        for line in proc.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            with tail_lock:
                tail.append(line.rstrip("\n"))
            if any(marker in line for marker in OOM_MARKERS):
                state["oom"] = True

    reader = threading.Thread(target=pump, daemon=True)
    reader.start()

    timed_out = False
    try:
        returncode = proc.wait(timeout=timeout_s)
    except subprocess.TimeoutExpired:
        timed_out = True
        returncode = None
        kill_tree(proc)
    reader.join(timeout=10)
    with tail_lock:
        captured = list(tail)
    return returncode, timed_out, state["oom"], captured


def oom_advice(saw_oom):
    if not saw_oom:
        return ""
    return (
        "\nThe CUDA allocator reported an out-of-memory condition. It logs that as a "
        "warning and carries on, so the job keeps the GPU without making progress and "
        "never fails on its own. Reduce binder_length or num_replicates, or run on a "
        "GPU with more memory."
    )


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
            # Only upstream's *_generated.cif files are designs. Anything else in that
            # tree (a copy of the resolved input, an intermediate) is skipped rather than
            # renamed into a design, which would put it in the metrics table as one.
            print(
                f"WARNING: ignoring '{base}' under inference_samples/: not a "
                "*_generated.cif file, so not treated as a design.",
                flush=True,
            )
            continue
        dest = os.path.join(dest_dir, base)
        if os.path.exists(dest):
            # Flattening a tree can collide; keep both rather than overwrite one silently.
            parent = os.path.basename(os.path.dirname(path)) or "dup"
            dest = os.path.join(dest_dir, f"{parent}_{base}")
            suffix = 2
            while os.path.exists(dest):
                dest = os.path.join(dest_dir, f"{parent}_{suffix}_{base}")
                suffix += 1
        shutil.copyfile(path, dest)
        collected.append(os.path.basename(dest))
    return collected


def baked_revisions():
    """Revisions of nvidia/RFD4-Proteina present in the image's HF cache."""
    pattern = os.path.join(
        os.environ.get("HF_HOME", "/models"),
        "hub",
        "models--nvidia--RFD4-Proteina",
        "snapshots",
        "*",
    )
    return sorted(os.path.basename(p) for p in glob.glob(pattern) if os.path.isdir(p))


def main():
    dry_run = "--dry-run" in sys.argv

    os.makedirs("./outputs", exist_ok=True)

    num_replicates = int_param("PARAM_NUM_REPLICATES", "num_replicates", 4)
    seed = int_param("PARAM_SEED", "seed", 42)
    disable_compile = bool_param("PARAM_DISABLE_COMPILE", "disable_compile", True)
    hf_ckpt_revision = os.environ.get(
        "PARAM_HF_CKPT_REVISION", "a06e2509c7033b010fc43614e9fd840303588502"
    )
    timeout_s = int_param("PARAM_GENERATE_TIMEOUT_S", "generate_timeout_s", 3600)
    batch_size = int_param("PARAM_BATCH_SIZE", "batch_size", 1)

    if num_replicates < 1:
        fail(f"num_replicates must be at least 1, got {num_replicates}")
    if batch_size < 1:
        fail(f"batch_size must be at least 1, got {batch_size}")
    if timeout_s < 1:
        fail(f"generate_timeout_s must be at least 1 second, got {timeout_s}")

    # The weights are baked in at one revision and HF_HUB_OFFLINE=1 keeps this node off the
    # network, so any other revision can only fail deep inside rfd4-generate as a cache
    # miss. Say so here instead.
    available = baked_revisions()
    if available and hf_ckpt_revision not in available:
        fail(
            f"hf_ckpt_revision '{hf_ckpt_revision}' is not baked into this image "
            f"(it has: {', '.join(available)}). This node runs offline, so another "
            "revision cannot be fetched — rebuild the image to change it."
        )

    spec_path = f"./inputs/{SPEC_NAME}"
    if not os.path.exists(spec_path):
        fail(f"Condition spec not found: {spec_path}")
    target_path = f"./inputs/{TARGET_NAME}"

    argv = build_argv(
        SPEC_NAME,
        "out",
        num_replicates,
        seed,
        disable_compile,
        hf_ckpt_revision,
        batch_size,
    )

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

    print(
        f"Generating {num_replicates} design(s) in {run_dir} "
        f"(budget {timeout_s}s)",
        flush=True,
    )
    started = time.time()
    returncode, timed_out, saw_oom, tail = run_generation(argv, run_dir, env, timeout_s)
    elapsed = round(time.time() - started, 1)

    # The scratch directory is left in place on every failure path below: it holds the
    # spec and the partial output tree, which is what a post-mortem needs. It only
    # survives long enough to be read when this node is run directly, as the README's
    # docker invocation does; under Silva the container exits with the job and takes the
    # directory with it, so the message says where it would be rather than promising it
    # is still there. The retained tail above is what reaches a Silva user.
    kept = f"\nRun directory, if this container outlives the job: {run_dir}"
    if timed_out:
        fail(
            f"rfd4-generate exceeded its {timeout_s}s budget and was killed after "
            f"{elapsed}s. Raise generate_timeout_s if the job was merely slow."
            + oom_advice(saw_oom)
            + f"\nLast output:\n  " + "\n  ".join(tail[-12:])
            + kept
        )
    if returncode != 0:
        fail(
            f"rfd4-generate exited {returncode} after {elapsed}s (output above)."
            + oom_advice(saw_oom)
            + kept
        )
    if saw_oom:
        print(
            "WARNING: the CUDA allocator reported an out-of-memory condition during this "
            "run. It completed anyway, but it was close to the limit of this GPU.",
            flush=True,
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

    target_chains = preserved_chain_ids(spec_path)
    if target_chains is None:
        target_chains = (
            target_chain_ids(target_path) if os.path.exists(target_path) else None
        )

    manifest = {
        "argv": argv,
        "num_replicates": num_replicates,
        "designs_collected": len(designs),
        "designs": designs,
        # The chains of the designed system that came in as target, taken from the
        # contig. Node 04 needs this to name the designed chain without reading anything
        # into chain order: the design is the chain that was not already there. The
        # contig does put the designed segment first, but that is a fact about node 02,
        # and a node two steps downstream should not be inferring it. Reading the whole
        # target file instead lists chains the contig dropped and the solvent's asym
        # ids, which can collide with the id the generator gives the design.
        "target_chains": target_chains,
        "seed": seed,
        "disable_compile": disable_compile,
        "batch_size": batch_size,
        # Upstream's configs/inference.yaml defaults torch_deterministic to false and this
        # node does not change it, so the seed is recorded for provenance and reproduces in
        # practice on fixed hardware rather than being guaranteed to.
        "torch_deterministic": False,
        "hf_ckpt_revision": hf_ckpt_revision,
        "rfd4_proteina_version": package_version(),
        "gpu": gpu_info(),
        "elapsed_seconds": elapsed,
        "timeout_seconds": timeout_s,
        "cuda_oom_reported": saw_oom,
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
