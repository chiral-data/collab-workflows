#!/usr/bin/env python3
"""Refold each designed sequence with a predictor that never saw the design. GPU.

Contract
--------
Reads   mono_*.fasta, prediction_inputs.json
Writes  <design>_scores_rank_*.json   flat, prefixed per design -- see below
        <design>_*rank_*.pdb          all five ranks; node 04 selects rank_001
        prediction_manifest.json      per design: status, output file names, chain_order,
                                      chain_lengths, argv, elapsed, failure reason

Run with --dry-run to resolve and record every argv without starting the predictor. That
is 026's pattern (02_predict/script.py:37) and it is what makes the offline and
model-type guarantees testable on a machine with no GPU: the flags are the thing being
asserted, and they are decided here.

Flat output, prefixed per design
--------------------------------
Silva's ports match file NAMES, not paths: `silva validate` refuses an input like
`predictions/*` because it "can never match". So a per-design subdirectory would be
flattened on the way into node 04's inputs anyway, and the separation between designs has
to live in the filename instead.

That is not a workaround, it is the lesson. Flattening per-design output into one
directory is exactly how workflow-026 conflates designs -- its 02_predict "copies all
colabfold output flat into ./outputs/, so N designs' files land in one directory
undistinguished", and its 03_analyze_confidence then globs across all of them and takes
job_id from whichever sorts first (818954e9). The defence is that ColabFold names its
outputs after the record id and node 02 writes exactly one record per design, so every
file already carries the design it belongs to.

Run colabfold_batch into a scratch directory per design, then move its files out under
that prefix. On a name collision, keep a disambiguating prefix rather than overwriting:
workflow-037 node 03 flattens the same way and had to learn this twice -- first that the
real layout nests one level deeper than it assumed, then that flattening a tree can
collide and the old code would have overwritten one file with another silently
(cb3f6d23).

The manifest carries chain_order and chain_lengths forward from node 02 so node 04 reads
them from its own direct dependency rather than taking an extra edge back to node 02.
Passing the fact forward is the rule; re-deriving it is the defect (CLAUDE.md).

argv
----
    colabfold_batch <input> <outdir>
      --model-type    alphafold2_ptm          (mono; alphafold2_multimer_v3 for Phase 2)
      --msa-mode      single_sequence         (PARAM_MSA_MODE)
      --num-models    5                       (PARAM_NUM_MODELS)
      --num-recycle   3                       (PARAM_NUM_RECYCLE)

--model-type is not optional. The multimer model is what produces iptm at all, and
running a single chain through it is not the monomer prediction anyone means. Measured:
the real monomer scores file carries ptm and NO iptm at all, and the real multimer file
carries iptm 0.2 -- the key does not merely differ in value, it is absent (tests/
reference/).

--host-url is NEVER passed. Workflow-026 passes it unconditionally, defaulting to
https://api.colabfold.com (02_predict/script.py:58), which is how a workflow reaches a
third-party service without anyone deciding that it should. The weights are baked into
the image -- workflow-026's Dockerfile runs `python -m colabfold.download` at build time,
5.3 GB of them, measured -- so single_sequence mode has no reason to touch the network.

What this node can and cannot assert about the network: it can guarantee the flag is
absent, and it records that in the manifest. It cannot prove no egress happened from
inside its own process. The proof is the harness running it under `docker run
--network none`, which is how the offline claim was measured in the first place: a
35-residue monomer prediction completed in 61 s with no network at all. The fixture
asserting it therefore belongs to the harness, and the manifest records the flags it
would have needed.

Phase 1 is monomer-only and single_sequence is CORRECT here rather than a compromise: a
de novo binder has no homologs, and designed proteins are idealised, which is why the
self-consistency convention works without an alignment. It is wrong for a complex, for a
reason that is about the target rather than the binder -- see validation_pipeline.md §4.

Timeout: kill the process group, not the child
----------------------------------------------
subprocess.run's own timeout is not enough. Workflow-037 found that generation spawns
workers which go on holding the GPU after the direct child is signalled: a run outlived
the harness that started it, the client was killed, the container was not, and the device
had to be freed by hand (755d4af2). So: Popen in its own process group, under
PARAM_PREDICT_TIMEOUT_S, with kill_tree() escalating SIGTERM then SIGKILL. Verified
against a child that installs SIG_IGN for SIGTERM, which waits out the grace period and
dies of SIGKILL at rc=-9 (55edd8fc).

A stalled prediction holds the device exactly the same way a stalled generation does. The
budget here is unmeasured for this workflow (§11.5): recycles multiply it, and
workflow-026 records 6-12 recycles as the recommendation for multimers against the 3 that
suits a monomer.

Retain the output
-----------------
Stream it and keep it, so a failure is diagnosed rather than inferred. Both halves of
that have cost someone a day already: 037 could only infer an allocator OOM from
warnings the CUDA caching allocator prints while carrying on (755d4af2), and 029
discarded a prediction's ColabFold output unless the design passed, so a failing funnel
could not afterwards distinguish a bad design from a target that was never modelled --
an ambiguity that produced one wrong diagnosis and had to be resolved by computing the
missing number by hand, outside the pipeline (2359ffa4).

So the log of every design, passing or failing, is written as <design>_predict.log and
kept. A design that produced no scores file keeps its log too -- that is the case the log
exists for.

Assertions
----------
  * the prediction came back with prediction_inputs.json's expected_chains. A Phase 2
    complex that returns one chain means the ':' was dropped somewhere and the result is
    a monomer masquerading as a complex -- which is 026's exact failure, where
    independent monomers are labelled a multimer and --pair-mode unpaired_paired is set,
    and the flag changes nothing because there is no complex. Fail loudly.
  * a rank_001 scores file and a rank_001 PDB exist per design, found by the design's own
    filename prefix. Never glob across designs.
  * the residue count in the returned structure matches expected_residues.

Partial failure
---------------
Recorded per design; the batch continues. A design with no prediction is named in the
manifest and must reach the report as a failed row with a stated reason, not vanish from
it (§8.5).

Fails when
----------
  * no GPU (fail at scheduling via [resources], and again at check_gpu() to be legible)
  * every design failed to predict
  * the chain-count assertion fails for any design
"""
import glob
import json
import math
import os
import re
import shutil
import signal
import subprocess
import sys
import threading
import time
from collections import deque

# Phase 1. alphafold2_multimer_v3 is Phase 2's, and the choice is not a tuning knob: it
# decides whether iptm exists in the output at all.
MODEL_TYPE_MONO = "alphafold2_ptm"
MODEL_TYPE_CPLX = "alphafold2_multimer_v3"

# Never passed, by decision. Named here so the absence is searchable and so a future
# change has to delete a line that says why.
HOST_URL_FLAG = "--host-url"

TAIL_LINES = 40
# What the CUDA caching allocator prints while carrying on, which is the failure mode that
# looks like a slow job (755d4af2).
OOM_MARKERS = (
    "out of memory",
    "RESOURCE_EXHAUSTED",
    "CUDA_ERROR_OUT_OF_MEMORY",
    "Allocator (GPU_0_bfc) ran out of memory",
)
# Seconds between SIGTERM and SIGKILL. 037 waits 20 for generation; a prediction holding
# the device is the same situation.
KILL_GRACE_S = 20

SCORES_RANK_RE = re.compile(r"_scores_rank_(\d+)")
PDB_RANK_RE = re.compile(r"_rank_(\d+)")


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def numeric_param(env_name, display_name, default, cast):
    raw = os.environ.get(env_name)
    if raw is None or raw.strip() == "":
        return default
    try:
        value = cast(raw)
    except ValueError:
        kind = "a whole number" if cast is int else "a number"
        fail(f"{display_name} must be {kind}, got '{raw}'.")
    if not math.isfinite(value):
        fail(
            f"{display_name} must be a finite number, got '{raw}', "
            f"which parses as {value}."
        )
    return value


def positive_param(env_name, display_name, default, cast):
    value = numeric_param(env_name, display_name, default, cast)
    if value <= 0:
        fail(f"{display_name} must be positive, got {value}.")
    return value


def check_gpu():
    """Say out loud that there is no device, rather than failing 60 seconds deeper in.

    job.toml declares [resources] gpu = "required", so the scheduler is the real guard
    (97a16f7c). This exists because a scheduling failure and a missing driver inside the
    container read completely differently to whoever is looking at the log.
    """
    try:
        out = subprocess.run(
            ["nvidia-smi", "--query-gpu=name,memory.total", "--format=csv,noheader"],
            capture_output=True, text=True, timeout=60,
        )
    except (FileNotFoundError, subprocess.SubprocessError) as exc:
        fail(
            f"No usable GPU: {exc}. This node declares [resources] gpu = \"required\" and "
            "use_gpu = true, so a scheduler that honours either one should not have "
            "placed it here. Run with --dry-run to resolve and record every argv without "
            "a device."
        )
    if out.returncode != 0 or not out.stdout.strip():
        fail(
            "nvidia-smi reported no device "
            f"(rc={out.returncode}): {(out.stderr or out.stdout).strip()[:300]}"
        )
    for line in out.stdout.strip().splitlines():
        print(f"GPU: {line.strip()}", flush=True)


def kill_tree(proc):
    """End the whole prediction process group, not just the process we started."""
    for sig in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(os.getpgid(proc.pid), sig)
        except (ProcessLookupError, PermissionError, OSError):
            proc.kill() if sig == signal.SIGKILL else proc.terminate()
        try:
            proc.wait(timeout=KILL_GRACE_S)
            return
        except subprocess.TimeoutExpired:
            continue


def run_prediction(argv, timeout_s, log_path):
    """Run colabfold_batch under a wall-clock budget, streaming and retaining its output.

    Not subprocess.run(argv, timeout=...), for the two reasons 037 arrived at (755d4af2):
    its timeout signals the direct child only, while the predictor's workers go on holding
    the GPU; and the output has to be kept, because an allocator OOM is printed as a
    warning and the process then stops making progress without raising, exiting or
    changing its status.

    Returns (returncode, timed_out, saw_oom, tail).
    """
    proc = subprocess.Popen(
        argv,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        start_new_session=True,
    )
    tail = deque(maxlen=TAIL_LINES)
    # The reader is a daemon and the join below is bounded, so the snapshot of the tail can
    # in principle run while a line is still being appended: a grandchild that survives
    # the kill keeps the pipe open and the thread alive. Copying a deque that is being
    # mutated raises, and the traceback would land exactly where the timeout explanation
    # belongs. 037 made this impossible rather than merely unlikely; so does this.
    tail_lock = threading.Lock()
    state = {"oom": False}

    with open(log_path, "w") as log:
        def pump():
            # A grandchild that survives the kill keeps the pipe open and this thread
            # alive past the bounded join below, after which the log file is closed. A
            # write to a closed file would raise in a daemon thread and print a traceback
            # exactly where the timeout explanation belongs, so the thread swallows its
            # own end rather than reporting it as the failure.
            try:
                for line in proc.stdout:
                    sys.stdout.write(line)
                    sys.stdout.flush()
                    log.write(line)
                    log.flush()
                    with tail_lock:
                        tail.append(line.rstrip("\n"))
                    if any(marker in line for marker in OOM_MARKERS):
                        state["oom"] = True
            except (ValueError, OSError):
                return

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


def predictor_tail(captured, outputs, out_dir, design):
    """The predictor's last lines, from its stdout or from its own log file.

    Returns (lines, source, saw_oom_in_fallback).

    MEASURED, in control C6: colabfold_batch's stdout is block-buffered because it is a
    pipe rather than a tty, so when the process is KILLED the buffer dies with it. A
    completed prediction leaves a 3 KB `<design>_predict.log`; one killed at 35s leaves
    **zero bytes**, and the timeout reason then reads "Last lines: " with nothing after it.
    The one thing a reader needs at that moment is the one thing that was lost.

    But it is not actually lost. ColabFold logs through the logging module to its own
    `<design>_log.txt`, which is written as it goes and which collect() has already copied
    into ./outputs by the time this is called. So the tail falls back to that file and says
    which source it came from, because "the predictor said nothing" and "the predictor's
    stdout was buffered away" are different facts and a reason that conflates them is
    worse than one that admits the gap (§8.7).

    The OOM scan has the same blind spot for the same reason -- it only ever saw stdout --
    so the fallback file is scanned for the markers too. An OOM followed by a timeout was
    previously invisible: the allocator warns, the job keeps the GPU without progressing,
    the budget kills it, and every line that said so was in the lost buffer.
    """
    if captured:
        return captured, "the predictor's stdout", False

    # ColabFold's own log, named after the record rather than the design.
    candidates = [name for name in outputs if name.endswith("_log.txt")]
    for name in sorted(candidates):
        path = os.path.join(out_dir, name)
        try:
            with open(path, errors="replace") as handle:
                lines = [line.rstrip("\n") for line in handle if line.strip()]
        except OSError:
            continue
        if not lines:
            continue
        saw_oom = any(marker in line for line in lines for marker in OOM_MARKERS)
        return (
            lines[-TAIL_LINES:],
            f"colabfold's own {name} -- its stdout was block-buffered and lost when the "
            f"process was killed (C6)",
            saw_oom,
        )

    return [], (
        f"nothing: stdout was buffered away when the process was killed, and no "
        f"*_log.txt was collected for {design}"
    ), False


def oom_advice(saw_oom):
    if not saw_oom:
        return ""
    return (
        " The CUDA allocator reported an out-of-memory condition. It logs that as a "
        "warning and carries on, so the job keeps the GPU without making progress and "
        "never fails on its own. Reduce num_models, reduce num_recycle, or run on a GPU "
        "with more memory."
    )


def build_argv(record, in_dir, out_dir, num_models, num_recycle, msa_mode):
    """One invocation, one design. The flags are the contract; see the module docstring."""
    model_type = MODEL_TYPE_CPLX if record["mode"] == "cplx" else MODEL_TYPE_MONO
    argv = [
        "colabfold_batch",
        in_dir,
        out_dir,
        "--model-type", model_type,
        "--msa-mode", msa_mode,
        "--num-models", str(num_models),
        "--num-recycle", str(num_recycle),
    ]
    # Asserted here rather than trusted, because the whole offline claim rests on it and
    # the flag is one string away from being added back by someone copying from 026.
    assert HOST_URL_FLAG not in argv, "host-url must never be passed"
    return argv


def cross_check_cplx(fasta_path, a3m_path, record):
    """The fasta and the a3m must describe the same complex. Returns a reason or None.

    Node 02 writes two representations of one fact on the complex path: cplx_<d>.fasta
    (TARGET:BINDER, the declaration, where the ':' rule is enforceable) and
    pairedmsa_<d>.a3m (what the predictor consumes). Two representations are only better
    than one if they are checked against each other, which is what this does -- and it
    does it before any GPU time is spent, for 37c23bb5's reason.

    Each disagreement below produces a PREDICTION rather than an error if it goes
    unchecked, which is the §8.3 case: unserialize_msa slices the query row by the
    DECLARED lengths, so a wrong declaration returns two chains of the wrong sizes, and
    node 04 then averages a plddt slice that does not describe the binder.
    """
    try:
        with open(fasta_path) as handle:
            body = "".join(
                line.strip() for line in handle if line.strip()
                and not line.startswith(">")
            )
    except OSError as error:
        return f"could not read {os.path.basename(fasta_path)}: {error}"

    if body.count(":") != 1:
        return (
            f"{os.path.basename(fasta_path)} holds {body.count(':')} colon(s), not 1. "
            "Node 02 verified this on write, so the file has changed since."
        )
    target_seq, binder_seq = body.split(":")

    try:
        with open(a3m_path) as handle:
            lines = [line.rstrip("\n") for line in handle if line.strip()]
    except OSError as error:
        return f"could not read {os.path.basename(a3m_path)}: {error}"

    if not lines or not lines[0].startswith("#"):
        return f"{os.path.basename(a3m_path)} has no '#<len1>,<len2>' header"
    declared = lines[0][1:].split("\t")[0]
    expect = f"{len(target_seq)},{len(binder_seq)}"
    if declared != expect:
        return (
            f"{os.path.basename(a3m_path)} declares chain lengths {declared!r} but "
            f"{os.path.basename(fasta_path)}'s halves are {expect!r}. unserialize_msa "
            "slices the query by the declared lengths, so this returns two chains of the "
            "wrong sizes and no error."
        )

    recorded = ",".join(str(n) for n in record["chain_lengths"])
    if recorded != expect:
        return (
            f"node 02 recorded chain_lengths {recorded!r} but its own fasta's halves are "
            f"{expect!r}. Node 04 slices by the recorded value, so the two must agree."
        )

    rows = [lines[i + 1] for i in range(1, len(lines) - 1, 2) if lines[i].startswith(">")]
    if not rows:
        return f"{os.path.basename(a3m_path)} holds no sequence rows"
    query = rows[0]
    if query != target_seq + binder_seq:
        return (
            f"{os.path.basename(a3m_path)}'s query row is not the fasta's "
            "TARGET+BINDER concatenation. The two files describe different complexes, and "
            "the a3m is the one the predictor reads."
        )
    if record["chain_order"] != ["target", "binder"]:
        return (
            f"node 02 recorded chain_order {record['chain_order']!r} for a complex "
            "record. The a3m is written target-first, so any other order means node 04 "
            "would slice the wrong half (029 slices plddt[:binder_length])."
        )
    return None


def check_cplx_msa_mode(msa_mode, a3m_name):
    """Refuse the combination that silently discards a precomputed alignment.

    MEASURED, in tests/probe_msa_clamp.sh. colabfold/batch.py:137-142:

        if msa_mode == "single_sequence":
            num_seqs = 1
            ...
            max_seq = min(num_seqs, max_seq)

    and batch.py:258-270 shows that when a3m_lines is not None the alignment is read by
    unserialize_msa() and msa_mode is never consulted for MSA GENERATION at all. Put
    together: an a3m handed in with --msa-mode single_sequence is parsed, is re-saved as
    <job>.a3m, and gets a coverage.png -- and the model is fed ONE sequence.

    Every artifact that looks like proof an alignment was used is written from
    unserialize_msa()'s output BEFORE the clamp (batch.py:284-287), which is why §11.1's
    spike concluded "coverage.png written, so an MSA was really consumed" and could not
    have known otherwise. Measured on the same a3m, same seed, offline, both arms:

        msa_mode             iptm    ptm    max_pae   mean plddt
        single_sequence      0.270   0.510  23.27     75.3
        mmseqs2_uniref_env   0.120   0.440  27.84     66.2

    The numbers differ, so the clamp is real. (Which arm scores HIGHER says nothing about
    alignments: that reference file's homologs are synthetic variants built to probe the
    format, and a tiny synthetic heavily-gapped MSA can be worse than none. The finding is
    the difference, not its sign.)

    So a Phase 2 run at single_sequence would report interface numbers measured at
    effectively single_sequence, while config.json, the re-saved a3m and coverage.png all
    said an alignment was used. That is §4's exact failure mode and it is invisible --
    which is why it is refused by name here rather than warned about.

    mmseqs2_uniref_env is safe to pass ALONGSIDE an a3m: use_msa and use_env only feed
    get_msa_and_templates, which is not called when a3m_lines is not None. The caller
    guarantees the a3m is what is being passed -- see the fail() below it.
    """
    if msa_mode == "single_sequence":
        fail(
            f"the complex path was asked for with msa_mode='single_sequence' and a "
            f"precomputed alignment ({a3m_name}). That combination silently throws the "
            "alignment away: colabfold sets max_seq = 1 when msa_mode is "
            "single_sequence (batch.py:137-142), so the model is fed one sequence, while "
            "the re-saved a3m, config.json and coverage.png all record that an alignment "
            "was used. Measured in tests/probe_msa_clamp.sh: the same a3m gives iptm "
            "0.270 at single_sequence and 0.120 at mmseqs2_uniref_env, so the clamp "
            "changes what the model sees.\n"
            "  §4's whole decision is that interface metrics measured without a target "
            "alignment are about the prediction mode rather than the design (029: iPTM "
            "0.100 against 0.860). Reporting them under a Phase 2 heading is the failure "
            "this refusal exists to prevent.\n"
            "  Set PARAM_MSA_MODE_CPLX to mmseqs2_uniref_env (the default), which does "
            "NOT reach the network when an a3m is supplied, or run Phase 1."
        )


def read_pdb_chains(path):
    """Return {chain_id: residue_count} from a predicted PDB, by CA atoms.

    Deliberately a different parser from node 01's mmCIF reader: this file comes from the
    predictor and the one node 01 read came from the generator, and sharing a parser
    between the thing being measured and the thing it is measured against is how a check
    passes for the wrong reason (§8.4).
    """
    counts = {}
    with open(path) as handle:
        for line in handle:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            if line[12:16].strip() != "CA":
                continue
            chain = line[21].strip() or "_"
            counts[chain] = counts.get(chain, 0) + 1
    return counts


def collect(scratch, out_dir, design):
    """Move the predictor's files out flat, prefixed by the design, without overwriting.

    Returns the list of basenames written. ColabFold already names its outputs after the
    record id, which is the design name, so most files need no prefix at all -- but a file
    that does not carry it gets one, and a collision keeps a disambiguating parent-
    directory prefix rather than silently replacing a file (cb3f6d23).
    """
    written = []
    for src in sorted(glob.glob(os.path.join(scratch, "**", "*"), recursive=True)):
        if not os.path.isfile(src):
            continue
        base = os.path.basename(src)
        name = base if base.startswith(f"{design}_") else f"{design}_{base}"
        dest = os.path.join(out_dir, name)
        if os.path.exists(dest):
            parent = os.path.basename(os.path.dirname(src)) or "dir"
            name = f"{design}_{parent}_{base}"
            dest = os.path.join(out_dir, name)
            print(
                f"NOTE: '{base}' would have overwritten an earlier file; kept as "
                f"'{name}'. Flattening a tree can collide, and overwriting silently is "
                "what cb3f6d23 fixed.",
                flush=True,
            )
            if os.path.exists(dest):
                print(
                    f"WARNING: '{name}' also exists. Skipping '{src}' rather than "
                    "overwriting it; the predictor's own output is kept in the scratch "
                    "directory.",
                    flush=True,
                )
                continue
        shutil.move(src, dest)
        written.append(name)
    return written


def rank_one(names, pattern, suffix, prefer=None):
    """The rank_001 file among these, by the design's OWN names. Never a global glob.

    `prefer` breaks a tie towards a substring. It exists for one case: with --amber,
    ColabFold writes both relaxed and unrelaxed structures at the same rank, and
    'relaxed' sorts before 'unrelaxed', so a plain lexicographic tie-break would quietly
    hand node 04 the force-field-minimised coordinates. The self-consistency RMSD is
    about what the predictor produced from the sequence, not about what a minimiser did
    to it afterwards. This node passes no --amber, so the case should not arise -- which
    is exactly why it would not be noticed if it did.
    """
    candidates = []
    for name in names:
        if not name.endswith(suffix):
            continue
        match = pattern.search(name)
        if match and int(match.group(1)) == 1:
            candidates.append(name)
    if not candidates:
        return None
    if prefer:
        preferred = [name for name in candidates if prefer in name]
        if preferred:
            return sorted(preferred)[0]
    return sorted(candidates)[0]


def main():
    dry_run = "--dry-run" in sys.argv

    num_models = positive_param("PARAM_NUM_MODELS", "num_models", 5, int)
    num_recycle = positive_param("PARAM_NUM_RECYCLE", "num_recycle", 3, int)
    timeout_s = positive_param(
        "PARAM_PREDICT_TIMEOUT_S", "predict_timeout_s", 1800, int
    )
    # The MONO path's msa_mode. Still locked to single_sequence, and still correct rather
    # than a compromise: a de novo binder predicted alone has no homologs to align (§4),
    # and anything else would need the network, which this node never has.
    msa_mode = (os.environ.get("PARAM_MSA_MODE") or "single_sequence").strip()
    if msa_mode != "single_sequence":
        fail(
            f"msa_mode '{msa_mode}' is not available for the self-consistency prediction, "
            "which predicts the binder alone -- a de novo sequence has no homologs to "
            "align (§4). Anything else would also need the network, and this node never "
            "passes --host-url. The complex path has its own setting, "
            "PARAM_MSA_MODE_CPLX, because it consumes a precomputed alignment instead."
        )

    # The CPLX path's msa_mode, separate because the two paths are different predictions
    # and sharing the value is how one silently gets the other's. See
    # check_cplx_msa_mode() for the measurement that makes single_sequence unusable here
    # even though it sounds like the safe choice.
    msa_mode_cplx = (
        os.environ.get("PARAM_MSA_MODE_CPLX") or "mmseqs2_uniref_env"
    ).strip()
    # Not inherited. 3 recycles is right for a monomer; 026 records 6-12 for multimers,
    # and inheriting is how a multimer silently gets a monomer's budget.
    num_recycle_cplx = positive_param(
        "PARAM_NUM_RECYCLE_CPLX", "num_recycle_cplx", 6, int
    )

    os.makedirs("./outputs", exist_ok=True)

    inputs_path = "./inputs/prediction_inputs.json"
    if not os.path.exists(inputs_path):
        fail(
            "No 'prediction_inputs.json' in ./inputs/. This node runs one invocation per "
            "record that node 02 declared, and it does not glob ./inputs/ for FASTA files "
            "instead: globbing is how a stale file from an earlier run gets predicted, and "
            "the manifest is the list of what this run is actually for."
        )
    with open(inputs_path) as handle:
        declared = json.load(handle)
    records = declared.get("inputs") or []
    if not records:
        fail(
            "node 02 declared no prediction inputs. Every design was skipped, or node 02 "
            "did not run. Either way there is nothing to predict, and reporting success "
            "here would make an empty report look like a clean one."
        )

    if not dry_run:
        check_gpu()
    else:
        print(
            "DRY RUN: resolving every argv and writing the manifest, without starting the "
            "predictor and without touching a device.",
            flush=True,
        )

    results = []
    for record in records:
        design = record["design"]
        fasta = os.path.join("./inputs", record["path"])
        if not os.path.exists(fasta):
            results.append({
                **{k: record[k] for k in
                   ("design", "mode", "chain_order", "chain_lengths",
                    "expected_chains", "expected_residues")},
                "status": "failed",
                "reason": (
                    f"node 02 declared '{record['path']}' but it is not in ./inputs/. The "
                    "declaration and the file disagree, so the design is recorded as "
                    "failed rather than skipped -- it was meant to be predicted."
                ),
                "argv": None,
                "elapsed_s": None,
                "outputs": [],
            })
            print(f"  {design}: FAILED, declared input missing", flush=True)
            continue

        # One input directory per design, holding exactly one file. colabfold_batch takes
        # a directory, and a directory with two designs in it is 026's conflation bug.
        scratch_in = os.path.join("./cf", design, "in")
        scratch_out = os.path.join("./cf", design, "out")
        os.makedirs(scratch_in, exist_ok=True)
        os.makedirs(scratch_out, exist_ok=True)

        this_msa_mode = msa_mode
        this_recycle = num_recycle
        msa_from = "msa_mode"

        if record["mode"] == "cplx":
            # The a3m is the input, NOT the fasta, and the fasta is deliberately not
            # copied in beside it: colabfold_batch takes a directory, so two readable
            # inputs would be two predictions.
            #
            # And there is no fallback to the fasta if the a3m is missing. A bare
            # TARGET:BINDER fasta under mmseqs2_uniref_env would make colabfold align BOTH
            # chains -- sending the designed binder to the MSA server, which is the one
            # thing §4 guarantees never happens. So an absent a3m is a failure, and that
            # is also what makes the msa_mode below safe: mmseqs2_uniref_env is only ever
            # passed when the input really is a precomputed alignment.
            a3m_name = record.get("paired_msa")
            if not a3m_name:
                results.append({
                    **{k: record[k] for k in
                       ("design", "mode", "chain_order", "chain_lengths",
                        "expected_chains", "expected_residues")},
                    "status": "failed",
                    "reason": (
                        "node 02 declared a complex record with no paired_msa. The "
                        "complex path needs the precomputed alignment; predicting the "
                        "TARGET:BINDER fasta instead would align the designed binder as "
                        "well, which §4 forbids, so this is refused rather than "
                        "downgraded."
                    ),
                    "argv": None, "elapsed_s": None, "outputs": [],
                })
                print(f"  {design}: FAILED, complex record with no paired a3m", flush=True)
                continue

            a3m_path = os.path.join("./inputs", a3m_name)
            if not os.path.exists(a3m_path):
                results.append({
                    **{k: record[k] for k in
                       ("design", "mode", "chain_order", "chain_lengths",
                        "expected_chains", "expected_residues")},
                    "status": "failed",
                    "reason": (
                        f"node 02 declared '{a3m_name}' but it is not in ./inputs/. Not "
                        "falling back to the fasta: that would align the designed binder."
                    ),
                    "argv": None, "elapsed_s": None, "outputs": [],
                })
                print(f"  {design}: FAILED, paired a3m missing", flush=True)
                continue

            # Cross-check the two representations against each other before spending GPU
            # time. The a3m's declared lengths are what unserialize_msa slices by, and the
            # fasta's halves are what node 02 said it was asking for; if they disagree,
            # the prediction comes back as two chains of the wrong sizes with no error,
            # and node 04 then slices a plddt array by numbers that do not describe it.
            problem = cross_check_cplx(fasta, a3m_path, record)
            if problem:
                results.append({
                    **{k: record[k] for k in
                       ("design", "mode", "chain_order", "chain_lengths",
                        "expected_chains", "expected_residues")},
                    "status": "failed", "reason": problem,
                    "argv": None, "elapsed_s": None, "outputs": [],
                })
                print(f"  {design}: FAILED, {problem[:80]}", flush=True)
                continue

            this_msa_mode = msa_mode_cplx
            this_recycle = num_recycle_cplx
            msa_from = f"the precomputed alignment {a3m_name}"
            check_cplx_msa_mode(this_msa_mode, a3m_name)
            shutil.copyfile(a3m_path, os.path.join(scratch_in, a3m_name))
        else:
            shutil.copyfile(fasta, os.path.join(scratch_in, record["path"]))

        argv = build_argv(
            record, scratch_in, scratch_out, num_models, this_recycle, this_msa_mode
        )
        print(f"  {design}: {' '.join(argv)}", flush=True)
        if record["mode"] == "cplx":
            print(
                f"    the alignment comes from {msa_from}, not from --msa-mode. "
                f"--msa-mode is passed as {this_msa_mode!r} because single_sequence "
                "would clamp max_seq to 1 and discard it (batch.py:137-142, measured in "
                "tests/probe_msa_clamp.sh); it reaches no network when an a3m is the "
                "input. The binder carries no alignment of its own.",
                flush=True,
            )

        entry = {
            **{k: record[k] for k in
               ("design", "mode", "chain_order", "chain_lengths",
                "expected_chains", "expected_residues")},
            "record_id": record["record_id"],
            "argv": argv,
            "model_type": argv[argv.index("--model-type") + 1],
            "msa_mode": this_msa_mode,
            "msa_source": msa_from,
            "num_models": num_models,
            "num_recycle": this_recycle,
            "host_url_passed": False,
            "outputs": [],
            "log": None,
            "elapsed_s": None,
        }
        if record["mode"] == "cplx":
            entry["paired_msa"] = record.get("paired_msa")
            entry["n_target_homologs_in_msa"] = record.get("n_target_homologs_in_msa")
            entry["binder_aligned"] = False

        if dry_run:
            entry["status"] = "dry-run"
            entry["reason"] = None
            results.append(entry)
            continue

        log_name = f"{design}_predict.log"
        started = time.monotonic()
        returncode, timed_out, saw_oom, tail = run_prediction(
            argv, timeout_s, os.path.join("./outputs", log_name)
        )
        entry["elapsed_s"] = round(time.monotonic() - started, 1)
        entry["log"] = log_name
        entry["saw_oom"] = saw_oom

        # Collected whether or not the run succeeded: 029 discarded ColabFold output
        # unless the design passed, and a failing funnel then could not be diagnosed
        # afterwards (2359ffa4).
        entry["outputs"] = collect(scratch_out, "./outputs", design)

        if timed_out or returncode != 0:
            # Collected first, above, so colabfold's own log is already on disk and can
            # stand in for the stdout that a kill throws away. See predictor_tail().
            tail, tail_source, oom_in_log = predictor_tail(
                tail, entry["outputs"], "./outputs", design
            )
            saw_oom = saw_oom or oom_in_log
            entry["saw_oom"] = saw_oom
            entry["tail_source"] = tail_source

        if timed_out:
            entry["status"] = "failed"
            entry["reason"] = (
                f"timed out after {timeout_s}s and the whole process group was killed."
                + oom_advice(saw_oom)
                + f" Last lines, from {tail_source}: " + " | ".join(tail[-5:])
            )
        elif returncode != 0:
            entry["status"] = "failed"
            entry["reason"] = (
                f"colabfold_batch exited {returncode}." + oom_advice(saw_oom)
                + f" Last lines, from {tail_source}: " + " | ".join(tail[-5:])
            )
        else:
            scores = rank_one(entry["outputs"], SCORES_RANK_RE, ".json")
            pdb = rank_one(entry["outputs"], PDB_RANK_RE, ".pdb", prefer="unrelaxed")
            if scores is None or pdb is None:
                entry["status"] = "failed"
                entry["reason"] = (
                    "colabfold_batch exited 0 but produced "
                    f"{'no rank_001 scores file' if scores is None else ''}"
                    f"{' and ' if scores is None and pdb is None else ''}"
                    f"{'no rank_001 PDB' if pdb is None else ''}"
                    f" among {len(entry['outputs'])} file(s). Exit codes see nothing "
                    "(§8.3): the check is on the artifact, not the status."
                )
            else:
                entry["scores_rank_001"] = scores
                entry["pdb_rank_001"] = pdb
                chains = read_pdb_chains(os.path.join("./outputs", pdb))
                entry["returned_chains"] = chains
                total = sum(chains.values())
                if len(chains) != record["expected_chains"]:
                    # Loud and fatal, not per-design. A complex that came back as one
                    # chain means the ':' was dropped, and every interface number
                    # computed from it would be a monomer's confidence reported as an
                    # interface -- 026's exact failure.
                    fail(
                        f"'{design}' was predicted as {len(chains)} chain(s) "
                        f"{sorted(chains)} but node 02 declared "
                        f"{record['expected_chains']}. For a cplx_ input this means the "
                        "':' was dropped somewhere and the result is a monomer "
                        "masquerading as a complex: iptm and pae_interaction computed "
                        "from it would describe one chain's confidence and be reported as "
                        "an interface. Refusing the whole run rather than this design, "
                        "because the cause is in the inputs and the next design has it "
                        "too."
                    )
                if total != record["expected_residues"]:
                    entry["status"] = "failed"
                    entry["reason"] = (
                        f"the prediction came back with {total} residue(s) but node 02 "
                        f"declared {record['expected_residues']}. Node 04 maps the two "
                        "structures by POSITION, so a length mismatch cannot be scored -- "
                        "truncating to the shorter one would report an RMSD over a "
                        "different set of residues than the one named (§8.1)."
                    )
                else:
                    entry["status"] = "ok"
                    entry["reason"] = None

        results.append(entry)
        state = entry["status"].upper()
        print(
            f"  {design}: {state} in {entry['elapsed_s']}s, "
            f"{len(entry['outputs'])} file(s)"
            + (f" -- {entry['reason'][:160]}" if entry.get("reason") else ""),
            flush=True,
        )

    ok = [e for e in results if e["status"] in ("ok", "dry-run")]
    failed = [e for e in results if e["status"] == "failed"]

    manifest = {
        "n_declared": len(records),
        "n_ok": len([e for e in results if e["status"] == "ok"]),
        "n_failed": len(failed),
        "dry_run": dry_run,
        "phase": declared.get("phase", 1),
        # Per path, because there are now two kinds of prediction in one run and a single
        # pair of workflow-level values would describe only one of them. Node 05 reads
        # these to say which question each row answers.
        "model_type": MODEL_TYPE_MONO,
        "msa_mode": msa_mode,
        "num_recycle": num_recycle,
        "cplx": None if not any(e["mode"] == "cplx" for e in results) else {
            "model_type": MODEL_TYPE_CPLX,
            "msa_mode": msa_mode_cplx,
            "num_recycle": num_recycle_cplx,
            "msa_from": "the precomputed alignment written by 02 from 01b's target a3m",
            "binder_aligned": False,
            # The reason the msa_mode above is not single_sequence, recorded where a
            # reader of the manifest will see it rather than only in the source.
            "why_not_single_sequence": (
                "colabfold clamps max_seq to 1 when msa_mode is single_sequence "
                "(batch.py:137-142), which discards a precomputed alignment while still "
                "re-saving the a3m and writing coverage.png. Measured in "
                "tests/probe_msa_clamp.sh: same a3m, iptm 0.270 at single_sequence "
                "against 0.120 at mmseqs2_uniref_env. mmseqs2_uniref_env reaches no "
                "network when an a3m is the input (batch.py:258-270)."
            ),
        },
        "num_models": num_models,
        "predict_timeout_s": timeout_s,
        # The offline guarantee, as far as this node can state it: the flag is never
        # passed, and the proof of no egress is the harness's --network none.
        "host_url_passed": False,
        "network_expected": "none",
        "weights": "baked into the image at build time (026 Dockerfile, measured 5.3 GB)",
        # Carried forward from node 02 so node 04 reads the chain order from its own
        # direct dependency instead of taking an extra edge back (CLAUDE.md).
        "chain_order_source": "02-build-prediction-inputs",
        "designed_chain_sources": declared.get("designed_chain_sources"),
        "securedna_screened": False,
        "predictions": results,
    }
    with open("./outputs/prediction_manifest.json", "w") as handle:
        json.dump(manifest, handle, indent=2)

    if failed:
        print(
            f"WARNING: {len(failed)} of {len(records)} design(s) have no usable "
            f"prediction: {', '.join(e['design'] for e in failed)}. Each is recorded in "
            "prediction_manifest.json with its reason and its retained log, and node 05 "
            "renders it as a failed row rather than omitting it (§8.5). The batch "
            "continued.",
            flush=True,
        )
    if not ok:
        fail(
            f"All {len(records)} design(s) failed to predict. The per-design reasons and "
            "logs are in ./outputs/. Reporting success with nothing scored would make an "
            "empty report look like a clean one."
        )

    print(
        f"Wrote prediction_manifest.json: "
        f"{manifest['n_ok']} ok, {manifest['n_failed']} failed"
        + (" (dry run)" if dry_run else ""),
        flush=True,
    )


if __name__ == "__main__":
    main()
