#!/usr/bin/env python3
"""Stage and run every CI fixture in fixtures.json, then hand the tree to check_fixtures.py.

This file does the WORK; it decides nothing. What each case is, what it needs staged and
what must come out of it all live in fixtures.json, and the judgement lives in
check_fixtures.py. Three files rather than one because of how 037's suite went wrong: a
runner that also judged could only report pass or fail, and when a case failed there was
nothing to look at afterwards. So every step here leaves its inputs, its outputs, its log
and its exit code on disk, and the checker reads that tree rather than watching the run.

What it replicates
------------------
Silva copies each dependency's declared outputs into the dependent's ./inputs/ and runs
the node with ./inputs/ and ./outputs/ relative to its own directory. That is all a node
sees, so that is all this gives it. Output globs come from each node's own job.toml rather
than being restated here: a harness that carried its own copy of the port list would keep
passing after job.toml changed, which is the shape of defect this suite exists to catch.

The complex path
----------------
Node 01b runs on every case, as it does in the DAG; on the Phase 1 path it writes a
target_alignment.json saying it was a no-op, which is what makes "did 01b run?" a question
with an answer (§8.3). A case takes the complex path when its own `params` set
PARAM_SCORE_INTERFACE, and there is no harness token for it: a second place to say so is a
second place for it to disagree with the node. When it is set, the harness builds a
stand-in target alignment out of the inventory node 01 has just written
(make_target_msa.py) and 01b reads it from a file, because the suite does not touch the
network. That tool refuses nothing -- a monomer, or a batch disagreeing about its target,
is a refusal belonging to 01b, and a harness that pre-empted it would leave the node's
guard never shown refusing (§8.13).

Container, or not
-----------------
The gate is the suite running in `colabfold:2026_10_08`, because that is where the
workflow runs and biopython 1.85 in particular behaves differently from other versions on
a FASTA comment line (see node 05). If the image is present this file uses it. If it is
not, it falls back to the host interpreter and says so in the banner, in run.json and in
the checker's header, every time -- a fallback that does not announce itself is how a
suite ends up certifying something it never ran.

Stand-ins on PATH
-----------------
Two predict cases need a predictor that misbehaves on purpose: one that ignores SIGTERM
and leaves a worker behind, and one that fails for a single design. Node 03 resolves
`colabfold_batch` from PATH and that is deliberately not overridable by an environment
variable, so the harness writes a stand-in named `colabfold_batch` into the case's own
bin/ and prepends it. `nvidia-smi` is stood in for the same way, since those two cases are
about the subprocess contract rather than about device detection. Both stand-ins print
what they are, into the log the checker reads.

Run
---
  python3 tests/run_fixtures.py [--out DIR] [--only ID,...] [--section S] [--no-docker]

Default DIR is ./fixture-runs/<timestamp>. The last run is also linked as
./fixture-runs/latest.
"""
import datetime
import glob
import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WF = os.path.dirname(HERE)                       # workflows/workflow-040
WORKFLOWS = os.path.dirname(WF)                  # workflows/
IMAGE = "colabfold:2026_10_08"
GENERATOR = os.path.join(WORKFLOWS, "workflow-037", "tests", "make_designs.py")
METRICS_NODE = os.path.join(WORKFLOWS, "workflow-037", "04-analyze", "script.py")
PREDICTIONS = os.path.join(HERE, "make_predictions.py")
TARGET_MSA = os.path.join(HERE, "make_target_msa.py")

NODES = {
    "n1": "01-ingest-designs",
    # Not "n1b": that step name is already taken by also_run()'s second tree, and a step
    # name meaning two things is how a harness comes to stage the wrong directory.
    "nalign": "01b-align-target",
    "n2": "02-build-prediction-inputs",
    "n3": "03-predict",
    "n4": "04-score",
    "n5": "05-filter-and-report",
}

# Which step feeds which, mirroring .chiral/workflow.toml's [dependencies]. Read from the
# toml would be better still, but a toml parser is not guaranteed in every interpreter
# this may run under, so it is restated here and asserted against the toml at startup.
EDGES = {
    "n1": [],
    "nalign": ["n1"],
    "n2": ["n1", "nalign"],
    "n3": ["n2"],
    "n4": ["n1", "n3"],
    "n5": ["n4"],
}

# A predictor that ignores SIGTERM and leaves a worker behind. This is 037's measured
# failure (755d4af2): subprocess's own timeout signals the direct child while the workers
# carry on, so the client died and the GPU stayed busy. Both processes install `trap ''
# TERM`, so only the SIGKILL escalation can end them. The worker writes its pid where the
# checker can find it and look for it afterwards.
#
# NOT `setsid`. The first version of this stand-in put the worker in a new session, and it
# survived -- correctly, because a process-group kill cannot reach a process that has left
# the group. That is a real limit of the fix rather than a defect in it, and chasing it
# would mean walking /proc or running the predictor in a cgroup. What 037 measured, and
# what this stand-in has to reproduce, is a worker that inherits the group and ignores the
# polite signal. The limit is written down in 03-predict/script.py and in
# validation_pipeline.md §8.14 instead of being hidden behind a stand-in nobody could pass.
STANDIN_SIGTERM = """#!/bin/bash
echo "STAND-IN colabfold_batch: ignores SIGTERM and leaves a worker in the same group."
echo "argv: $*"
trap '' TERM
bash -c 'echo $$ > "$STANDIN_PIDFILE"; trap "" TERM; sleep 600' &
echo "worker started"
sleep 600
echo "STAND-IN: slept the whole way, which should never be reached"
"""

# A predictor that says nothing on stdout and writes its progress only to its own log
# file, then hangs until it is killed.
#
# This is the REAL colabfold_batch's behaviour under a kill, measured in control C6: its
# stdout is a pipe rather than a tty, so it is block-buffered, and a process killed before
# the buffer fills takes every line with it. A completed prediction leaves a 3 KB
# `<design>_predict.log`; one killed at 35s leaves zero bytes, and the timeout reason read
# "Last lines: " with nothing after it.
#
# STANDIN_SIGTERM prints to stdout and flushes line by line, so it could never have shown
# this -- the Tier 1 check passed on output the real predictor does not produce at that
# moment. That is §8.4 with the stand-in as the culprit: the check was right about the
# stand-in and wrong about the world. This one reproduces the world.
STANDIN_SILENT_LOGFILE = """#!/bin/bash
# No stdout at all, deliberately: everything goes to the log file colabfold writes.
out_dir="$2"
mkdir -p "$out_dir"
name="$(head -1 "$1"/*.fasta | sed 's/^>//' | tr -d '\\r' | awk '{print $1}')"
{
  echo "2026-11-21 10:46:04,178 Running colabfold 1.6.1"
  echo "2026-11-21 10:46:05,897 Running on GPU"
  echo "2026-11-21 10:46:05,942 Query 1/1: $name (length 60)"
  echo "STAND-IN logfile: this line exists only in the predictor's own log"
} > "$out_dir/log.txt"
trap '' TERM
sleep 600
"""

# A predictor that fails for exactly one design and succeeds for the rest. It writes a
# plausible scores file and PDB for the ones it does not fail, so the batch has real
# output to carry forward -- 029 threw that output away and then could not tell a bad
# design from a target that was never modelled (2359ffa4).
STANDIN_PARTIAL = """#!/bin/bash
echo "STAND-IN colabfold_batch: fails for $STANDIN_FAIL_FOR, succeeds otherwise."
echo "argv: $*"
in_dir="$1"; out_dir="$2"
# The RECORD ID, read from the header -- not the filename. ColabFold names its outputs
# after the record id, and node 02's file is mono_<design>.fasta while the record inside
# it is <design>, so a stand-in that used the filename would be testing a name the real
# predictor never produces.
name="$(head -1 "$in_dir"/*.fasta | sed 's/^>//' | tr -d '\\r' | awk '{print $1}')"
if [ "$name" = "$STANDIN_FAIL_FOR" ]; then
  echo "STAND-IN: refusing to predict $name" >&2
  exit 3
fi
mkdir -p "$out_dir"
suffix="rank_001_alphafold2_ptm_model_1_seed_000"
n=$(grep -v '^>' "$in_dir"/*.fasta | tr -d '\\n' | wc -c)
python3 - "$out_dir/${name}_scores_${suffix}.json" "$n" <<'PY'
import json, sys
path, n = sys.argv[1], int(sys.argv[2])
json.dump({"plddt": [88.0] * n,
           "pae": [[1.0] * n for _ in range(n)],
           "max_pae": 1.0, "ptm": 0.8}, open(path, "w"))
PY
i=1
: > "$out_dir/${name}_unrelaxed_${suffix}.pdb"
while [ "$i" -le "$n" ]; do
  printf 'ATOM  %5d  CA  ALA A%4d    %8.3f%8.3f%8.3f  1.00 50.00           C\\n' \\
    "$i" "$i" "$(echo "$i" | awk '{print $1*3.8}')" 0 0 \\
    >> "$out_dir/${name}_unrelaxed_${suffix}.pdb"
  i=$((i+1))
done
echo "END" >> "$out_dir/${name}_unrelaxed_${suffix}.pdb"
echo "STAND-IN: wrote $n-residue output for $name"
"""

STANDIN_NVIDIA_SMI = """#!/bin/bash
# Stand-in. These cases are about the subprocess contract -- the process-group kill and
# one design failing without taking the batch down -- not about device detection, and
# node 03 resolves nvidia-smi from PATH exactly as it resolves the predictor.
echo "STAND-IN nvidia-smi (no device present), 0 MiB"
"""


def say(message):
    print(message, flush=True)


def tokens(harness):
    """Parse a fixture's `harness` field. mp= takes the rest of the string."""
    out = {}
    text = (harness or "").strip()
    cut = re.search(r"(^|\s)mp=", text)
    if cut:
        out["mp"] = text[cut.end():].strip()
        text = text[:cut.start()]
    for part in text.split():
        if "=" not in part:
            raise SystemExit(f"run_fixtures: harness token '{part}' has no '='")
        key, value = part.split("=", 1)
        out[key] = value
    return out


def job_outputs(node):
    """The output globs this node declares, read from its own job.toml."""
    path = os.path.join(WF, node, ".chiral", "job.toml")
    text = open(path).read()
    block = re.search(r"^outputs\s*=\s*\[(.*?)\]", text, re.S | re.M)
    if not block:
        raise SystemExit(f"run_fixtures: no outputs list in {path}")
    return re.findall(r'"([^"]+)"', block.group(1))


def check_edges():
    """Fail loudly if EDGES has drifted from .chiral/workflow.toml."""
    text = open(os.path.join(WF, ".chiral", "workflow.toml")).read()
    block = re.search(r"\[dependencies\](.*?)(\n\[|\Z)", text, re.S)
    declared = {}
    for line in (block.group(1) if block else "").splitlines():
        match = re.match(r'\s*([\w-]+)\s*=\s*\[(.*)\]', line)
        if match:
            declared[match.group(1)] = re.findall(r'"([^"]+)"', match.group(2))
    by_node = {NODES[k]: [NODES[p] for p in v] for k, v in EDGES.items()}
    for node, parents in declared.items():
        if sorted(by_node.get(node, [])) != sorted(parents):
            raise SystemExit(
                f"run_fixtures: EDGES says {node} depends on {by_node.get(node)} but "
                f".chiral/workflow.toml says {parents}. The harness is staging the wrong "
                "inputs; fix EDGES rather than this check."
            )
    return declared


class Runner:
    """Runs a python script with ./inputs and ./outputs, in the image when there is one.

    Inside the container the repo and the run tree are mounted at THEIR OWN host paths,
    so an absolute path means the same thing on both sides. Mounting the case at /work
    instead would mean every path this harness hands to a tool as an argument was valid
    in one mode and not the other, and the suite would then only really work in whichever
    mode it was last debugged in.
    """

    def __init__(self, use_docker, root):
        self.use_docker = use_docker
        self.root = os.path.abspath(root)

    def run(self, script, cwd, env, argv=(), path_prepend=None, timeout=900):
        if self.use_docker:
            command = [
                "docker", "run", "--rm", "--network", "none",
                "-v", f"{WORKFLOWS}:{WORKFLOWS}:ro",
                "-v", f"{self.root}:{self.root}",
                "-w", os.path.abspath(cwd),
            ]
            for key, value in sorted(env.items()):
                command += ["-e", f"{key}={value}"]
            if path_prepend:
                standin = os.path.abspath(path_prepend)
                command += ["-e",
                            f"PATH={standin}:/usr/local/bin:/usr/bin:/bin"]
            command += [IMAGE, "python", script] + list(argv)
            shell_env = dict(os.environ)
        else:
            command = [sys.executable, script] + list(argv)
            shell_env = dict(os.environ)
            shell_env.update({k: str(v) for k, v in env.items()})
            if path_prepend:
                shell_env["PATH"] = (
                    os.path.abspath(path_prepend) + os.pathsep + shell_env["PATH"]
                )
        try:
            done = subprocess.run(
                command, cwd=cwd, capture_output=True, text=True,
                timeout=timeout, env=shell_env,
            )
            return done.returncode, done.stdout, done.stderr, " ".join(command)
        except subprocess.TimeoutExpired as exc:
            return (
                124,
                (exc.stdout or b"").decode(errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or ""),
                f"run_fixtures: the harness itself timed out after {timeout}s. The node "
                "under test has its own budget and should have ended first.",
                " ".join(command),
            )


class Case:
    def __init__(self, runner, root, section, entry):
        self.runner = runner
        self.section = section
        self.entry = entry
        self.id = entry["id"]
        self.dir = os.path.join(root, section, self.id)
        self.tok = tokens(entry.get("harness"))
        self.params = {k: str(v) for k, v in (entry.get("params") or {}).items()}
        self.steps = []
        os.makedirs(self.dir, exist_ok=True)

    # -- plumbing ---------------------------------------------------------------

    def path(self, *parts):
        return os.path.join(self.dir, *parts)

    def step_dir(self, name):
        for sub in ("inputs", "outputs"):
            os.makedirs(self.path(name, sub), exist_ok=True)
        return self.path(name)

    def copy_in(self, step, patterns, source):
        """Copy files matching patterns from source into step/inputs/, flat, as Silva does."""
        taken = []
        for pattern in patterns:
            for found in sorted(glob.glob(os.path.join(source, pattern))):
                shutil.copy2(found, self.path(step, "inputs", os.path.basename(found)))
                taken.append(os.path.basename(found))
        return taken

    def record(self, name, rc, out, err, command, note=None):
        log = self.path(f"{name}.log")
        with open(log, "w") as handle:
            handle.write(f"$ {command}\n\n--- stdout ---\n{out}\n--- stderr ---\n{err}\n")
        self.steps.append({
            "step": name, "rc": rc, "log": os.path.relpath(log, self.dir),
            "stdout": out, "stderr": err, "command": command, "note": note,
        })
        return rc

    def tool(self, name, script, argv, cwd=None, env=None):
        rc, out, err, command = self.runner.run(
            script, cwd or self.dir, env or {}, argv,
        )
        return self.record(name, rc, out, err, command)

    def node(self, name, argv=(), path_prepend=None, extra_env=None, timeout=900,
             step=None, as_name=None):
        """Run node `name`.

        `step` names the directory, for a node that runs over a second tree in the same
        case; `as_name` names the log, for a node that runs twice over the SAME directory
        and whose two logs must not overwrite each other.
        """
        node = NODES[name]
        env = dict(self.params)
        env.update(extra_env or {})
        rc, out, err, command = self.runner.run(
            os.path.join(WF, node, "script.py"), self.step_dir(step or name), env,
            argv, path_prepend, timeout,
        )
        return self.record(as_name or step or name, rc, out, err, command)

    def ok(self, name):
        for step in self.steps:
            if step["step"] == name:
                return step["rc"] == 0
        return False

    # -- staging ----------------------------------------------------------------

    def generate(self, into="g", metrics_into="m", shape=None, count=None):
        shape = shape or self.tok.get("shape", "compact")
        count = count or self.tok.get("n", "1")
        length = self.tok.get("len", "10")
        os.makedirs(self.path(into), exist_ok=True)
        rc = self.tool(into, GENERATOR,
                       [self.path(into), count, shape, length, "--manifest"])
        if rc != 0:
            return False
        if self.tok.get("metrics", "node04") == "none":
            return True
        # The real producer, not a synthesised file. design_metrics.json is
        # workflow-037 node 04's output and the whole designed_chain contract rests on
        # what it actually writes, so the suite runs it rather than writing a stand-in
        # that would agree with the spec instead of with 037.
        self.step_dir(metrics_into)
        self.copy_in(metrics_into, ["*_generated.cif", "generation_manifest.json"],
                     self.path(into))
        return self.tool(metrics_into, METRICS_NODE, [],
                         cwd=self.path(metrics_into)) == 0

    def also_run(self):
        """A second, independent tree in the same case, for a claim that needs two shapes.

        `nonstandard-residues` is the one: a sequence that is PARTLY unknown has to be
        folded with the unknowns counted, and one that is ENTIRELY unknown has nothing to
        fold and must be skipped. Those are two different inputs and so two different
        runs; the fixture says so, and this is where it happens.
        """
        shape = self.tok.get("also")
        if not shape:
            return
        if not self.generate(into="g2", metrics_into="m2", shape=shape, count="1"):
            return
        self.step_dir("n1b")
        self.copy_in("n1b", ["*_generated.cif"], self.path("g2"))
        if self.tok.get("metrics", "node04") != "none":
            self.copy_in("n1b", ["design_metrics.json"], self.path("m2", "outputs"))
        if self.node("n1", step="n1b") != 0:
            return
        self.step_dir("n2b")
        self.copy_in("n2b", job_outputs(NODES["n1"]), self.path("n1b", "outputs"))
        self.node("n2", step="n2b")

    def feed(self, name):
        """Stage one node's ./inputs/ from its parents' declared outputs."""
        self.step_dir(name)
        if name == "n1":
            self.copy_in("n1", ["*_generated.cif"], self.path("g"))
            if self.tok.get("metrics", "node04") != "none":
                self.copy_in("n1", ["design_metrics.json"], self.path("m", "outputs"))
            return
        for parent in EDGES[name]:
            source = (self.path("pred") if parent == "n3" and self.has_pred
                      else self.path(parent, "outputs"))
            self.copy_in(name, job_outputs(NODES[parent]), source)

    # -- mutations --------------------------------------------------------------

    def mutate(self, when):
        """Apply the named edit, if this is its moment. Returns a note or None."""
        name = self.tok.get("mutate")
        if not name or MUTATIONS[name][0] != when:
            return None
        note = MUTATIONS[name][1](self)
        self.steps.append({"step": f"mutate:{name}", "rc": 0, "log": None,
                           "stdout": note, "stderr": "", "command": f"mutate {name}",
                           "note": note})
        return note

    # -- stages -----------------------------------------------------------------

    has_pred = False

    def wants_interface(self):
        """Does this case ask for the complex path? The ENTRY says so, not the harness.

        PARAM_SCORE_INTERFACE lives in the fixture's `params`, so there is one place it is
        stated and no way for two to disagree. A harness that staged a complex the node
        was never told to build would be a case testing the harness.
        """
        return self.params.get("PARAM_SCORE_INTERFACE", "").strip().lower() == "true"

    def align(self):
        """Run node 01b, as Silva does on EVERY run. Returns its exit code.

        01b is a dependency of node 02 in .chiral/workflow.toml, so it runs whether or not
        an interface is being scored; on the Phase 1 path it writes a target_alignment.json
        saying it was a no-op. Skipping it when score_interface is false would mean the
        Phase 1 cases ran a DAG the workflow does not have, and would leave the no-op path
        -- the one that makes "did 01b run?" answerable at all (§8.3) -- never executed.

        On the complex path 01b also needs an alignment to read. It comes from the
        inventory node 01 just wrote, via make_target_msa.py, which writes nothing when
        there is no single agreed target -- so 01b is reached and 01b does the refusing
        (§8.13).
        """
        self.feed("nalign")
        if self.wants_interface():
            self.tool("target-msa", TARGET_MSA, [
                self.path("nalign", "inputs", "design_inventory.json"),
                self.path("nalign", "inputs", "target_msa.a3m"),
            ])
        return self.node("nalign")

    def stage_ingest(self):
        if not self.generate():
            return
        self.feed("n1")
        self.mutate("n1-inputs")
        if self.node("n1") != 0:
            return
        if self.align() != 0:
            return
        self.feed("n2")
        self.mutate("n2-inputs")
        self.node("n2")

    def stage_rerun(self):
        if not self.generate():
            return
        self.feed("n1")
        if self.node("n1") != 0:
            return
        if self.align() != 0:
            return
        self.feed("n2")
        if self.node("n2") != 0:
            return
        before = snapshot(self.path("n2", "outputs"))
        self.node("n2", as_name="n2-again")
        after = snapshot(self.path("n2", "outputs"))
        self.steps[-1]["note"] = json.dumps({
            "identical": before == after,
            "before": before, "after": after,
        })
        # Now a file from a design this run knows nothing about. 029's defect fired on a
        # directory left by an EARLIER run (4d1a523e), so the stale file is the case.
        stale = self.path("n2", "outputs", "mono_design_99.fasta")
        with open(stale, "w") as handle:
            handle.write(">design_99\nAAAAAAAAAA\n")
        stale_bytes = open(stale, "rb").read()
        self.node("n2", as_name="n2-stale")
        self.steps[-1]["note"] = json.dumps({
            "stale_survived": os.path.exists(stale)
            and open(stale, "rb").read() == stale_bytes,
        })

    def _to_predict(self):
        if not self.generate():
            return False
        self.feed("n1")
        if self.node("n1") != 0:
            return False
        if self.align() != 0:
            return False
        self.feed("n2")
        if self.node("n2") != 0:
            return False
        self.feed("n3")
        self.mutate("n3-inputs")
        return True

    def stage_predict_dry(self):
        if not self._to_predict():
            return
        # What node 02 actually declared, recorded before node 03 runs. On the complex
        # path this is how a check can tell "the multimer model was not selected" from
        # "no record ever asked for it", which an argv assertion alone cannot.
        inputs = self.path("n3", "inputs", "prediction_inputs.json")
        modes = []
        if os.path.exists(inputs):
            declared = json.load(open(inputs))
            modes = sorted({r.get("mode") for r in declared.get("inputs") or []})
        self.steps.append({
            "step": "cplx-probe", "rc": 0, "log": None,
            "stdout": json.dumps({"modes": modes}), "stderr": "",
            "command": "inspect prediction_inputs.json", "note": json.dumps(modes),
        })
        self.node("n3", argv=["--dry-run"])

    def _standin(self, body, extra_env):
        binary = self.path("bin")
        os.makedirs(binary, exist_ok=True)
        for name, text in (("colabfold_batch", body),
                           ("nvidia-smi", STANDIN_NVIDIA_SMI)):
            target = os.path.join(binary, name)
            with open(target, "w") as handle:
                handle.write(text)
            os.chmod(target, 0o755)
        return binary, extra_env

    def stage_predict_timeout(self):
        """Did the worker outlive the node? Observable on the host, not in a container.

        The worker records its own pid and the harness looks for /proc/<pid> afterwards.
        In the container that number is a pid in the container's OWN namespace, so looking
        for it on the host finds nothing whatever happened -- and the check would have
        passed for exactly the reason §8.4 says is worse than no check. Running the probe
        inside the container is no better: the container exits when node 03 returns and
        takes everything in it down, so the worker dies regardless of what node 03 did.

        So the harness records whether it could observe anything at all, and the check
        reports a gap rather than a pass when it could not. The claim is established by
        the host-interpreter run, where the pid means what it says, and that is stated in
        the check's own detail rather than left for someone to notice.
        """
        if not self._to_predict():
            return
        pidfile = self.path("worker.pid")
        binary, env = self._standin(
            STANDIN_SIGTERM, {"STANDIN_PIDFILE": pidfile},
        )
        self.node("n3", path_prepend=binary, extra_env=env, timeout=240)
        observable = not self.runner.use_docker
        alive, pid = None, None
        if os.path.exists(pidfile):
            pid = (open(pidfile).read() or "").strip()
            if pid.isdigit() and observable:
                alive = os.path.exists(f"/proc/{pid}")
                if alive:
                    # Do not leave it running whatever the verdict is; the check reads
                    # `alive`, which was sampled before this.
                    try:
                        os.kill(int(pid), 9)
                    except OSError:
                        pass
        self.steps.append({
            "step": "worker", "rc": 0, "log": None,
            "stdout": json.dumps({"pid": pid, "alive_after_node_returned": alive,
                                  "observable": observable}),
            "stderr": "", "command": (f"check /proc/{pid}" if observable else
                                      "not checked: the pid is container-namespaced"),
            "note": json.dumps({
                "pid": pid, "alive": alive, "observable": observable,
                "why": None if observable else (
                    "node 03 ran in a container, so the recorded pid is in the "
                    "container's pid namespace and means nothing on this host; and the "
                    "container's own exit would have ended the worker regardless, so "
                    "running the probe inside it would answer a different question."
                ),
            }),
        })

    def stage_predict_timeout_buffered(self):
        """A predictor killed before its buffered stdout ever reached node 03.

        Nothing to sample afterwards and nothing to observe about pids: the claim here is
        purely about the recorded reason. When stdout is empty, node 03 has to fall back to
        the predictor's own log file -- which it has already collected -- and SAY that is
        where the lines came from, because "the predictor said nothing" and "the
        predictor's output was buffered away" are different facts (§8.7).
        """
        if not self._to_predict():
            return
        binary, env = self._standin(STANDIN_SILENT_LOGFILE, {})
        self.node("n3", path_prepend=binary, extra_env=env, timeout=240)

    def stage_predict_partial(self):
        if not self._to_predict():
            return
        binary, env = self._standin(
            STANDIN_PARTIAL, {"STANDIN_FAIL_FOR": "design_1"},
        )
        self.node("n3", path_prepend=binary, extra_env=env, timeout=300)

    def stage_full(self):
        if not self.generate():
            return
        self.feed("n1")
        self.mutate("n1-inputs")
        if self.node("n1") != 0:
            return
        if self.align() != 0:
            return
        self.feed("n2")
        self.mutate("n2-inputs")
        if self.node("n2") != 0:
            return
        os.makedirs(self.path("pred"), exist_ok=True)
        self.has_pred = True
        rc = self.tool("pred", PREDICTIONS, [
            self.path("n1", "outputs", "design_inventory.json"),
            self.path("n2", "outputs", "prediction_inputs.json"),
            self.path("pred"),
        ] + (self.tok.get("mp", "").split() if self.tok.get("mp") else []))
        if rc != 0:
            return
        self.feed("n4")
        self.mutate("n4-inputs")
        if self.node("n4") != 0:
            return
        self.feed("n5")
        self.node("n5")

    def go(self):
        stage = self.tok.get("stage", "ingest")
        {
            "ingest": self.stage_ingest,
            "rerun": self.stage_rerun,
            "predict-dry": self.stage_predict_dry,
            "predict-timeout": self.stage_predict_timeout,
            "predict-timeout-buffered": self.stage_predict_timeout_buffered,
            "predict-partial": self.stage_predict_partial,
            "full": self.stage_full,
        }[stage]()
        self.also_run()
        return {
            "id": self.id, "section": self.section, "phase": self.entry.get("phase"),
            "harness": self.entry.get("harness"), "tokens": self.tok,
            "params": self.params, "dir": self.dir, "steps": self.steps,
        }


def snapshot(directory):
    """{name: sha256} for every file, so idempotence is a comparison and not a feeling."""
    import hashlib

    out = {}
    for name in sorted(os.listdir(directory)):
        full = os.path.join(directory, name)
        if os.path.isfile(full):
            out[name] = hashlib.sha256(open(full, "rb").read()).hexdigest()[:16]
    return out


# -- the named mutations ---------------------------------------------------------
#
# Each is (when, function). `when` names the moment: after a step's ./inputs/ is staged
# and before that step runs. A mutation applied anywhere else would be testing a different
# claim, so the moment is part of the definition rather than a detail of the stage code.

def _mutate_chain_c(case):
    path = case.path("n1", "inputs", "design_metrics.json")
    payload = json.load(open(path))
    for row in payload.get("designs") or []:
        row["designed_chain"] = "C"
    json.dump(payload, open(path, "w"), indent=2)
    return "design_metrics.json rewritten to name chain C, which make_designs.py never writes"


def _mutate_bad_letter(case):
    path = case.path("n2", "inputs", "design_inventory.json")
    payload = json.load(open(path))
    entry = payload["designs"][0]
    sequence = entry["binder_sequence"]
    entry["binder_sequence"] = "Z" + sequence[1:]
    json.dump(payload, open(path, "w"), indent=2)
    return (f"binder_sequence[0] of {entry['design']} changed from {sequence[0]} to Z; "
            "no generator shape can reach this because node 01 maps every unknown "
            "residue to X and X is in the alphabet")


def _mutate_drop_ca(case):
    path = case.path("n2", "inputs", "design_inventory.json")
    payload = json.load(open(path))
    entry = payload["designs"][0]
    before = len(entry["binder_ca"])
    entry["binder_ca"] = entry["binder_ca"][:-1]
    json.dump(payload, open(path, "w"), indent=2)
    return (f"binder_ca of {entry['design']} shortened from {before} to {before - 1} "
            f"while binder_sequence stayed at {len(entry['binder_sequence'])}")


def _mutate_reverse_manifest(case):
    path = case.path("n4", "inputs", "prediction_manifest.json")
    payload = json.load(open(path))
    names = [p["design"] for p in payload["predictions"]]
    payload["predictions"] = list(reversed(payload["predictions"]))
    json.dump(payload, open(path, "w"), indent=2)
    return (f"prediction_manifest.json records reversed: {names} -> "
            f"{list(reversed(names))}. A node 04 that paired scores with designs by "
            "POSITION now pairs every design with another design's numbers")


def _mutate_strip_colon(case):
    """Remove the colon from every cplx_ fasta, after node 02 wrote and verified it.

    Node 02 writes the colon itself and reads the bytes back before returning, so no
    input it is given can reach node 03 with the colon missing -- which means node 03's
    own colon guard could never be SHOWN refusing, and §8.13 says a guard that cannot be
    shown refusing is not a guard. This is the corruption that reaches it: a file changed
    between the two nodes, which is exactly what node 03's refusal says has happened.

    A missing colon is the defect that matters most here because it does not fail. It
    predicts a single chain and reports its confidence as an interface, which is 026's
    bug (§8.6).
    """
    touched = {}
    directory = case.path("n3", "inputs")
    for name in sorted(os.listdir(directory) if os.path.isdir(directory) else []):
        if not (name.startswith("cplx_") and name.endswith(".fasta")):
            continue
        path = os.path.join(directory, name)
        before = open(path).read()
        after = before.replace(":", "")
        with open(path, "w") as handle:
            handle.write(after)
        touched[name] = [before.count(":"), after.count(":")]
    if not touched:
        return ("NO cplx_ FASTA TO MUTATE in n3/inputs, so node 03 was handed nothing "
                "that could be missing a colon. The case cannot make its claim.")
    return (f"colons stripped from {len(touched)} cplx_ fasta file(s) after node 02 "
            f"verified them: {touched}")


MUTATIONS = {
    "chain-c": ("n1-inputs", _mutate_chain_c),
    "bad-letter": ("n2-inputs", _mutate_bad_letter),
    "drop-ca": ("n2-inputs", _mutate_drop_ca),
    "strip-colon": ("n3-inputs", _mutate_strip_colon),
    "reverse-manifest": ("n4-inputs", _mutate_reverse_manifest),
}


def docker_available(disabled):
    if disabled:
        return False, "disabled with --no-docker"
    if not shutil.which("docker"):
        return False, "docker is not on PATH"
    env = dict(os.environ)
    # The GPU host's own ~/.bashrc exports a DOCKER_HOST that does not exist, and a
    # non-login shell then reports no images rather than an error. Point at the rootless
    # socket when it is there.
    socket = f"/run/user/{os.getuid()}/docker.sock" if hasattr(os, "getuid") else None
    if socket and os.path.exists(socket):
        env["DOCKER_HOST"] = f"unix://{socket}"
        os.environ["DOCKER_HOST"] = env["DOCKER_HOST"]
    probe = subprocess.run(["docker", "image", "inspect", IMAGE],
                           capture_output=True, text=True, env=env)
    if probe.returncode != 0:
        return False, f"`docker image inspect {IMAGE}` failed: " + (
            (probe.stderr or "").strip().splitlines() or ["no output"]
        )[0][:160]
    return True, f"{IMAGE} present"


def main():
    argv = sys.argv[1:]

    def flag(name, default=None):
        if name in argv:
            return argv[argv.index(name) + 1]
        return default

    only = {x for x in (flag("--only") or "").split(",") if x}
    section_filter = flag("--section")
    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    root = os.path.abspath(flag("--out") or os.path.join(HERE, "fixture-runs", stamp))

    declared = check_edges()
    use_docker, why = docker_available("--no-docker" in argv)
    os.makedirs(root, exist_ok=True)

    say("=" * 78)
    if use_docker:
        say(f"RUNNING IN {IMAGE} -- this is the gate.")
    else:
        say("NOT THE GATE: running on the host interpreter, not in "
            f"{IMAGE}.\n  Reason: {why}.\n  Node behaviour that depends on the image's "
            "own library versions is NOT established by\n  this run. biopython 1.85 in "
            "particular only warns on a FASTA comment line where\n  other versions "
            "raise, which is the hazard node 05's no-';'-banner rule exists for.")
    say(f"out: {root}")
    say(f"edges: {json.dumps(declared)}")
    say("=" * 78)

    payload = json.load(open(os.path.join(HERE, "fixtures.json")))

    # An --only id that matches nothing is a typo, and a typo here is a run that tests
    # less than whoever typed it believes. The first version skipped silently: asking for
    # three cases and getting two printed no complaint and exited 0, so the missing one
    # looked like a case that had passed. Refuse instead, and name what is close.
    if only:
        known = {entry["id"] for section in ("ingest", "predict", "scoring", "rejections")
                 for entry in payload[section]} | {e["id"] for e in payload["controls"]}
        unknown = sorted(only - known)
        if unknown:
            import difflib
            lines = []
            for name in unknown:
                near = difflib.get_close_matches(name, sorted(known), n=3, cutoff=0.5)
                lines.append(f"  {name!r}" + (f" -- did you mean {near}?" if near else ""))
            raise SystemExit(
                "run_fixtures: --only names {} id(s) that are not in fixtures.json:\n{}\n"
                "Nothing was run. A silently ignored id is a case you think you tested."
                .format(len(unknown), "\n".join(lines))
            )
        gpu_only = sorted(only & {e["id"] for e in payload["controls"]})
        if gpu_only:
            say(f"note: {gpu_only} are tier 'gpu' controls and are not run here. "
                "See tests/run_controls.sh and tests/run_controls_phase2.sh.")

    runner = Runner(use_docker, root)
    results = []
    for section in ("ingest", "predict", "scoring", "rejections"):
        for entry in payload[section]:
            if entry.get("tier") != "ci":
                continue
            if only and entry["id"] not in only:
                continue
            if section_filter and section != section_filter:
                continue
            say(f"[{section}/{entry['id']}] {entry.get('harness')}")
            case = Case(runner, root, section, entry)
            try:
                results.append(case.go())
            except Exception as exc:          # a broken harness is not a failed fixture
                import traceback
                results.append({
                    "id": entry["id"], "section": section,
                    "phase": entry.get("phase"), "harness": entry.get("harness"),
                    "tokens": case.tok, "params": case.params, "dir": case.dir,
                    "steps": case.steps,
                    "harness_error": f"{type(exc).__name__}: {exc}",
                    "traceback": traceback.format_exc(),
                })
                say(f"    HARNESS ERROR: {type(exc).__name__}: {exc}")
                continue
            tail = ", ".join(
                f"{s['step']}={s['rc']}" for s in results[-1]["steps"]
                if not s["step"].startswith(("mutate:", "worker", "cplx-probe"))
            )
            say(f"    {tail}")

    run = {
        "generated": datetime.datetime.now().isoformat(timespec="seconds"),
        "root": root,
        "image": IMAGE if use_docker else None,
        "is_gate": use_docker,
        "runner": "docker" if use_docker else f"host {sys.executable}",
        "runner_reason": why,
        "edges": declared,
        "n_cases": len(results),
        "cases": results,
    }
    with open(os.path.join(root, "run.json"), "w") as handle:
        json.dump(run, handle, indent=2)

    link = os.path.join(HERE, "fixture-runs", "latest")
    try:
        if os.path.islink(link) or os.path.exists(link):
            os.remove(link)
        os.symlink(root, link)
    except OSError:
        pass

    say(f"\n{len(results)} case(s) staged into {root}")
    checker = os.path.join(HERE, "check_fixtures.py")
    if not os.path.exists(checker):
        say("check_fixtures.py is not here, so nothing has been judged.")
        return 1
    return subprocess.run([sys.executable, checker, root]).returncode


if __name__ == "__main__":
    sys.exit(main())
