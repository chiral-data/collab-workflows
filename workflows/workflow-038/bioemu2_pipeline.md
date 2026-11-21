# workflow-038 — BioEmu2 Equilibrium Ensemble Emulation (design doc)

Status: **run end to end on a GPU, with a fixture suite over seven RCSB entries**
(`tests/README.md`). This document is the design rationale and the record of decisions; it also
holds the draft `[Workflow]` issue body and PR description for whenever the work is filed.

## Context

`collab-workflows` predicts *one* structure per sequence (workflow-026 ColabFold,
workflow-012 Boltz-2) and designs new ones (workflow-013, -033, -037). Nothing in the repo
samples the **conformational ensemble** of an existing protein. The MD workflows
(-030 … -036) are polymer/materials pipelines on GROMACS, not biomolecular sampling, and they
cost hours of wall clock per system.

BioEmu2 (Microsoft Research AI for Science, public preview 2026-10-20, MIT) fills exactly that
gap: a ~500M-parameter all-atom diffusion model that emits *many* conformations per system and
can propagate them forward in time, in minutes rather than the days an equivalent MD campaign
would take. It is the first tool in this repo that answers "how does this protein *move*?"

## Scope decisions

| Decision | Choice | Consequence |
|---|---|---|
| Headline capability | **Equilibrium ensemble sampling** | The differentiator vs workflow-026. `predict-structure` is available via `mode` but is not the reason this workflow exists. |
| Mode coverage | **All three**, via `mode` + `propagate` params | `sample-equilibrium` / `predict-structure` in node 03; time-lagged propagation in node 04. |
| CUDA | **CUDA 12.6 (`[cuda12]` extra)** | See [CUDA choice](#the-cuda-choice-is-forced). Forgoes the CUDA 13 cuEquivariance kernels. |
| MSA | **Its own CPU node, writing A3M to disk** | See [why 02-fetch-msa is its own node](#why-02-fetch-msa-is-its-own-node). |
| Steering | **Exposed as `sampler`, default `edm` (unsteered)** | `edm_gd`/`edm_fk`/`edm_fkc` selectable; custom potential YAML is deliberately *not* exposed. |
| Context parallelism | **Not supported** | Needs ≥4 GPUs and a `torchrun` launcher; Silva gives one GPU per job. Systems that do not fit on one GPU are out of scope. |
| Physicality | **Reported, never filtered** | Matches upstream: all frames are kept, `physicality.json` labels them, node 05 computes the pass rate and node 06 shows it. |

## Upstream facts this design is built on

Every item here was read from upstream source or the HF API on 2026-11-11, not recalled.

- PyPI `bioemu2`, latest **`0.1.0rc0`** (released 2026-10-20) — a **pre-release**, so `pip`
  needs `--pre` or an exact `==0.1.0rc0` pin. MIT. `requires-python = ">=3.12"`;
  classifiers list 3.12/3.13/3.14, but only **3.12** has evaluated benchmark metrics.
- **Linux + NVIDIA GPU only.** macOS and Windows unsupported.
- Two mutually exclusive CUDA extras (`[tool.uv] conflicts` declares them so):
  - `cuda` → `cuequivariance-ops-cu13==0.11.1`, i.e. **CUDA 13**, compute capability **≥ 7.5**.
  - `cuda12` → `torch==2.13.0+cu126` + `cuequivariance-ops-cu12==0.11.1`, needs
    `--extra-index-url https://download.pytorch.org/whl/cu126`. Covers cc **7.0** (V100) and
    upward; upstream explicitly tested T4 (7.5) on this path.
- Core deps include `openmm`, `rdkit`, `mdtraj`, `biotite>=1.5`, `biopython>=1.88`, `numba`,
  `triton`, `torch>=2.11,<2.14`. **`openmm` is the install risk** — see open item 1.
- Entry point: `bioemu2 = "bioemu2.cli:main"`. Commands: `sample-equilibrium`,
  `predict-structure`, `propagate`, `convert-boltz`, plus `bioemu2.cp` under `torchrun`.
- Weights: `microsoft/bioemu2` on HF — `checkpoints/bioemu2-preview/checkpoint.ckpt`
  (**2,047,661,181 B ≈ 1.91 GiB**) + `config.yaml` (1672 B). Current revision
  `dd4b0c1b69984b68edcb6200a544213be4d75242` (`lastModified` 2026-10-19). Default checkpoint
  name is `bioemu2-preview`. `HF_HOME` / `HF_HUB_CACHE` relocate the cache;
  `--checkpoint-cache-dir` overrides per run.
- MSA: on a cache miss `get_msa()` **queries the public ColabFold server** and needs network.
  Cache defaults to `~/.cache/bioemu2/colabfold_msa`, relocatable with **`BIOEMU2_MSA_CACHE`**.
  Local A3M is accepted (`source=` a path, or a `{chain_id: path}` dict for heteromers); the
  first A3M record must be the query and must match the system's sequence after stripping
  gaps/insertions. Upstream asks users **not to send bulk workloads to the public server**.
- Output of every sampling command in `--out`: `samples.xtc`, `topology.cif`, `system.yaml`,
  `command.txt`, and `physicality.json` (one entry per frame: zero-based `frame`, `status` ∈
  `passed`/`rejected`/`not_checked`, `reason`). Completed batches land in `<out>/partial`,
  which the CLI removes only on full success — so an interrupted run leaves it behind.
- `save_frames()` (Python API) writes `samples.xtc` + `topology.cif`.
- Mode semantics: `predict-structure` ⇒ `step_scale=1.5` with sharpening;
  `sample-equilibrium` ⇒ `step_scale=1.0` unsharpened. Returned frames are **independent
  samples, not a time series**.
- `propagate` takes `--x0-cif` / `--x0-xtc` (one endpoint per input frame) **or** `--x0-cif`
  with `--num-samples`; exactly one of `--x0-xtc` and `--num-samples`. `--lag-time-ns`
  recommended range **1–1000 ns**.
- System YAML schema is `entities:` (each with `entity_type`, `sequence`, `asym_ids`, `mods`)
  plus `bonds:`. `asym_ids` must be unique and ≤ 4 characters.
- `--pdb-id` loads **biological assembly 1 including all non-water ligands**, and the CLI
  exposes no assembly/chain/ligand selection — that requires a system YAML.
- Model-card limits: quality degrades on shallow/poor MSAs; training data over-represents
  well-studied families; **sampling frequency is not a calibrated Boltzmann probability**;
  binding-affinity prediction out of scope; >10k tokens may not be physically plausible.
  Smallest device upstream tested is an 80 GB A100.

## The CUDA choice is forced

The `cuda` extra pulls `cuequivariance-ops-cu13`, which needs a CUDA 13 driver. This repo's
estate does not have one — `apps/applications.json` tops out at **12.6** (workflow-037's
`rfd4-proteina`), with everything else at 12.1 or below:

```
11.0 ×1   11.3 ×1   11.6 ×1   11.8 ×3   12.0 ×1   12.1 ×3   12.6 ×1
```

So `[cuda12]` it is. This is the same conclusion workflow-037 reached for RFD4-Proteina, and it
has a bonus: `cuda12` is the only path that supports compute capability 7.0, so this workflow
runs on V100 as well as T4/A100/L40S/H100. The cost is the CUDA 13 cuEquivariance kernels,
which are a speed optimisation, not a capability.

## Why `02-fetch-msa` is its own node

The MSA is the one step that needs the public internet, is the slowest non-GPU step, and is the
most likely thing to fail for reasons that have nothing to do with the model — a rate limit, a
server outage, a sequence the server cannot align. Folding it into node 03 would mean paying for
a GPU to sit idle while an HTTP request to `api.colabfold.com` retries, and would turn a network
failure into a failed GPU job.

Splitting it has three concrete payoffs:

1. **Failure is cheap and legible.** A bad sequence or an unreachable server fails in a CPU node
   in seconds.
2. **Node 03 becomes network-free.** Node 02 writes `msa/<chain>.a3m` into its outputs, node 03
   loads those files, and the checkpoint is already baked into the image — so node 03 runs with
   `HF_HUB_OFFLINE=1` and never reaches for the network at all.
3. **The alignment becomes an inspectable, reusable artifact.** `msa/` and
   `msa_summary.json` (depth per chain) are exactly what you want to look at when an ensemble
   comes out suspiciously narrow, because shallow MSA is upstream's stated first cause of
   degraded quality. The user can also drop their own A3M into the inputs and skip the server
   entirely — which is also the polite thing to do for repeated runs, per upstream's request
   not to send bulk workloads to the public server.

## Why propagation is a separate node that always writes output

Propagation doubles GPU cost and most runs will not want it, but Silva has no notion of a
conditionally-skipped node, and a declared output that is not written is pitfall #5 in
`SILVA_MIGRATION_GUIDE.md` §9 (silent failure).

**Decision:** node 04 always writes its declared outputs. When `propagate = false` it copies the
node 03 ensemble to `propagated_samples.xtc` / `propagated_topology.cif` and records
`{"enabled": false, "lag_time_ns": 0.0}` in `propagation_manifest.json` — an identity
propagation. Node 05 reads that manifest and, when `enabled` is false, reports the ensemble as
equilibrium samples only and emits no time-lagged section; node 06 does the same. Nothing
downstream ever presents identity frames as propagated. This keeps the DAG static and the
declared outputs honest at the cost of one file copy.

## Why nodes 03–05 use the Python API, and never parse `topology.cif`

Two decisions that only became clear from reading `bioemu2/cli.py` and `inference.py`, after
the rest of this design was settled:

**The CLI cannot accept a pre-computed MSA.** `_common_options` in `cli.py` offers
`--system`/`--pdb-id`/`--checkpoint`/`--checkpoint-path`/`--checkpoint-cache-dir`/`--out`/
`--batch-size`/`--seed`/`--num-steps`/`--sampler`/`--potentials`/`--physicality-filter` plus
three speed flags — and nothing for an alignment. Node 02 exists precisely so node 03 needs no
network, and only `EmuSession(system, model, mode, msa=...)` can take the result. So nodes 03
and 04 drive the Python API instead. A bonus: `save_frames()` writes exactly `samples.xtc`,
`topology.cif` and `physicality.json` into the directory it is handed, so node 03 points it at
`./outputs/` and copies nothing.

**`--num-steps` does have a Python equivalent, and this design originally said it did not.**
`EmuSession.sample()` takes `num_samples`, `batch_size`, `sampler`, `seed`, `progress` and
`partial_dir`, and the step count does belong to the sampler object — but `sampler` is typed
`SamplerPreset | Sampler`, so a configured sampler can be passed straight in, and
`sampler_from_preset(preset, mode=…, num_steps=…)` builds one with the mode's step scale
already applied. The step count is exposed as `num_steps`, defaulting to 0, which means the
preset's own default of 200. Only a non-zero value takes the sampler-object path, so a default
run stays on the string-preset API.

The lesson is narrower than "read the source": the original claim was drawn from
`sample()`'s parameter *names* without reading its parameter *types*, and the type was the
answer.

**The topology comes from the System, not from the CIF.** `samples.xtc` carries coordinates
only, and MDTraj's PDBx/mmCIF support is the weak link in that chain. But `System` already
exposes `system.topology` as an MDTraj topology — it is what `map_cif_atom_order()` itself
builds its trajectories against — so nodes 04 and 05 bind the ensemble with
`mdtraj.load_xtc(path, top=system.topology)`. Atom order then matches by construction on both
the sampling and the analysis side, and no mmCIF is parsed anywhere in this workflow.
`topology.cif` stays a declared output because it is the portable artifact for everything
*outside* the workflow.

## Node structure

| Node | Name | Inputs | Outputs | GPU |
|------|------|--------|---------|-----|
| 01 | `01-prepare-system` | (params; optional upload) | `system.yaml`, `system_inventory.json` | No |
| 02 | `02-fetch-msa` | `system.yaml`, `system_inventory.json` | `msa/` (unpaired + paired A3M), `msa_summary.json`, `system.yaml`, `system_inventory.json` | No |
| 03 | `03-sample` | `system.yaml`, `msa_summary.json`, `msa/` | `samples.xtc`, `topology.cif`, `physicality.json`, `system.yaml`, `sampling_manifest.json`, `msa/` | **Yes** |
| 04 | `04-propagate` | node 03 ensemble, `physicality.json`, `sampling_manifest.json`, `msa/` | `propagated_samples.xtc`, `propagated_topology.cif`, `propagated_physicality.json`, `propagation_manifest.json` | **Yes** |
| 05 | `05-analyze` | nodes 03 + 04, minus `topology.cif` which nothing reads | `ensemble_metrics.csv`, `ensemble_summary.json`, `ensemble_top.pdb` | No |
| 06 | `06-report` | nodes 03 + 04 + 05 | `report.html` | No |

```text
params / uploaded system or structure
   │
   ▼
[01-prepare-system] ──► system.yaml, system_inventory.json
   │
   ▼
[02-fetch-msa]  (network) ──► msa/, msa_summary.json, system.yaml
   │
   ▼
[03-sample]     (GPU) ──► samples.xtc, topology.cif, physicality.json, msa/
   │                      sampling_manifest.json, system.yaml
   ├───────────────┐
   │               ▼
   │        [04-propagate] (GPU, identity when disabled)
   │               │   └──► propagated_samples.xtc, propagated_physicality.json,
   │               │        propagation_manifest.json
   ├───────────────┼──────────► [05-analyze] ──► ensemble_metrics.csv,
   │               │                │            ensemble_summary.json, ensemble_top.pdb
   └───────────────┴────────────────┴──────────► [06-report] ──► report.html
```

Six nodes rather than five: `02-fetch-msa` earns its place for the reasons above, and analysis
is kept out of the report node so the metrics table is a reusable artifact rather than a
side effect of rendering HTML.

## Conventions followed

Parameters are declared in `.chiral/workflow.toml` and **duplicated verbatim into the
`job.toml` of each node that consumes them** (the repo's established pattern — Silva auto-maps
`[params.x]` → `PARAM_X`; no `env` field, no `depends_on`).

Scripts follow workflow-026, the template-conformant reference. `run.sh` is a 5-line bash
entrypoint with `set -e` that calls `python script.py "$@"`; all parameter reading happens in
Python. Every `script.py`: `os.makedirs("./outputs", exist_ok=True)` first;
`os.environ.get("PARAM_X", "<default identical to job.toml>")`; reads literal
`./inputs/<name>` or `glob.glob("./inputs/<pattern>")`; `print(..., flush=True)`;
`sys.exit(1)` after an `ERROR:`-prefixed message; no argparse, no absolute paths, no runtime
`pip install`.

Both GPU nodes support `--dry-run`, which validates every parameter, prints the fully resolved
settings as JSON, and exits 0 without loading the model or touching the GPU — the testability
pattern workflow-026's `02_predict` uses. Every parameter is validated against an explicit
allowed set before the model is loaded, so a typo in `sampler` or `mode` fails in milliseconds
rather than after a trunk embedding.

Both also set `HF_HUB_OFFLINE=1` **before importing bioemu2**, so a mis-built image turns a
silent multi-gigabyte download into an immediate, legible error.

Image name `bioemu2:2026_11_20`. The repo's shared-image CI
(`build-changed-shared-workflow-images.yml`) greps the tag out of `job.toml` with
`grep -A20 '\[container\]' | grep '^image'`, so **`image = "..."` must sit at column 0 within
20 lines of `[container]`** in every `job.toml`.

## Open items

1. **Closed: openmm installs from PyPI, no conda base image needed.** bioemu2 requires it
   unpinned, and PyPI publishes cp312 manylinux wheels. manylinux_2_28 needs glibc >= 2.28 and
   `python:3.12-slim` (Debian 12) ships 2.36. Reproduce without building anything:
   `pip download openmm --no-deps --python-version 312 --only-binary=:all: --platform
   manylinux_2_28_x86_64`. The `continuumio/miniconda3` fallback is dropped from the Dockerfile.
   Measured in the built image (2026-11-20): openmm 8.6.1, mdtraj 1.11.2, torch 2.13.0+cu126,
   numpy 2.4.6, biotite 2.0.0, bioemu2 0.1.0rc0.
2. **Verify `0.1.0rc0` is still the newest pre-release** at build time
   (`pip index versions bioemu2 --pre`) and re-pin if not. It is a *preview*; the API may move.
3. **Closed: `physicality.json` is a bare JSON list.** `PhysicalityReport.save()` in
   `bioemu2/physicality.py` writes `[{"frame": i, "status": ..., "reason": ...}]` with status in
   `passed` / `rejected` / `not_checked`, and `not_checked` whenever `checked` is false or the
   mode is `no_filter`. The `{"frames": [...]}` envelope node 05 also accepted does not exist,
   so that branch was dead code tolerating a format upstream never writes. Node 05 now parses
   the real shape and fails on anything else instead of degrading to blank statuses, which is
   what let a dropped file look like a run where no check was requested.
4. **Half closed: the compute-capability floor is measured, the VRAM ceiling is not.**
   `torch.cuda.get_arch_list()` in the built image returns
   `['sm_50', 'sm_60', 'sm_70', 'sm_75', 'sm_80', 'sm_86', 'sm_90']`, so `sm_70` is present —
   the 7.0 floor and the V100 claim are real rather than inherited — and `sm_100` is absent,
   which is the B200 exclusion. The pipeline has been run end to end on an 8 GB RTX 3060 Ti
   (sm_86, driver 580.178.04, 7.66 GiB usable), 8 samples in one batch: 1UBQ (76 tokens,
   31.2 s), 5PTI (65), 2ZTA (72), 1BNA (69) and 1CLL (155) all complete; **1AKE (542 tokens)
   and 1A3N (746) run out of memory**, with node 03's own OOM message and advice. The useful
   finding is the gap between that and upstream's estimate, which node 03 prints: it gave 33
   for 542 tokens and 19 for 746, both at least four times optimistic on this card. The formula
   is linear in tokens and upstream's tested device is an 80 GB A100, so node 03, the README and
   the troubleshooting section now all say the estimate is an upper bound rather than a ceiling.
   **Now closed, with the opposite answer to the one expected.** Re-run at `BATCH_SIZE=1`, both
   systems failed on the *same* allocation as at batch 8 — 1.06 GiB for 542 tokens, 1.46 GiB for
   746 — so on this card the batch is not what the memory is going on: the peak is set by the
   token count. Upstream's formula is not merely optimistic here, it scales the variable that
   does not bind. Node 03's OOM message now says to try a lower batch but not to expect much,
   and names the token count as what to come down on; the README states the ceiling as a size
   range (8 GB holds somewhere above 155 and below 542 tokens) rather than as a batch to tune.
   `requirements.memory_gb: 16` in `apps/applications.json` remains an estimate and is labelled
   as one — it is a VRAM-adjacent guess for the host, and nothing here has measured it.

9. **The root README's numbering depends on merge order.** The list jumps from 36 to 38,
   because workflow-037 is on the `workflow/rfd4-proteina` branch and not yet merged. Markdown
   renumbers an ordered list sequentially, so until 037 lands this entry renders as "37". This
   is a merge-order note for the reviewer, not something to work around in the source — the
   numbers come right on their own once both branches are in.
5. **No resume.** `EmuSession.sample()` only writes resumable batches when it is given
   `partial_dir`, and node 03 does not pass one — so nothing pollutes `./outputs/`, but a
   pre-empted production run restarts from zero. If large runs start getting pre-empted, pass a
   `partial_dir` outside `./outputs/` and have Silva's retry reuse it; note that resuming
   changes the RNG story (batch *i* uses `seed + i` with a seed, new randomness without one).
6. **Trp-cage MSA depth.** The 20-residue test system may return a very shallow alignment,
   which is fine for a smoke test but makes the ensemble scientifically uninteresting. The GPU
   campaign therefore runs seven systems, not one: 1UBQ (ubiquitin, the well-behaved control),
   1AKE (adenylate kinase, upstream's own example and a genuine two-state system), 5PTI (BPTI,
   three disulfides, so the `bad_bonds_sulfur` check is actually exercised), 2ZTA (the GCN4
   leucine zipper, a homodimer, so the two-chain alignment and the per-chain RMSF split are
   exercised), 1CLL (calmodulin, four structural calcium ions, every one an atom named `CA`),
   1A3N (haemoglobin, four chains over two sequences, so a **paired** alignment is due) and
   1BNA (B-DNA, no protein chain, so no alignment is due and none must be demanded).
   `tests/README.md` says what each one catches. Trp-cage stays as the bundled smoke-test
   input.
7. **A user-supplied A3M cannot reach node 02 through Silva.** Node 02 reads `./inputs/*.a3m`
   and uses it verbatim, which works when the node is run directly with Docker — but Silva's
   only file-typed parameter (`input_file`) routes to node 01, and node 02's inputs come from
   node 01's outputs. So the documented "bring your own alignment" path is Docker-only today.
   Fix by either adding a second `file` param bound to node 02, or having node 01 pass uploaded
   A3M files through to its outputs; the second needs open item 8 resolved first. The READMEs
   now state the limitation rather than promising the feature.
10. **Closed: every non-standard component needs the network, from inside `system.topology`.**
   An ion, a cofactor, a terminal cap or a modified residue is built from its CCD definition,
   which bioemu2 downloads from `files.rcsb.org` the first time it meets the code
   (`bioemu2/structure/templates.py`, `_CCD_URL`), caching under `BIOEMU2_CCD_CACHE` and
   `BIOEMU2_TEMPLATE_CACHE`. That call is reached from `system.topology`, which nodes 01, 03, 04
   and 05 all make, so an air-gapped run of any holo system failed four frames inside the
   library with a DNS error. Found by running the fixture suite with `--network none`: five of
   the seven entries hit it. The image now pre-warms 126 codes (~2 MB), each node names the
   component and the cache when a fetch fails, and node 01 fails rather than warning, because
   the nodes after it build the same structure. Two of the names were themselves findings:
   5PTI's `DOD` (deuterated water — upstream omits waters written `HOH`, not those written
   `DOD`) and its `UNX` (density nobody identified). `UNL` cannot be built at all; node 01
   lists the components and names `exclude_ligands`.

8. **Closed, and it was the other half of the question that mattered.** An *output* glob
   matching nothing is harmless: `collect_output_files` runs its patterns under
   `shopt -s nullglob`, logs "matched no files", and exits 0. The fatal half was the *input*
   pattern. Silva matches an `inputs` entry against the **basename only**, so `msa/*.a3m` could
   never match anything, the directory was never staged, and node 03 sampled with no alignment
   while reporting the one node 02 fetched. Both sides now declare the directory name `msa`,
   following merged workflow-029's `inputs = ["backbones"]`. `silva validate` names this class
   of error directly, which is the lesson: run it before anything else.


## Verification plan

Per `SILVA_MIGRATION_GUIDE.md` §6 — all must pass before this workflow is considered done.
**Run on 2026-11-20** against the built image on an RTX 3060 Ti, plus the offline fixture
suite in `tests/` over seven RCSB entries. `tests/README.md` records what passed, what each
fixture catches, and what the suite still cannot establish.

Build-time checks (these also resolve open items 1–2):

```bash
docker build -t bioemu2:2026_11_20 workflows/workflow-038/
docker run --rm bioemu2:2026_11_20 bioemu2 --help
docker run --rm bioemu2:2026_11_20 python -c "import bioemu2, openmm, mdtraj, biotite; print('OK')"
docker run --rm bioemu2:2026_11_20 ls -l /models/hub/models--microsoft--bioemu2
docker run --rm --gpus all bioemu2:2026_11_20 \
  python -c "import torch; print(torch.cuda.is_available(), torch.version.cuda)"
docker run --rm --gpus all bioemu2:2026_11_20 \
  nvidia-smi --query-gpu=name,compute_cap --format=csv
```

Compute capability must be ≥ 7.0. Below that the workflow cannot run at all.

Per-node loop (`AI_WORKFLOW_GUIDE.md` §5), node 03 shown:

```bash
mkdir -p /tmp/wf038/inputs /tmp/wf038/outputs
cp -r /tmp/wf038-node02/outputs/* /tmp/wf038/inputs/
docker run --rm --gpus all \
  -v /tmp/wf038/inputs:/workspace/inputs \
  -v /tmp/wf038/outputs:/workspace/outputs \
  -v $(pwd)/workflows/workflow-038/03-sample:/workspace \
  -w /workspace \
  -e PARAM_MODE=equilibrium -e PARAM_NUM_SAMPLES=8 -e PARAM_SEED=42 \
  bioemu2:2026_11_20 bash run.sh --dry-run    # then again without --dry-run
```

Chain + end-to-end: `cp -r /tmp/wf038-nodeN/outputs/* /tmp/wf038-nodeN+1/inputs/` through all
six nodes, confirm `report.html` renders with an animating 3Dmol ensemble, then run under Silva
with `export SILVA_WORKFLOW_HOME=…` and capture the green-checkmark screenshot the §6 checklist
requires.

**Correctness gates**, in increasing strength:

1. `samples.xtc` holds exactly `num_samples` frames and `physicality.json` has one entry per
   frame, with `frame` zero-based and contiguous.
2. The ensemble has not collapsed. "Non-zero spread" is **not** the test: two independently
   sampled frames always differ, so that condition is met by any output whatsoever. For a
   folded protein in `equilibrium` mode the maximum pairwise C-α RMSD has to clear roughly
   0.1 nm, the scale of ordinary backbone fluctuation; below that the ensemble is one structure
   plus sampling noise. Node 05's own collapse warning fires at 1e-3 nm, which catches only the
   degenerate case, so this gate is the stricter one and is checked by eye against
   `ensemble_summary.json` until the first run gives a measured per-system baseline.
3. `mode=structure` with `num_samples=1` on Trp-cage (1L2Y) returns a compact frame. Whether
   it passes the physicality checks is to be **measured, not asserted**: nothing here has run
   on a GPU yet, and a 20-residue miniprotein with a single tryptophan cage is not obviously an
   easy case for a bond-geometry check. Record the verdict and the first failed check if any,
   and only then decide whether a pass belongs in this gate.
4. The ground truth worth having is upstream's own example: **adenylate kinase (1AKE)** in
   equilibrium mode should sample both the open and closed lobe conformations. If the ensemble
   shows only one, sampling is broken in a way no smaller test would reveal.

   **Not achieved on the available hardware.** 1AKE is 542 tokens and does not fit the 8 GB
   reference card at any batch size, so no 1AKE ensemble exists to judge. The substitute the
   campaign actually measured is **calmodulin (1CLL)**, 155 tokens: Rg 1.5677 nm with max
   pairwise C-α RMSD 1.3682 nm over 8 clusters, against ubiquitin's 1.1597 nm over a single
   0.2204 nm cluster. A two-lobed protein that pivots about a flexible central helix coming out
   multi-state while a rigid fold comes out single-state is the same kind of evidence 1AKE was
   wanted for, one system smaller. It is a weaker control — nothing here identifies the clusters
   as the open and closed states — and the 1AKE check stays worth running on a larger card.

## Scientific caveats this workflow must not paper over

BioEmu2 is a **research preview emulator, not a simulator**, and the report says so:

- **Sampling frequency is not a calibrated Boltzmann probability.** Frame counts are not free
  energies. The report must not present a population fraction as a ΔG, and node 05 computes no
  thermodynamic quantity.
- **Frames are independent samples, not a trajectory.** Only node 04's output is time-related,
  and then only as a lagged endpoint per start frame — not a continuous path.
- **Quality tracks MSA depth.** `msa_summary.json` carries the depth per chain, and the report
  surfaces it, so a narrow ensemble can be attributed rather than believed.
- **Physicality is reported, not enforced.** Rejected frames remain in the XTC by design; the
  report shows the pass rate and marks rejected frames in the per-frame table.
- Microsoft's model card calls it research software and not a validated commercial tool.

No biosafety screening question arises: this workflow samples conformations of a sequence the
user already has. It designs nothing new.

## Process note

`SILVA_MIGRATION_GUIDE.md` §8 says to open a `[Workflow]` issue and wait for approval from
Rongjun before starting implementation. Rongjun left on 2026-10-19 and the repo now has a sole
maintainer, so that gate no longer applies. The drafts below are here to file whenever desired.

### Draft issue body

> **[Workflow] BioEmu2 equilibrium ensemble emulation**
>
> Add a conformational-ensemble workflow built on BioEmu2 (Microsoft Research AI for Science,
> MIT, public preview 2026-10-20). The repo can predict a single structure and design new ones,
> but has nothing that samples how an existing protein moves; BioEmu2 emits ensembles in minutes
> rather than the days an equivalent MD campaign costs.
>
> Six nodes: prepare system → fetch MSA (CPU, network) → sample (GPU) → propagate (GPU,
> optional) → analyze → report. Covers all three BioEmu2 modes — `predict-structure`,
> `sample-equilibrium`, and time-lagged `propagate` — selected by parameter.
>
> CUDA 12.6 via the `[cuda12]` extra: the `[cuda]` extra needs a CUDA 13 driver and the estate
> tops out at 12.6. `[cuda12]` is also the only path that supports compute capability 7.0, so
> this runs on V100 as well as T4/A100/L40S/H100. Context-parallel (≥4 GPU) sampling is out of
> scope.

### Draft PR description

> Adds `workflows/workflow-038/` — BioEmu2 equilibrium ensemble emulation, six nodes on one
> shared image (`bioemu2:2026_11_20`, `python:3.12-slim` base, the 1.91 GiB `bioemu2-preview`
> checkpoint baked at a pinned HF revision, plus the CCD definitions for 126 common components).
> Nothing in the repo samples how an existing protein *moves*; this does.
>
> - `01-prepare-system` — resolve a system from an upload, a structure file, a raw sequence or an
>   RCSB accession; inventory chains, tokens and every non-polymer component by CCD code
> - `02-fetch-msa` — ColabFold alignment to disk as A3M, unpaired and, for a heteromer, paired,
>   so node 03 needs no network
> - `03-sample` — `equilibrium` or `structure` mode on pinned offline weights, with a
>   `--dry-run` that checks the staging, the chemistry and the alignment's identity on a CPU
> - `04-propagate` — time-lagged endpoints; an identity pass-through when disabled, so the
>   declared outputs always exist
> - `05-analyze` — per-frame and ensemble geometry, clusters, per-chain RMSF, and a multi-model
>   PDB of representative frames
> - `06-report` — one self-contained `report.html` with an animating 3Dmol ensemble
>
> **Verified on a GPU, not only reviewed.** On an RTX 3060 Ti (8 GB, sm_86): 1UBQ samples 8
> frames in 31.2 s with 8/8 physical and Rg 1.1597 nm against an experimental 1.17 nm; 5PTI,
> 2ZTA, 1CLL and 1BNA (no protein chain) run through all six nodes; and the time-lagged path
> runs for real on 1UBQ at a 1 ns lag, giving 8 physical endpoints at Rg 1.1628 nm against
> 1.1597 before it. 1CLL is the ensemble worth looking at: max pairwise C-α RMSD 1.3682 nm over
> 8 clusters on a 148-C-α selection that correctly excludes its four calcium ions, against
> ubiquitin's single 0.2204 nm cluster — a two-lobed protein coming out multi-state and a rigid
> one coming out single-state.
>
> 1AKE (542 tokens) and 1A3N (746) do not fit the card, at batch 8 or at batch 1, both times on
> the same allocation — so the ceiling is the token count rather than the batch, which is what
> the README now says and what node 03's OOM message now advises. The image's compute-capability
> floor is measured rather than inherited: `torch.cuda.get_arch_list()` returns `sm_50 … sm_90`,
> so V100 is in and B200 is out.
>
> `tests/` holds a fixture suite over seven RCSB entries, chosen by measurement — `probe_entry.py`
> parses `_atom_site` by hand so it shares no parser with the code under test — split into an
> offline driver and the GPU campaign. 1CLL is in it because a calcium ion is an atom named `CA`
> and mdtraj's `name CA` matches atom names, which silently puts 148 atoms in a 144-C-α
> selection. The offline driver runs every node with `--network none`, which is what turned a
> class of silent failures into named ones: it is how the CCD network dependency, the
> `n_tokens: null` on every holo system, and 5PTI's `DOD` and `UNX` components were found.
> `tests/README.md` says what each entry catches, what the suite cannot establish, and the one
> gap that stays open — no heteromer finishes the pipeline on an 8 GB card.
>
> Registry: root `README.md` entry #38, `workflows/CHANGELOG.md`, `apps/applications.json`.
>
> Sampling frequencies are **not** calibrated Boltzmann probabilities and frames are independent
> samples rather than a trajectory; cluster sizes are frame counts, not populations. The README
> and the report both lead with that.

## References

- [microsoft/bioemu2](https://github.com/microsoft/bioemu2) — MIT. README, `docs/`, `cli.py`,
  `pyproject.toml` read at 2026-11-11.
- [microsoft/bioemu2 model card](https://huggingface.co/microsoft/bioemu2) — revision
  `dd4b0c1b69984b68edcb6200a544213be4d75242`, checkpoint 2,047,661,181 B.
- [bioemu2 on PyPI](https://pypi.org/project/bioemu2/) — `0.1.0rc0`, 2026-10-20.
- BioEmu (v1): *Scalable emulation of protein equilibrium ensembles with generative deep
  learning* — Microsoft Research AI for Science. <!-- TODO: BioEmu2 preprint/DOI once published -->
