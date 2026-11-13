# workflow-037 — RFD4-Proteina De Novo Design (design doc)

Status: **image builds and node 01 runs in it; generation not yet exercised.** Built on the GPU
node amd-c12-1 (RTX 3060 Ti, 8 GB) as `rfd4-proteina:2026_11_04`, 10.1 GB. The CPU nodes were
developed against local Python 3.12; node 01 is confirmed in-image, node 02's validation path was
blocked by a missing system library (now fixed, rerun pending), and node 03 has only run
`--dry-run`. This document is the design rationale and the record
of decisions; it also holds the draft `[Workflow]` issue body and PR description for whenever the
work is filed.

## Context

`collab-workflows` has binder design (workflow-013 BoltzGen; workflow-033 RFdiffusion →
ProteinMPNN → Boltz2) and structure prediction (workflow-026 ColabFold), but nothing on the
RFD4 line. RFD4-Proteina (RosettaCommons + NVIDIA GenAIR, 2026-10-23, BSD-3) supersedes
RFdiffusion 1/2/3: one model covering monomers, binders, motif scaffolding, antibodies and
enzymes, with rich conditioning.

## Scope decisions

| Decision | Choice | Consequence |
|---|---|---|
| Scoring | **Generation only** (`~metrics`) | No AF2/OpenFold3 oracle env; small image, one GPU node. Refolding is a follow-up workflow chaining into workflow-026. |
| Task scope | **Multi-task** via `task` param | binder / motif_scaffold / monomer / binder_redesign, each a condition-spec template. |
| Inputs | **Params → generated YAML** | User supplies a structure (upload or PDB ID) + guided params; node 02 renders the condition spec. |
| CUDA | **CUDA 12** | Repo estate tops out at 12.1 (`grep cuda_version apps/applications.json`); no CUDA 13 anywhere. Uses `[all-cuda12]` + `--torch-backend cu126`. Forgoes the BioIR backend, which is cu13-only. |
| Screening | **No screening node** | Decided explicitly. See [Biosafety](#biosafety). |

## Upstream facts this design is built on

- PyPI `rfd4-proteina`; **Linux x86_64 + Python 3.12 only**; CLI `rfd4-generate`.
- Upstream's own Dockerfile uses a **slim Python base, not a CUDA base image**, with
  `RFD4_ACCELERATOR` ∈ `cu13` (default) / `cu12` / `cpu` — so `python:3.12-slim` is the right base
  and matches existing repo Dockerfiles.
- CUDA 12 install: `uv pip install --no-sources --torch-backend cu126 'rfd4-proteina[all-cuda12]'`.
- Weights: `nvidia/RFD4-Proteina` on HF, **~1.74 GB**, flow + autoencoder as a *matched pair*
  (`release-2026-10-19`), pinned revision `a06e2509c7033b010fc43614e9fd840303588502`. `HF_HOME`
  controls the cache.
- Input is a YAML/JSON **condition spec**: `input`, `contig`, `conditions` (`C_CRDGRP`, `C_SEQ`,
  `C_HOT`, `C_CRD`) with `sel('chain/res_name/res_id/atom_name')` selectors.
- Output lands in `<out_dir>/inference_samples/` as `*_generated.cif`.
- Supported GPUs: A100 / L40S / RTX6000Ada / H100 / H200 / Blackwell.
- **The CLI and Python API do not run SecureDNA screening** — only the HF Spaces apps do.

## Why `02-build-spec` is its own node

`sel()` selector and `contig` errors are the most common RFD4 failure mode, and the docs give an
explicit escape hatch: `scripts/write_condition_spec_cifs.py` renders a spec to an annotated CIF
without running the model. Running that as a separate CPU node turns a 20-minute GPU failure into a
10-second CPU failure, and `spec_preview.cif` is itself a useful artifact — the user can open it and
see exactly which atoms are conditioned.

## The numbering footgun `01-acquire-target` exists to defuse

The docs are explicit: `res_id` in `contig` and `conditions` refers to the original mmCIF
`label_seq_id`, **not** author numbering. RCSB PDB files are author-numbered, so a user typing
"hotspot A45" from a paper will usually be wrong. Node 01 emits both numberings in
`target_inventory.json`, node 02 accepts `residue_numbering` (`auth` | `label`) and translates, and
the report shows which was used.

## Conventions followed

Parameters are declared in `.chiral/workflow.toml` and **duplicated verbatim into the `job.toml` of
each node that consumes them** (the repo's established pattern — Silva auto-maps `[params.x]` →
`PARAM_X`; no `env` field, no `depends_on`).

Scripts follow workflow-026, the template-conformant reference. `run.sh` is a 5-line bash entrypoint
with `set -e` that calls `python script.py`; all parameter reading happens in Python. Every
`script.py`: `os.makedirs("./outputs", exist_ok=True)` first; `os.environ.get("PARAM_X", "<default
identical to job.toml>")`; reads literal `./inputs/<name>` or `glob.glob("./inputs/<pattern>")`;
`print(..., flush=True)`; `sys.exit(1)` after an `ERROR:`-prefixed message; no argparse, no absolute
paths, no runtime `pip install`.

Deliberately **not** copied from workflow-013: `cp inputs/* .`, bash-side param parsing, argparse
CLI flags, writing to cwd, and `pip install -r requirements.txt` at runtime (pitfall #4 in
`SILVA_MIGRATION_GUIDE.md` §9).

Node 03 builds the `rfd4-generate` argv as a **Python list passed to `subprocess.run` without
`shell=True`**. The `~metrics` override is a shell-expansion hazard (upstream docs quote it
specifically to stop `~` becoming `$HOME`); bypassing the shell removes that class of bug entirely.
It also sets `HF_HUB_OFFLINE=1`, so a baked-in-weights image never silently reaches for the network.

Image name `rfd4-proteina:2026_11_04`. The repo's shared-image CI
(`build-changed-shared-workflow-images.yml`) greps the tag out of `job.toml` with
`grep -A20 '\[container\]' | grep '^image'`, so **`image = "..."` must sit at column 0 within 20
lines of `[container]`** in every `job.toml`.

## What has and has not been verified

Verified on the development host with local Python 3.12 (Biopython, jinja2, numpy, requests):

- Node 01 on all three paths — uploaded mmCIF, uploaded PDB, and RCSB fetch of 5O45. Confirms the footgun this design
  is built around: on 5O45, chain A is author 17–145 and label 1–129, so author 45 is label 29.
- Node 02 for all four tasks, the `contig` override, both numbering modes, and every validation
  failure (unknown task, bad length, out-of-range residue, missing motif, redesign == target).
  Rendered specs parse as YAML.
- Node 03 `--dry-run` (argv resolution), missing-spec and bad-replicate failures, and that the real
  path fails cleanly with "No GPU available" on a CPU host.
- Node 04 against stand-in designs built by reordering 5O45's chain B ahead of chain A, and node 05
  over its metrics; the report's inline JS parses under node, and the chain test 01 → 02 → 03(dry) →
  04 → 05 completes end to end.

Not verified, and needing the image:

- Everything in the Dockerfile: the `[all-cuda12]` extra resolving, the version pin (still absent),
  the weight download, cu126 against the host driver.
- `rfd4-generate`'s real CLI surface. Node 03's overrides (`--config`, `out_dir`,
  `inference.num_replicates`, `inference.seed`, `~metrics`, `compile_lengths`) are a documented
  guess, built in one function (`build_argv`) so a correction is one edit.
- The `sel()` wildcard spelling. Node 02 emits `sel('A//38,40/')` for a chain-and-residues
  selection, with empty fields as wildcards per the documented `chain/res_name/res_id/atom_name`
  order. If upstream wants `*` or omitted fields instead, the fix is in `sel()`.
- Whether `inference_samples/*.cif` is really the output layout (node 03 warns and falls back).
- Whether the designed chain is really first in the output (node 04's documented assumption; it
  reports all chains so a mismatch is visible).
- The upstream `write_condition_spec_cifs.py` path — node 02 probes three locations and
  `RFD4_SPEC_CIF_SCRIPT`, and warns rather than failing when it finds none.

### Fixed during review

- **Node 01, PDB path (silent wrongness).** Biopython's mmCIF writer opens a new `label_asym_id` at
  each non-standard residue, so 5O45's chain B came out as 14 separate one- and two-residue chains.
  `target_chain=B` would then resolve to a one-residue fragment and the contig would design against
  it without complaint. Node 01 now canonicalises label numbering after conversion: one
  `label_asym_id` per author chain, `label_seq_id` sequential over standard residues.
- **Node 02, empty `loop_`.** A spec with no residue-level conditions (binder without hotspots) made
  `spec_preview.cif` carry a `loop_` with zero rows, which is invalid mmCIF — Biopython tolerates it,
  stricter parsers do not. It now writes a comment instead, and whole-chain conditions such as
  `C_CRD sel('A///')` are recorded with `label_seq_id .` rather than being dropped.
- **Node 02, unresolved gaps.** A contig span like `B1-14` over a chain with only 8 resolved residues
  now warns and names the holes, at CPU time.
- Node 01 also warns on repeated author numbers (insertion codes) instead of silently overwriting
  the translation map; node 03 keeps its scratch directory on failure and names it; node 04 computes
  Rg once; node 05's job.toml no longer declares an input it never reads.

## Open items

Most of the original list was settled by reading the real package and building the image on the
GPU node (2026-11-12/13). What the pre-build guesses got wrong is recorded in
[Corrected after reading the package](#corrected-after-reading-the-package).

1. **VRAM is the live question.** The available card is an RTX 3060 Ti, 8 GB, compute 8.6 — Ampere,
   so the kernels load, but far below upstream's smallest tested device (a 40 GB A100) and not on
   their supported list. Whether a 50-120 residue binder against a 129-residue target fits at
   `num_replicates=1` is unknown until it runs.
2. **`monomer` works, but its output quality is the open question.** The scaffold is now generated
   (see below), so the task runs; the one test design was the requested length but extended and
   poly-glutamate-tailed. Whether that is the sampler at defaults, the absence of any scoring, or
   simply one sample of one, is unknown. Upstream's own recipe uses `num_replicates=16` with
   refolding metrics, which this workflow deliberately does not have.
3. **Nodes 04 and 05 have only seen stand-in structures**, built by reordering 5O45's chains. Their
   designed-chain-is-first assumption is unverified against real RFD4 output.
4. **Open Babel still prints its plugin noise** even with libxrender1 present, and `CCD_MIRROR_PATH`
   / `PDB_MIRROR_PATH` warnings appear on every invocation. Harmless, but they bury real output;
   setting those two variables in the image would quieten it.

## The monomer scaffold is synthesised, not vendored

`examples/monomer.cif` in the model repository is 3 KB of annotation with no real atoms: a single
sentinel `HETATM` at NaN coordinates, and the segment length written three ways — a component named
`100-100`, and `annotation_expsegmin_1_residue` / `annotation_expsegmax_1_residue`. Node 02 writes
the same file at the requested `binder_length`, so monomer works at any length instead of only the
100 upstream ships.

The format is undocumented, so the implementation is inferred and was therefore proven by running
it: the spec validates through RFD4's own parser with `mon_id 60-90`, and generation returns a
70-residue single chain in 36.9s. If an upstream version changes those field names this breaks, and
it breaks loudly at node 02's validation rather than silently.

The scaffold replaces `target.cif` rather than travelling as a file of its own. A file present for
one task in four would need conditional entries in two `job.toml` contracts, and Silva's behaviour
on a declared-but-absent file is unverified; writing it as the target keeps one contract across all
four tasks and the spec's `input:` constant.

## Corrected after reading the package

Before building, the package itself was read: PyPI metadata, the wheel's `entry_points.txt` and
`METADATA`, the Hydra configs, `rfproteina/utils/hub.py`, `rfproteina/utils/mask_spec.py`, and
upstream's README plus the condition-spec guide and cookbook. Six things the plan asserted from
documentation were wrong.

| Assumed | Actual |
|---|---|
| `--config=spec.yaml` | `inputs=spec.yaml` (Hydra, `main_hydra` on `configs/inference.yaml`) |
| `inference.num_replicates`, `inference.seed` | top-level `num_replicates`, `seed` |
| import name is `rfd4_proteina` *or* `rfproteina` | both exist, plus `script_utils`; a false dilemma |
| `write_condition_spec_cifs.py` script | no such script; `resolve_condition_spec_inputs()`, a function |
| spec is a bare `{input, contig, conditions}` mapping | keyed by an example id, conditions keyed by condition *value* with `{select: ...}` entries |
| HF repo is ~1.74 GB | 82.57 GB: AlphaFold2 (51.6) and OpenFold3 (21.2) refolding runtimes, native `.ckpt` copies, and the 1.74 GB bundle. The loader itself asks only for `config.json` + `model.safetensors` |

Two assumptions held: `res_id` does resolve against mmCIF `label_seq_id` ("before contig renumbers
the output"), so the whole node 01 numbering design stands; and designs do land in
`<out_dir>/inference_samples/`.

Worse than any of the above, a `binder` run set **no conditions at all** unless hotspots were
given, so the target was never held in place. Every task now applies `C_CRDGRP` over backbone atoms
plus `C_SEQ`, which is upstream's own binder recipe.

Build-time facts now established: `rfd4-proteina==0.0.3` is the newest release (requires-python
`>=3.12,<3.13`); the `[all-cuda12]` extra resolves, pulling `tmol[cu12]` and the cuEquivariance
kernels; torch lands as 2.14.1+cu126 and reports `cuda.is_available() True` on driver 580.178.04;
the image is 10.1 GB with 1.7 GB of weights at `/models`. And `libxrender1`/`libxext6`/`libsm6` are
mandatory, not cosmetic: without them Open Babel's plugin load fails, `GetSupportedInputFormat()`
returns malformed data, and `rfproteina.inference_utils.condition_spec` cannot be imported at all.

## Verification plan

Per `SILVA_MIGRATION_GUIDE.md` §6 — all must pass before the workflow is considered done.

Build-time checks (these also resolve open items 1–3):

```bash
docker build -t rfd4-proteina:2026_11_04 workflows/workflow-037/
docker run --rm rfd4-proteina:2026_11_04 rfd4-generate --help          # CLI present
docker run --rm rfd4-proteina:2026_11_04 python -c "import rfproteina" # resolve import name
docker run --rm --gpus all rfd4-proteina:2026_11_04 \
  python -c "import torch; print(torch.cuda.is_available(), torch.version.cuda)"
```

Per-node loop (`AI_WORKFLOW_GUIDE.md` §5), node 03 shown:

```bash
mkdir -p /tmp/wf037/inputs /tmp/wf037/outputs
cp workflows/workflow-037/input_files/5o45_cropped.cif.zst /tmp/wf037/inputs/
docker run --rm --gpus all \
  -v /tmp/wf037/inputs:/workspace/inputs \
  -v /tmp/wf037/outputs:/workspace/outputs \
  -v $(pwd)/workflows/workflow-037/03-generate:/workspace \
  -w /workspace \
  -e PARAM_NUM_REPLICATES=2 -e PARAM_SEED=42 -e PARAM_DISABLE_COMPILE=true \
  rfd4-proteina:2026_11_04 bash run.sh
```

Node 03 first runs with `--dry-run` (prints resolved argv, exits 0, no GPU) — the testability
pattern workflow-026's `02_predict` uses.

Chain + end-to-end: `cp /tmp/wf037-nodeN/outputs/* /tmp/wf037-nodeN+1/inputs/` through all five
nodes, confirm `report.html` renders with a populated Mol* viewer, then run under Silva with
`export SILVA_WORKFLOW_HOME=…` and capture the green-checkmark screenshot the §6 checklist requires.

**Correctness gate:** the ground truth is the documented PD-L1 quickstart —
`contig: "50-120, /0, A17-131"` against `5o45_cropped.cif.zst` must produce `num_replicates` designs
whose binder chain length falls in [50, 120] and which contact the target. Anything else means the
spec builder is mis-rendering.

## Biosafety

RFD4's HF Spaces screen every input and design with SecureDNA; the model card is explicit that the
**CLI/Python path does not**. This workflow uses the CLI path, so designs come out unscreened.
**Decided: no screening node.** Instead, a prominent notice in the README and in `report.html`
states that designs are unscreened and require biosecurity review before synthesis.

## Process note

`SILVA_MIGRATION_GUIDE.md` §8 says to open a `[Workflow]` issue and wait for approval from Rongjun
before starting implementation. Rongjun left on 2026-10-19 and the repo now has a sole maintainer,
so that gate no longer applies. Nothing in this workflow is pushed; the drafts below are here to
file whenever desired.

### Draft issue body

> **[Workflow] RFD4-Proteina de novo protein design**
>
> Add a multi-task generative design workflow built on RFD4-Proteina (RosettaCommons + NVIDIA
> GenAIR, BSD-3), which supersedes RFdiffusion 1/2/3. Five nodes: acquire target → build and
> validate condition spec (CPU) → generate (GPU) → analyze → report. Supports binder,
> motif_scaffold, monomer and binder_redesign tasks from an uploaded structure or an RCSB
> accession.
>
> Generation only — no refolding oracle; scoring is deferred to a follow-up workflow chaining into
> workflow-026. CUDA 12 (`[all-cuda12]`, `--torch-backend cu126`) to match the existing GPU estate.
> Designs are **not** SecureDNA-screened, since the RFD4 CLI path does no screening; the README and
> report say so prominently.

### Draft PR description

> Adds `workflows/workflow-037/` — RFD4-Proteina de novo design, 5 nodes on one shared image
> (`rfd4-proteina:2026_11_04`, `python:3.12-slim` base, weights baked at a pinned HF revision).
>
> - `01-acquire-target` — upload or RCSB fetch, normalise to mmCIF, inventory both residue numberings
> - `02-build-spec` — render the task condition spec and validate it to an annotated CIF on CPU
> - `03-generate` — `rfd4-generate` with `~metrics`, offline pinned weights, `--dry-run` support
> - `04-analyze` — geometric per-design metrics table
> - `05-report` — self-contained `report.html` with Mol* viewer
>
> Registry: root `README.md` entry #37, `workflows/CHANGELOG.md`, `apps/applications.json`.
>
> Designs are unscreened — the RFD4 CLI/Python path does not run SecureDNA. Flagged in README and
> report; no screening node by design.
