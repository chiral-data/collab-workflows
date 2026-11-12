# workflow-037 — RFD4-Proteina De Novo Design (design doc)

Status: **scaffold committed, nodes not implemented.** This document is the design rationale and
the record of decisions; it also holds the draft `[Workflow]` issue body and PR description for
whenever the work is filed.

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

## Open items

1. **Import name is unresolved upstream.** The GitHub README shows `import rfd4_proteina as rfd4`;
   the HF model card shows `from rfproteina.proteina import Proteina`. The docs decline to say.
   Resolve empirically in the built image before any Python here imports the package.
2. **`[all-cuda12]` extra must actually resolve** against `tmol-ops-cu12`. If it does not, fall back
   to the no-extra install (generation still works, minus the fused cuEquivariance kernels).
3. **Driver must support cu126.** workflow-026 already runs `colabfold:1.6.1-cuda12` on this estate,
   so CUDA 12 is proven; cu126 specifically needs driver ≥ 525.
4. **Version pin is still a placeholder** in the Dockerfile — pin `rfd4-proteina==<version>` at
   first build.
5. **Sample target not vendored** — `input_files/5o45_cropped.cif.zst` still needs fetching.

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
