# workflow-040 — Binder Design Validation

**Scaffold only. No node is implemented; every one exits 1.**

This README is a skeleton. It gets filled in node by node, in the same commit as the node
it describes — a README that describes behaviour the code does not have is a bug in this
repo, not a follow-up (`CLAUDE.md`).

- Design spec and decisions: [validation_pipeline.md](validation_pipeline.md)
- Test suite and merge gate: [tests/README.md](tests/README.md),
  [tests/fixtures.json](tests/fixtures.json)

## What this workflow is for

Workflow-037 generates binder candidates and stops. Every metric it reports is geometric
by declaration — contact counts, radius of gyration, net charge, composition — and every
one of them is computable on a design that folds into nothing and binds nothing. On such
a design they all look ordinary. 037 says so itself, and its report deliberately ranks
nothing as a result.

This workflow refolds each designed sequence with a predictor that never saw the design,
and asks whether the sequence encodes the structure it was designed as. That is
*self-consistency*, and it is the first filter that can tell a candidate from a
plausible-looking string.

It reads 037's **structures**, not its exported FASTA. A sequence alone cannot be compared
against a backbone that stayed behind, which is why the handoff is the `*_generated.cif`
files plus `design_metrics.json`.

## When not to use this

- **It does not tell you a design binds.** A sequence can recover its own backbone and
  bind nothing. Interface scoring is Phase 2 and is not built; see
  [§12](validation_pipeline.md) for the phasing and
  [§4](validation_pipeline.md) for why the complex path needs a decision that Phase 1
  does not.
- **There is no experimental data anywhere in this repo.** Everything here is in silico.
  A design that passes every threshold is a better hypothesis, not a result.
- **The thresholds are uncalibrated.** They are literature conventions until the Phase 3
  controls run. `tests/README.md` says what the suite does and does not establish.
- **Designs are unscreened.** 037's designs carry no SecureDNA screening, and passing
  validation makes synthesis *more* likely, not less. The notice travels into every
  artifact here.

## Nodes

| Node | GPU | Reads | Writes |
|---|---|---|---|
| `01-ingest-designs` | no | `*_generated.cif`, `design_metrics.json` | `design_inventory.json`, `design_metrics.json` (passed through) |
| `02-build-prediction-inputs` | no | `design_inventory.json` | `mono_<design>.fasta` (one record each), `prediction_inputs.json` |
| `03-predict` | **yes** | `mono_*.fasta`, `prediction_inputs.json` | `<design>_scores_rank_*.json`, `<design>_*rank_*.pdb`, `prediction_manifest.json` |
| `04-score` | no | inventory, predictions, manifests | `design_scores.json`, `design_scores.csv` — one row per design, failures included |
| `05-filter-and-report` | no | `design_scores.json` | `passing_sequences.fasta`, `validation_metrics.csv`, `report.html` (static; needs no JavaScript to be complete) |

```
[037] *_generated.cif + design_metrics.json
        │
  [01-ingest-designs] ──► design_inventory.json  (binder chain, sequence,
        │                                          residue ids and Ca coords IN ORDER)
  [02-build-prediction-inputs] ──► mono_<design>.fasta, one record each
        │
  [03-predict]  (GPU) ──► <design>_scores_rank_*.json, <design>_*rank_*.pdb,
        │                 prediction_manifest.json      (flat, prefixed per design)
        │
  [04-score] ──► design_scores.json   (sc_rmsd_A, plddt_binder)
        │
  [05-filter-and-report] ──► passing_sequences.fasta, report.html
```

`01b-align-target` is specced but not scaffolded: it belongs to Phase 2. The a3m question
that gated it is [closed by measurement](tests/reference/) — `colabfold_batch` does accept
an alignment for the target with the binder present only as its own query — so what is
left is building it, not deciding whether it can exist. The letter suffix is so that adding
it later does not renumber nodes 02–05, which would invalidate every saved `params.json`
and every cross-reference.

Node 01 also **refuses every bad parameter**, including the ones it never reads. It is the
only node guaranteed to run before the GPU node, and a `nan` threshold first noticed by
node 05 has already cost the GPU its minutes. Nodes 04 and 05 parse their own thresholds
with the same helper, so a node run on its own is still safe.

## Image

`colabfold:2026_10_08`, built from
[workflow-026's Dockerfile](../workflow-026/Dockerfile). This workflow defines no image of
its own, for three reasons:

- It is the only ColabFold image in this repo, and workflow-029's node 04 already borrows
  it rather than forking one.
- Its base is `ghcr.io/sokrypton/colabfold:1.6.1-cuda12` — the same ColabFold version that
  produced `workflow-026/sample_outputs/confidence_summary.json`, which is the measured
  reference for the scores-file shape node 04 parses.
- It bakes the AlphaFold2 weights at build time (`python -m colabfold.download`), so
  `single_sequence` mode has no reason to reach the network.

It is a **local tag, not in the `ghcr.io/chiral-data` registry** — there is no
`apps/c/colabfold_*` — so it exists only where someone has built it from 026's Dockerfile.
And the `chown -R root:root /cache` in that Dockerfile has to stay in the same `RUN` as the
download: the upstream weight archives carry owner IDs above 65535, which hosts with a
limited UID/GID range cannot apply, and a later `RUN` would copy the weights into a new
layer.

Measured in the image on 2026-11-21: ColabFold 1.6.1, numpy 2.5.3, biopython 1.85, scipy
1.17.1, jax 0.5.3 resolving to `CudaDevice(id=0)`, and 5.3 GB of baked weights in `/cache`.
Biopython 1.85 is the version that warns rather than raises on a FASTA comment line, which
is the hazard behind node 05's "no `;` banner" rule.

### Working on the GPU host

`amd-c12-1` already carries `colabfold:2026_10_08` and `rfd4-proteina:2026_11_04`, so
neither needs building. Two things about that host cost an hour to find and are worth
writing down:

- **Docker is rootless**, and its socket is `/run/user/1009/docker.sock`.
- **`DOCKER_HOST` exported in `~/.bashrc` points at `~/.docker/run/docker.sock`, which
  does not exist.** Any non-login shell — which is what `ssh host 'cmd'` and any script
  run with `bash -s` gets — therefore fails with "Cannot connect to the Docker daemon",
  and `docker images` reports nothing at all rather than reporting an error. That is how
  this workflow briefly concluded the host had no images on it.

So scripts here set `DOCKER_HOST=unix:///run/user/1009/docker.sock` explicitly rather than
relying on the profile.

## Parameters

Declared once, at workflow level. See [`.chiral/workflow.toml`](.chiral/workflow.toml) for
the hints, which carry the reasoning.

| Parameter | Type | Default | Notes |
|---|---|---|---|
| `max_sc_rmsd_a` | float | `2.0` | self-consistency cutoff, Å. Scale-dependent; a default to calibrate |
| `min_plddt_binder` | float | `80.0` | fold confidence, not binding |
| `num_models` | integer | `5` | scores read from `rank_001` |
| `num_recycle` | integer | `3` | right for a monomer; 026 records 6–12 for multimers, so Phase 2 should raise it |
| `msa_mode` | string | `single_sequence` | correct for a lone de novo binder, wrong for a complex — [§4](validation_pipeline.md) |
| `predict_timeout_s` | integer | `1800` | unmeasured for this workflow — [§11.5](validation_pipeline.md) |
| `top_n` | integer | `10` | size of the shortlist the report highlights; every design is still in the table |
| `score_interface` | boolean | `false` | Phase 2, not built. Setting it true is **refused**, not ignored |

## Status

| | |
|---|---|
| Scaffold | done |
| `01-ingest-designs` | **implemented** |
| `02-build-prediction-inputs` | **implemented** (mono path; `cplx_` is Phase 2) |
| `03-predict` | **implemented**, with `--dry-run` for CI |
| `04-score` | **implemented** (Phase 1 metrics; `iptm` and `pae_interaction` read when present) |
| `05-filter-and-report` | **implemented**; report is static HTML, no structure viewer in Phase 1 |
| Test generator | `workflow-037/tests/make_designs.py` extended in place: six new shapes and `--manifest`; `tests/make_predictions.py` written |
| Test harness | **built** — `run_fixtures.sh` stages and runs, `check_fixtures.py` judges and falsifies |
| Tier 1 suite | **47 pass, 0 FAIL, 5 gap** of 52 CI cases; 55/55 checks implemented, invoked and falsifiable. Run on the host, **not yet in the image**, so Gate A is met but not satisfied — see [tests/README.md](tests/README.md) |
| Tier 2 controls | not run. Needs the GPU host |
| Phase 2 | specced, §11.1 closed by measurement. The 5 gaps are all the complex path |
| Phase 3 | thresholds uncalibrated |

```
bash tests/run_fixtures.sh          # stage, run, judge, falsify  (~1 min)
```

It prints a loud banner when it is not running in `colabfold:2026_10_08`, and says in its
verdict that such a run is a dry run. A suite that quietly certified something it never
ran in the right place would be worse than no suite.

Prerequisites from [§11](validation_pipeline.md) and [tests/README.md](tests/README.md):

1. ~~Confirm `colabfold:2026_10_08` is present on the GPU host~~ — **done**, already built
   there, 14.7 GB.
2. ~~Confirm `colabfold_batch` does not reach the network in `single_sequence` mode~~ —
   **done**, a prediction completed under `docker run --network none` in 61 s.
3. ~~Capture a real monomer prediction, and the complex equivalent~~ — **done**, both in
   [tests/reference/](tests/reference/) with provenance. This also closed the §11.1 a3m
   gate, so Phase 2 proceeds as designed.
4. ~~Extend `workflow-037/tests/make_designs.py` with the shapes four `ingest` cases
   need~~ — **done**, in place, six shapes plus `--manifest`; the original six verified
   byte-identical at three lengths. See [tests/README.md](tests/README.md).

Nothing is outstanding on that list. What remains is a Tier 1 run inside the image, the
Tier 2 GPU controls, and Phase 2.
