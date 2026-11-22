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
  bind nothing. Interface scoring is built, off by default, and **uncalibrated**: `iptm`
  and `pae_interaction` are reported and do not gate, because no control has yet been run
  against either threshold and 115 alanines already score 95.81 pLDDT. Turning
  `score_interface` on buys you a number to look at, not a verdict. See
  [§12](validation_pipeline.md) for the phasing and [§4](validation_pipeline.md) for the
  MSA decision the complex path rests on.
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
| `num_recycle` | integer | `3` | right for a monomer; the complex path has its own `num_recycle_cplx`, default 6, after 026's 6–12 for multimers |
| `msa_mode` | string | `single_sequence` | correct for a lone de novo binder, wrong for a complex — [§4](validation_pipeline.md). The complex path uses `msa_mode_cplx` and **refuses** `single_sequence`, which discards a precomputed alignment while still re-saving it |
| `predict_timeout_s` | integer | `1800` | unmeasured for this workflow — [§11.5](validation_pipeline.md) |
| `top_n` | integer | `10` | size of the shortlist the report highlights; every design is still in the table |
| `score_interface` | boolean | `false` | Phase 2. Setting it true adds a second prediction per design, through the multimer model against the aligned target, and requires `01b-align-target` to have produced an alignment |

## Status

| | |
|---|---|
| Scaffold | done |
| `01-ingest-designs` | **implemented** |
| `01b-align-target` | **implemented** — Phase 2. Aligns the target once per run; the only node permitted egress, and the designed sequences are not inputs to it |
| `02-build-prediction-inputs` | **implemented**, both paths: `mono_<d>.fasta` always, plus `cplx_<d>.fasta` and `pairedmsa_<d>.a3m` when `score_interface` is on |
| `03-predict` | **implemented**, both paths, with `--dry-run` for CI. The complex path has its own `msa_mode` and **refuses `single_sequence`**, which silently discards a precomputed alignment |
| `04-score` | **implemented**, both paths. Self-consistency from the monomer, `iptm` / `pae_interaction` / `plddt_target` / `plddt_binder_in_complex` from the complex. Interface metrics are **reported and do not gate by default** |
| `05-filter-and-report` | **implemented**; report is static HTML, no structure viewer. Says per row whether the interface was scored and whether it gated |
| Test generator | `workflow-037/tests/make_designs.py` extended in place: six new shapes and `--manifest`; `tests/make_predictions.py` written |
| Test harness | **built** — `run_fixtures.sh` stages and runs, `check_fixtures.py` judges and falsifies; `run_controls.sh`, `report_controls.py` and `falsify_controls.sh` do the same for the GPU tier |
| Tier 1 suite | **54 pass, 0 FAIL, 0 gap** of 54 CI cases on the host; 57/57 checks implemented, invoked and falsifiable, 55 falsified by a corruption and 0 vacuous. Every case reaches the thing it is about, which is new as of 2026-11-21 |
| **Gate A** | **met in `colabfold:2026_10_08` for Phase 1** — 47 pass, 0 FAIL, 6 gap of 53, 56/56 checks, 51 falsified, 0 vacuous. That run tested the commit before the Phase 2 cases landed, so it does not speak for the complex path; the image run against this code is in flight. One gap in the image is by design: the container-pid claim the image refuses to pass (§8.4) — see [tests/README.md](tests/README.md) |
| Tier 2 controls, Phase 1 | **Gate B met**, measured twice on the GPU with `--network none`. The `C12` designed positive passes at 0.93 Å / 95.05; scrambled and poly-alanine fail; `sc_rmsd_A` separates them by 17.9 Å. Judge falsified 6/6. Found that the native positive controls were the wrong controls and that `plddt_binder` cannot gate alone — [tests/reference/tier2-phase1-controls.md](tests/reference/tier2-phase1-controls.md) |
| Tier 2 controls, `C6` | **met** — a real `colabfold_batch` holding 2357 MiB killed at a 35 s budget, nothing surviving in a container that outlived the node, device back to 1 MiB. Found a node 03 defect: the predictor's stdout is block-buffered, so a kill lost it and the timeout reason was empty |
| Tier 2 controls, Phase 2 | `C7`–`C11` **not run.** No GPU has yet folded a complex through this workflow. `C5` (real 037 designs) is the only Phase 1 control still outstanding |
| Phase 2 | **built, and validated in CI only.** `01b` plus the `cplx_` path through 02–05, with all five of the suite's complex-path cases now passing and falsifiable. Making them run found four defects, one of which would have degraded every complex to a monomer on a real Silva run. What is still missing is a GPU: no interface metric in this workflow has been produced from a real multimer prediction |
| Phase 3 | thresholds uncalibrated. `iptm` and `pae_interaction` have had no control run against them, which is why they do not gate |

```
bash tests/run_fixtures.sh                     # Tier 1: stage, run, judge, falsify (~1 min)
bash tests/run_controls.sh                     # Tier 2: GPU, offline, ~6 min, needs the image
bash tests/falsify_controls.sh <control-run>   # show the Gate B judge refusing
bash tests/run_c6.sh <control-run> <out> 35     # C6: kill a real prediction at a budget
bash tests/run_controls_phase2.sh              # Tier 2 Phase 2: C8/C9, the interface pair
bash tests/falsify_controls_phase2.sh <run>    # show the Phase 2 judge refusing
python3 tests/audit_ports.py                   # declared ports vs what the code touches
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

Nothing is outstanding on that list. What remains is a Tier 1 run inside the image, `C5`,
and Phase 2.

Two things the Tier 2 controls changed, both recorded in
[validation_pipeline.md](validation_pipeline.md) beside the decisions they revise:

- **`plddt_binder` may reject a design, never admit one.** A 115-residue poly-alanine
  scores 95.81, above the real design, and clears the 80.0 threshold. `sc_rmsd_A` is the
  metric doing the discriminating (§8.6).
- **A native complex is not a usable positive control for this workflow.** §4 predicts at
  `single_sequence`; natives need the alignment they are not being given, and score 49.33
  and 43.09 pLDDT. They are diagnostics now (§8.4, §9).
- **A killed prediction loses its stdout**, because `colabfold_batch` block-buffers into a
  pipe. Node 03 now falls back to colabfold's own log file for the timeout reason and says
  which source it used; the OOM scan had the same blind spot and is fixed with it (§6).
- **`--msa-mode single_sequence` silently discards a precomputed alignment.** ColabFold
  clamps `max_seq` to 1 in that mode, so the model sees one sequence — while the re-saved
  a3m, `config.json` and `coverage.png` all record that an alignment was used. One clause
  of §11.1's closure rested on `coverage.png`, and it was wrong: that plot is written from
  the *parsed* alignment, before the clamp. The complex path refuses `single_sequence` by
  name (§11.1, §8.4).
