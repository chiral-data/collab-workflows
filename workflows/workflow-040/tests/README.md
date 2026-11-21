# workflow-040 test suite — plan and merge gate

Status: plan. Nothing here runs yet, because no node exists yet.

This document and `fixtures.json` beside it are the **acceptance criteria** for
workflow-040. Passing them is what makes it mergeable; `validation_pipeline.md` §8 is
what they are derived from.

## Why this suite is different from 037's

037's suite says, in writing, that it tests the plumbing and not the science:

> Passing means the pipeline does not silently corrupt, drop, mislabel or misnumber its
> input, and that its geometry flags fire where they should. It does not mean a design is
> good. (`5a2d9b3f`)

That was the honest claim for a generator. **It cannot be the standard here, because 040
is the scorer.** A 040 that passes every plumbing check and scores barnase–barstar badly
is worse than no 040 at all: it would reject good designs while looking like it worked.
029 shipped exactly that failure twice — once when a chain-order bug meant every design
scored iPTM ~0.2 (`d5ba7487`), once when a receptor trim moved the PAE floor out of reach
and two runs scored 2/16 then 0/16 (`8e4772a7`). Neither moved an exit code. Both looked
like "the designs are bad."

So the suite has two tiers, and **both** gate the merge.

## The two tiers

### Tier 1 — CI. No GPU, deterministic, runs on every change.

Plumbing, rejections, and scoring against synthetic predictions whose correct answer is
known by construction. This is the tier that catches the §8 defect classes.

It reaches nodes 04 and 05 without a GPU the same way 037 reached its nodes 04–06: those
nodes read structured output, not the predictor, so a hand-built prediction exercises them
exactly as a real one does, and deterministically, which a sampled prediction never could.

### Tier 2 — GPU. Run by hand on the GPU host, numbers recorded in the PR body.

The controls that can falsify the scoring itself. This mirrors how every 037 fix was
verified — each commit message carries its measured numbers, e.g. "Verified on the GPU
node: the clash is flagged at 0.49 Å and sorts last despite its 570 contacts"
(`304006fa`). The suite cannot run these in CI and does not pretend to.

## Merge gate

### Gate A — Tier 1, required for every phase

1. **Zero `FAIL`.** As in 037, known gaps carry `"status": "open"`, report as `GAP`, and
   do not fail the run.
2. **Every `GAP` is listed in "Known gaps" below**, with what it is and why it is not
   fixed. An unlisted gap is a fail.
3. **Every check demonstrated non-vacuous.** For each check, suppress the behaviour it
   asserts and confirm the check fails. Record this in the PR. This is not ceremony: 037
   shipped four checks that passed for the wrong reason, and a check that passes
   vacuously is worse than no check because it is also evidence (§8.4).
4. **Docs updated in the same commit** as the behaviour they describe (§8.15).

### Gate B — Tier 2, required per phase

**Phase 1 (self-consistency):** controls `C1`–`C4` and `C6`.

- `C1`, `C2` — native binders on their own backbones **must pass** the Phase 1
  thresholds.
- `C3`, `C4` — scrambled and poly-alanine sequences on those same backbones **must
  fail** them.
- The **separation margin** between the worst positive and the best negative is stated
  in the PR. If they do not separate, 040 is wrong and does not merge. This is the one
  check that can falsify the scoring rather than the plumbing.
- `C6` — a real prediction killed at a short timeout, with `nvidia-smi` showing the
  device freed afterwards, as `755d4af2` verified for 037.
- Every number recorded in the PR body, including both chains' pLDDT.

**Phase 2 (interface scoring):** adds `C7`–`C11`.

- `C7` is the three-way MSA diagnostic that settles §4 by measurement.
- Phase 2 does not merge if `C7` shows the target-only alignment failing to recover the
  full-MSA numbers — that is the §11.1 gate, and the answer is to re-make the §4 decision,
  not to patch around it.

**Phase 3 (calibration):** thresholds re-derived from `C1`–`C11` and written into the
README beside each number, with what was measured. Until then `pae_interaction` ships as
a provisional band, not a gate (§7).

## Harness

Four pieces. One already existed and was extended in place rather than forked.

| Piece | Status | Role |
|---|---|---|
| `make_designs.py` | **037's, extended in place** | builds the design CIFs that are 040's *input*. See below. |
| `make_predictions.py` | **to write** | builds the synthetic *predictions* — a `rank_001` PDB plus a scores JSON — with a known correct answer. |
| `probe_design.py` | **to write** | measures a design CIF independently: chains, lengths, sequences, residue ids, gaps. Must not share a parser with any node, per 037's `probe_entry.py`. |
| `check_fixtures.py` | **to write** | the `CHECKS` table. Each entry returns `(status, detail)`, and the detail string carries the numbers so a failure reads as evidence rather than as an assertion. |

### `make_designs.py` — extended, not forked

It already built six shapes that vary **the design**: `compact`, `extended`, `helix`,
`broken`, `clash`, `swap`, at any length, with exact 3.8 Å Cα spacing. 040 needs shapes
that vary **the target** and **the residue names** instead, because that is what its
ingest node has to get right, so six were added to the same file:

| Shape | What it is | Which case needs it |
|---|---|---|
| `monomer` | the design chain alone, no target | `monomer-single-chain`, `monomer-on-complex-path` |
| `multichain` | one design, two target chains **A and C** — skipping the design's own B, so nothing passes by assuming the targets are the letters after it | `multichain-target` |
| `gapped` | 4 residues absent from the target's numbering, second segment moved 12 Å so the junction is not a peptide bond either | `gapped-target` |
| `renumbered` | target numbered from 3, where barnase starts in 1BRS | `target-numbering-from-3` |
| `nonstandard` | 9KK, CCS, SAR among ordinary residues: some X, still foldable | `nonstandard-residues`, the count |
| `allnonstandard` | every design residue modified: the sequence is entirely X | `nonstandard-residues`, the skip |

A `--manifest` flag writes `generation_manifest.json` carrying the target chain ids the
generator **actually used**. 037's own harness writes that file inline with
`target_chains: ["A"]`, which is right for the six original shapes and wrong for
`multichain` and `monomer`; the flag is the general form, and it means the harness stops
restating a fact the generator decides.

Nothing the original six shapes produce changed. That was verified rather than asserted:
all six were generated from the pre-change file and the post-change file at lengths 10, 30
and 50, three designs each, and diffed — **byte-identical, 18 of 18**.

**`design_metrics.json` is not synthesised.** It is `designed_chain_source`'s only source,
and inventing a file to carry it would be inventing the answer. The harness runs
**workflow-037's node 04** over the generated structures and uses what comes out, so
`resolved-chain-good-path` asserts a value the real producer really produced. The negative
cases then mutate that file: `named-chain-absent` rewrites `designed_chain` to a chain the
CIF does not contain, and `no-metrics-file` withholds it.

### `make_predictions.py` — the part that carries the weight

It emits, per case, a prediction directory and a `truth.json`. Three constructions, chosen
so that the answer does not depend on the generator being right:

| Shape | Construction | True `sc_rmsd_A` | Why this one |
|---|---|---|---|
| `identity` | the design backbone, copied | **0, by construction** | needs no Kabsch to know the answer |
| `rotated` | the design backbone under a rigid rotation and a large translation | **0, by construction** | needs no Kabsch to know the answer, and it is the keystone test |
| `displaced:<d>` | alternate Cα displaced ±d perpendicular to the chain axis: centroid preserved, no rotation improves the fit | **d** | calibrates magnitude |

`identity` and `rotated` are the important two, and they are important because their
answer is independent of any superposition code. **A node 04 that forgets to superimpose
returns ~0 for `identity` and something enormous for `rotated`** — and that bug rejects
every good design while looking like a batch of bad designs, which is precisely how 029's
chain-order defect presented.

For `displaced:<d>`, the generator computes the truth with its **own** Kabsch and refuses
to emit the case if the measured value differs from `d` by more than 1e-6. That keeps the
case honest without pretending the number is free.

It also controls the scores JSON directly, so the per-residue `plddt` array, the `pae`
matrix blocks, and the presence or absence of `iptm` are all inputs to a case rather than
incidental.

### The scores-file shape is measured, not guessed

`workflow-026/sample_outputs/confidence_summary.json` holds a real ColabFold **1.6.1**
monomer run, 164 residues, five models. Measured from it:

- `plddt_per_residue` has length 164, exactly `chain_lengths[0]` — the array is
  per-residue and carries no chain boundaries.
- `pae` is 164×164, square and matching.
- `ptm` is present on all five models, 0.88–0.90.
- **`iptm` is absent.** So is `interface_pae_mean`. A monomer prediction has no interface
  quantity at all.

That last line is why §8.7 is a live hazard for Phase 1 and not a hypothetical:
029 reads `scores.get('iptm', 0.0)`, which on Phase 1's monomer-only path would report
**0.0 — a measured-looking bad interface — for every design in every run.** 026 does it
correctly, in one line: `iptm = float(data["iptm"]) if "iptm" in data else None`
(`026/03_analyze_confidence/script.py:86`). Copy 026, not 029.

Both of these were open when this plan was written, and both are now measured — see
[reference/](reference/) for the files and their provenance:

1. The **complex** scores file: captured. `plddt` is a flat `list[63]` spanning both
   chains, `pae` is `63×63`, `iptm` is present at 0.2. The predicted PDB put the target in
   chain A and the binder in chain B, so `TARGET:BINDER` makes the binder the **tail** of
   the array.
2. The **output filename prefix**: confirmed. A record named `mono_design_0` produced
   `mono_design_0_scores_rank_001_alphafold2_ptm_model_1_seed_000.json` and
   `mono_design_0_unrelaxed_rank_001_...pdb`. So the prefix is the record id exactly, and
   node 02's per-design filename flows straight through — which is what lets node 04 map
   output back to design by prefix rather than by listing order. Node 02 must still make
   ids filesystem-safe (`54bf21b6`).

One trap the capture exposed: the raw scores file and
`026/sample_outputs/confidence_summary.json` are **not the same shape**. The latter is
026's post-processed summary, with `plddt_per_residue`, `plddt_mean` and
`interface_pae_mean`. A parser written against those names fails on the real file. This is
`54dc7def` in miniature — it read the real package and found "several of the documented
guesses this workflow was built on were wrong." A fixture earns its place by measurement,
never by reputation.

## The experiments

Full detail, with params and expectations, is in `fixtures.json`. This is the map.

### Section `ingest` — node 01 and node 02, real and synthetic design CIFs

| id | What it catches | Pins |
|---|---|---|
| `resolved-chain-good-path` | `designed_chain_source == manifest-target-chains` on an ordinary input — that the **good path was taken**, not merely available | §8.2 |
| `reversed-chain-order` | the binder resolved by reading metrics, not by file order | #27 |
| `no-metrics-file` | absent `design_metrics.json` produces a stated fallback, not a silent guess | §8.2 |
| `monomer-single-chain` | `single-chain` accepted on the mono path and marked ineligible for the complex path | §5 |
| `named-chain-absent` | metrics naming a chain the CIF does not have is refused | §6 |
| `multichain-target` | all target chains recorded, in order, when the target is more than one chain | `53d1d460` |
| `gapped-target` | an internal numbering gap is measured and reported, so the fusion can be disclosed | §8.12 |
| `target-numbering-from-3` | a target whose numbering does not start at 1 maps by position, not by residue number | `ef74646e` |
| `nonstandard-residues` | `X` residues counted and reported; an all-`X` design skipped with a warning, not written as junk | `54bf21b6` |
| `mono-input-has-no-colon` | the `mono_` file contains no `:` | `d3d707ed`, inverted |
| `cplx-input-has-one-colon` | the `cplx_` file contains exactly one | same |
| `one-file-per-design` | N designs produce N single-record files and never a batch file | `818954e9` |
| `chain-order-recorded` | node 02 records the order node 04 will slice by | §8.1 |
| `rerun-over-populated-dir` | re-running is safe or refused, never half-overwritten | §8.11 |

### Section `predict` — node 03 under `--dry-run`

| id | What it catches | Pins |
|---|---|---|
| `argv-no-host-url` | `--host-url` never appears in argv | 026's unconditional default |
| `argv-model-type-mono` | `alphafold2_ptm` on the mono path | §6 |
| `argv-model-type-cplx` | `alphafold2_multimer_v3` on the cplx path | 029 |
| `argv-msa-mode-phase1` | `single_sequence` on the mono path | §4 |
| `timeout-kills-process-group` | a child that ignores SIGTERM is still killed, because generation spawns workers that outlive the direct child | `755d4af2` |
| `per-design-failure-continues` | one failed prediction is recorded and the batch continues | §6 |

### Section `scoring` — nodes 04 and 05, synthetic predictions

This is the core. Every case has an answer known before the node runs.

| id | What it catches | Pins |
|---|---|---|
| `rmsd-identity` | `sc_rmsd_A == 0` | — |
| **`rmsd-rotated-copy`** | **`sc_rmsd_A == 0` for a rotated, translated copy** — a missing superposition rejects every good design | keystone |
| `rmsd-known-displacement` | `sc_rmsd_A == d` for a construction whose truth is verified independently | — |
| `rmsd-length-mismatch` | a prediction of the wrong length is refused, not truncated and scored | §8.1 |
| `plddt-slice-binder-first` | the binder's pLDDT read from the head when the order is binder-first | §8.1 |
| **`plddt-slice-target-first`** | **the same, read from the tail when the order is target-first** — a node hardcoding 029's head-slice passes the previous case and fails this one | §8.1 |
| `plddt-length-mismatch` | `len(plddt) != n_target + n_binder` is refused | §8.1 |
| `missing-iptm-is-null` | a monomer's absent `iptm` becomes `null`, never 0.0 | §8.7 |
| `missing-pae-is-null` | an absent `pae` becomes `null`, never `inf` | §8.7 |
| `pae-both-blocks-averaged` | an asymmetric PAE matrix is averaged over **both** off-diagonal blocks | §7 |
| `plddt-target-recorded` | the target's pLDDT is recorded and never gates | `2359ffa4` |
| `merge-by-name-not-position` | a shuffled metrics file does not shift metrics onto the wrong design | `818954e9` |
| `no-score-file-is-failure` | a missing `rank_001` file is a recorded failure, not a zero | §8.7 |
| `cplx-returned-one-chain` | a complex prediction that came back single-chain fails rather than scoring | 026's bug |
| `failed-design-appears-in-report` | a design with no prediction is a row with a stated reason, not an omission | §8.5 |
| `metrics-are-columns` | every gating metric is a rendered column, read from the report's `DATA` payload and never from its prose | §8.4, `179003b9` |
| `source-shown-in-report` | `designed_chain_source` reaches the reader | §8.2 |
| `numeric-ordering-at-12` | `design_10` does not sit between `design_1` and `design_2` | #51 |
| `offline-report-keeps-provenance` | with the CDN unreachable, the provenance and biosafety block still renders | `cb418341` |
| `biosafety-tag-travels` | `passing_sequences.fasta` records carry the tag, with no `:` in it and no `;` banner | `7709809b` |
| `threshold-recorded-beside-value` | a reader can re-filter without re-running | §7 |

### Section `rejections` — every guard shown refusing, and one shown accepting

| id | What it catches |
|---|---|
| `threshold-nan` / `threshold-inf` | refused **by name**; neither is caught by a range check, because every comparison against `nan` is false (`55edd8fc`) |
| `threshold-non-numeric` | fails the same way every other bad parameter does (`975ab379`) |
| `num-models-non-integer` | a count read as a float: `float("2.5")` succeeds where `int("2.5")` raises, so the cast decides whether the typo is caught at all |
| `binder-too-short-for-rmsd` | fewer than 4 Cα makes the superposition exactly determined by its own inputs, so the RMSD is a *small* number without a meaning — it reads as an excellent design |
| `sequence-outside-alphabet` | a code 026 would refuse inside the predictor, refused here instead, before the GPU is allocated |
| `inventory-arrays-misaligned` | sequence and coordinate arrays of different lengths, which would shift every residue after the break |
| `threshold-negative` | refused |
| `monomer-on-complex-path` | refused, not scored against an absent target |
| `cplx-input-missing-colon` | refused at node 02, where it costs nothing, rather than at the predictor |
| `boundary-value-accepted` | a threshold exactly on its limit is **accepted** — so no guard can pass by refusing everything (`4b1a611a`, `b6b741e9`) |

Every parameter is refused at **node 01**. It is the only node guaranteed to run before
the GPU node, and a `nan` first noticed by node 05 has already cost the GPU its minutes.
Nodes 04 and 05 still parse their own thresholds with the same helper, so a node run
standalone is safe too — node 01 is the gate, not the only check.

Three of these entries asserted a message string before any node existed, and the built
node showed two of them were wrong about their own subject (§8.16, suspect the harness
first). They now carry a `measured` field saying so:

- `threshold-inf` expected the substring `inf`, which is in neither `1e400` nor "must be
  a finite number". The node now names the **parsed** value — "which parses as inf" —
  which is both better diagnostics and an assertion that can actually hold.
- `threshold-non-numeric` expected "whole number" for `min_plddt_binder`, which is a
  float; 037's helper says "a whole number" only for an `int` cast. The float path says
  "must be a number", and the `int` path got its own entry so both messages are pinned.

### Section `controls` — Tier 2, GPU

| id | Phase | Experiment | Must |
|---|---|---|---|
| `C1` | 1 | 4ZQK: PD-1's own sequence on its own backbone, mono path | score well |
| `C2` | 1 | 1BRS: barstar's own sequence on its own backbone, mono path | score well, and exercise barnase's label numbering starting at 3 |
| `C3` | 1 | the `C1` backbone with its sequence scrambled | score badly |
| `C4` | 1 | the `C1` backbone with a poly-alanine sequence | score badly |
| `C5` | 1 | real archived 037 designs | record the numbers — this is 040's first honest ranking |
| `C6` | 1 | a prediction killed at a short budget | the whole process tree dies and `nvidia-smi` shows the device freed |
| `C7` | 2 | **4ZQK complex, three ways**: `single_sequence`, target-only MSA, full MSA — reading `iptm`, `pae_interaction` and **both** chains' pLDDT each time | settle §4 by measurement; the §11.1 gate |
| `C8` | 2 | 1BRS as a native complex | good interface metrics |
| `C9` | 2 | the same complex, binder scrambled | bad interface metrics |
| `C10` | 2 | a real `cplx_` input with its `:` stripped | come back one chain and be refused |
| `C11` | 2 | a real target carrying an internal gap, folded whole | both chains' pLDDT recorded; the whole-chain fold not degraded as `8e4772a7`'s trim was |

`C7` is the experiment that earns its GPU time twice over: it is both the Phase 2 gate and
the calibration input for `pae_interaction`, and it is the reason `plddt_target` is a
first-class output rather than a diagnostic afterthought.

## Choosing a fixture

Same rule as 037, and it is the rule that matters most here:

> A fixture earns its place by **measurement**, never by reputation.

Every entry in `fixtures.json` carries an `evidence` string recording what was measured
and by what. Where 040 reuses an entry 037 already probed — 4ZQK, 1BRS, 12E8, 5O45 — the
evidence is **inherited from 037's `fixtures.json` rather than restated from memory**,
which is the same "don't re-derive an upstream fact" rule the workflow itself follows.

029's counter-example is why this is not pedantry: its README described its benchmark as
IL-7Rα and cited hotspots to a paper, when the file was the IL-2 receptor complex and the
hotspots were residue numbers from a different protein. They "passed validation because
they exist, not because they meant anything" (`11ddf044`). 040's controls set the
thresholds, so a mislabelled control miscalibrates the whole workflow.

## Scope

62 experiments, 58 distinct checks to implement.

| Section | Cases | Phase 1 | Phase 2 | Blocked |
|---|---|---|---|---|
| `ingest` | 14 | 11 | 3 | — |
| `predict` | 5 | 4 | 1 | — |
| `scoring` | 21 | 14 | 7 | — |
| `rejections` | 11 | 9 | 2 | — |
| `controls` (GPU) | 11 | 6 | 5 | 2 |

So Phase 1's gate is 38 CI cases plus 6 GPU controls. Nothing is blocked on the harness
any more: the six generator shapes four `ingest` cases were waiting on are built.

The count grew by four while nodes 01 and 02 were written, and every one came from
building rather than from planning:

- `num-models-non-integer`, when the built node showed `threshold-non-numeric` was
  asserting a message for the wrong cast.
- `binder-too-short-for-rmsd`, `sequence-outside-alphabet` and
  `inventory-arrays-misaligned`, which are three guards node 02 needed and nothing else in
  the suite could be seen refusing (§8.13). One of them, `sequence-outside-alphabet`,
  **cannot be reached from any generator shape at all**, because node 01 maps every
  unrecognised residue to `X` and `X` is in the alphabet; it mutates the inventory, and
  its `measured` field says so.

## Known gaps

Nothing is built, so nothing is reproduced-but-unfixed yet. When that happens it is
recorded here and marked `"status": "open"` in `fixtures.json`, so the suite reports it
without failing. Deleting the status field is what makes the suite enforce it, and a
regression then shows up as `FAIL` rather than as silence.

Carried forward as blocked at the start of the build. Each is a build task, not an
excuse — listed here so none of them becomes a silent omission.

**Needs a measurement before the case can exist:**

- ~~The complex scores-file shape is unmeasured.~~ **Done.** Both the monomer and the
  complex scores files are captured in [reference/](reference/), measured on the GPU host
  with `--network none`. `iptm` is absent from the monomer file and present at 0.2 in the
  complex one; neither carries a chain boundary. The working target-only-MSA a3m is
  captured there too, which closed §11.1.
- **`C11` needs a target with a measured internal gap.** 029 records chain D of 2B5I
  carrying a 38-residue gap whose ends are 33.9 Å apart (`11ddf044`), which is a starting
  point — but it must be re-probed rather than taken from a commit message, which is the
  rule in "Choosing a fixture" above.

**Needs inputs that are not in the repo:**

- **`C5` has no inputs.** 037's archived GPU designs are not committed. Either commit a
  few CIFs plus their `design_metrics.json`, or regenerate them on the GPU host.

**Needs a generator shape `make_designs.py` does not have:**

Four `ingest` cases carry a `needs` field and cannot run until it is built —
`multichain-target`, `gapped-target`, `target-numbering-from-3` and
`nonstandard-residues`. All four are extensions to 037's generator rather than new code:
it already writes mmCIF with exact Cα spacing, so each is a variation on its chain-writing
loop. Extend it in place, in 037's tree, so both workflows' suites get them.

## What this suite does not establish

Stating this plainly is part of the suite, for the same reason 037 states its own limits
(`5a2d9b3f`): reading the control set as proof of more than it shows would be the most
expensive mistake available here.

- **It does not establish that a passing design binds.** There is no experimental binding
  data anywhere in this repo. Everything here is in silico, and a design that passes every
  threshold is a better hypothesis, not a result.
- **The positive controls are native complexes, which are an easier problem than
  designs.** A native pair has co-evolutionary signal a designed one cannot have (§4). `C1`
  and `C2` establish that the scorer is not broken; they do not establish that it
  discriminates among designs.
- **The negative controls are scrambles, and a scramble is an easy negative.** The hard
  case is a design that looks fine and does not bind, and nothing in this repo can produce
  one without wet-lab data. So the suite can catch a scorer that is broken or merely
  generous; it cannot measure how well a working scorer ranks real candidates.
- **Five models, one seed, no variance estimate.** 029 recorded all five AF2 models
  agreeing on one native complex (iPTM 0.850–0.857). Nothing establishes that designs
  behave as tidily, and 037 found that its own seed reproducibility was empirical rather
  than guaranteed, since `torch_deterministic` is left false upstream (`11bb5149`).
- **Nothing here tests 037.** 040 consumes 037's output and `validation_pipeline.md` §10
  requires no changes to it. If 037 changes its contract, this suite notices only through
  the `ingest` section's fixtures.
