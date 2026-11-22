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

### Tier 2 — GPU. Run by hand on the GPU host, numbers recorded here.

The controls that can falsify the scoring itself. This mirrors how every 037 fix was
verified — each commit message carries its measured numbers, e.g. "Verified on the GPU
node: the clash is flagged at 0.49 Å and sorts last despite its 570 contacts"
(`304006fa`). The suite cannot run these in CI and does not pretend to.

**Phase 1 of this tier has run, twice, and it found more than Tier 1 did.** Five controls
through all five nodes on a real GPU with `--network none`: the scorer is right on a
designed binder (0.93 Å at 95 pLDDT, exported with its `UNSCREENED` tag), the native
positive controls were the wrong controls, and `plddt_binder` cannot gate alone because
115 alanines score 95.8. Numbers, provenance and reasoning in
[`reference/tier2-phase1-controls.md`](reference/tier2-phase1-controls.md).

## Merge gate

### Gate A — Tier 1, required for every phase

1. **Zero `FAIL`.** As in 037, known gaps carry `"status": "open"`, report as `GAP`, and
   do not fail the run.
2. **Every `GAP` is listed in "Known gaps" below**, with what it is and why it is not
   fixed. An unlisted gap is a fail.
3. **Every check demonstrated non-vacuous — mechanically, every run.** `check_fixtures.py`
   reads a tree on disk and nothing else, so corrupting that tree is a suppression of the
   behaviour. After judging, it copies each case's tree once per corruption, applies one
   edit a plausibly-wrong node would have produced, and re-runs every check that passed.
   A check that no corruption can break is reported `VACUOUS` and fails the gate. The
   handful of claims that do not live in an artifact at all — a process that no longer
   exists, three exit codes, the container's `--network none` — are named in
   `NOT_IN_THE_TREE` with what does falsify each, so "nothing broke it" cannot quietly
   mean two different things. This is not ceremony: 037 shipped four checks that passed
   for the wrong reason, and a check that passes vacuously is worse than no check because
   it is also evidence (§8.4).
4. **Docs updated in the same commit** as the behaviour they describe (§8.15).

Items 1–3 are enforced by the suite itself and printed at the end of every run. Item 1 and
3 also decide its exit code.

### Gate B — Tier 2, required per phase

**Phase 1 (self-consistency): MET, 2026-11-21.** Measured twice on the GPU host, offline.
Full numbers, provenance and findings in
[`reference/tier2-phase1-controls.md`](reference/tier2-phase1-controls.md).

| | role | len | `sc_rmsd_A` | `plddt_binder` | passes |
|---|---|---|---|---|---|
| `C12` designed helix | positive | 60 | **0.9272** | **95.05** | **yes** |
| `C3` scrambled | negative | 115 | 18.8308 | 44.94 | no |
| `C4` poly-alanine | negative | 115 | 49.0800 | **95.81** | no |
| `C1` 4ZQK PD-1 native | diagnostic | 115 | 16.5096 | 49.33 | no |
| `C2` 1BRS barnase native | diagnostic | 108 | 11.5737 | 43.09 | no |

- `C12` — a designed sequence on an idealised backbone **must pass**, and does: it
  recovers its own backbone to under 1 Å at 95 pLDDT, and is the one design exported,
  carrying the `UNSCREENED` tag.
- `C3`, `C4` — scrambled and poly-alanine sequences **must fail**, and do.
- The **separation margin** on `sc_rmsd_A` is **+17.90 Å**. If the positives and negatives
  do not separate, 040 is wrong and does not merge. This is the one check that can falsify
  the scoring rather than the plumbing.
- `plddt_binder` is **measured but does not gate**: `C4`, 115 alanines, scores 95.81 —
  above the real design, and past the 80.0 threshold. It may reject a design, never admit
  one. See finding 3.
- `C1`, `C2` — **downgraded from positives to diagnostics after being run.** They failed,
  and the control was wrong rather than the scorer: §4 predicts at `single_sequence`, which
  is the mode a native protein cannot be folded in. They are now the `single_sequence`
  baseline that Phase 2's MSA arms are measured against. See finding 2.
- The judge is shown refusing: `bash tests/falsify_controls.sh <run>` mutates one number
  at a time and requires a `NOT MET` each time. Six of six refused, and that pass found
  two defects in the judge itself.
- `C6` — a real prediction killed at a short timeout. **MET.** A 35 s budget against a
  `colabfold_batch` holding 2357 MiB: 85 predictor process lines sampled mid-prediction,
  nothing of it alive afterwards **inside a container that outlived the node**, the device
  back to 1 MiB, and all five designs recorded failed with the budget in the reason. It
  found a real defect in node 03 — the predictor's stdout is block-buffered, so a killed
  run lost it and the reason was empty; it now falls back to colabfold's own log file and
  names the source.
- Every number recorded here and in the reference file, including both chains' pLDDT.

**Phase 2 (interface scoring):** adds `C7`–`C11`. **The nodes are built; none of these
controls has been run.**

- `C7` is the three-way MSA diagnostic that settles §4 by measurement. One of its three
  arms is now in hand: `single_sequence` on a native chain gives 49.33 and 43.09 pLDDT,
  which is the number the other two arms have to beat.
- Phase 2 does not merge if `C7` shows the target-only alignment failing to recover the
  full-MSA numbers — that is the §11.1 gate, and the answer is to re-make the §4 decision,
  not to patch around it.

**Phase 3 (calibration):** thresholds re-derived from the controls and written into the
README beside each number, with what was measured. Until then `pae_interaction` ships as
a provisional band, not a gate (§7), and `plddt_binder` ships knowing a poly-alanine
clears it.

## Harness

Seventeen pieces: seven for Tier 1, ten for Tier 2. One already existed and was extended
in place rather than forked. Four more sit outside both gates and are listed at the end:
one static audit, two probes that answered a question once and are kept so the answer can
be re-taken, and one hand smoke test.

| Piece | Role |
|---|---|
| `make_designs.py` | 037's, **extended in place**. Builds the design CIFs that are 040's *input*. See below. |
| `make_predictions.py` | Builds the synthetic *predictions* — a `rank_001` PDB plus a scores JSON — with a known correct answer. |
| `make_target_msa.py` | Builds the target alignment node 01b reads on the complex path, with the query taken **out of the inventory node 01 just wrote** and the homologs built from it by point substitution. A hard-coded a3m would mismatch the query and 01b would correctly refuse it, which would test the refusal instead of the path. It refuses nothing itself: a monomer, or a batch disagreeing about its target, is 01b's refusal to make, and a harness that pre-empted it would leave the node's guard never shown refusing (§8.13). |
| `probe_design.py` | Measures a design CIF independently: chains, lengths, sequences, residue ids, gaps, and an md5 of the Cα coordinates in order. Shares no parser with any node, per 037's `probe_entry.py` (§8.4). |
| `run_fixtures.sh` | The entry point, at the name 037's suite uses. A wrapper. |
| `run_fixtures.py` | Stages and runs every CI case, and **decides nothing**. Leaves each step's inputs, outputs, log and exit code on disk. |
| `check_fixtures.py` | The `CHECKS` table and the falsification pass. Reads the tree; **runs nothing**. Each check returns `(status, detail)` and the detail carries the numbers, so a failure reads as evidence rather than as an assertion. |

Tier 2 has the same split -- one piece builds, one runs, one judges, one attacks the judge:

| Piece | Role |
|---|---|
| `make_controls.py` | Turns a deposited mmCIF into 040's input format: picks a chain by `auth_asym_id`, applies `native` / `scrambled` / `poly-alanine` to the **sequence only**, asserts the backbone did not move, and records what it did. |
| `make_c12_provenance.py` | Measures the provenance of `C12`, the one control that comes from 037's `make_designs.py` rather than the PDB, with the same reader and field names as the other four. |
| `run_controls.sh` | Fetches the entries, probes them, builds the controls, and runs 037 node 04 plus all five 040 nodes **in the image, on the GPU, with `--network none`**. Guards free space first and refuses below 8G. |
| `report_controls.py` | Judges Gate B Phase 1: per-control verdict, separation margin with its sign, provenance. Reads `design_scores.json`; recomputes nothing. |
| `falsify_controls.sh` | Mutates one number or flag in a copy of a control run and requires `report_controls.py` to refuse. A `MET` verdict means nothing without it (§8.13). |
| `run_c6.sh` | Control `C6`: starts a container that **outlives the node**, runs node 03 in it under a short budget against a real prediction, samples what is alive mid-run, and asks afterwards what survived. |
| `report_c6.py` | Judges `C6`'s five claims, and refuses a process listing it cannot show to be a listing. |

And the same four-way split again for Phase 2, because the interface metrics are a
different claim and share nothing with the self-consistency ones but the tree they are
read from:

| Piece | Role |
|---|---|
| `run_controls_phase2.sh` | `C8` and `C9` in **one run, one target, two binders**: barnase native and barnase shuffled on barstar. Identical coordinates, identical composition, identical numbering; the sequence is the only difference, which is what makes it a control rather than two experiments. Runs node 01b, both prediction inputs, the multimer model on the GPU, and nodes 04–05 — and shows node 03 refusing `single_sequence` on the complex path before spending GPU time on the real one. |
| `report_controls_phase2.py` | Judges Gate B Phase 2. `iptm` gates; `pae_interaction`, `plddt_target` and `plddt_binder_in_complex` are printed as diagnostics and do not, because no control has calibrated them. |
| `falsify_controls_phase2.sh` | Seven mutations, each of which must be refused: the positive collapsing, the negative scoring as well as the positive, the separation inverting, the gating metric going missing, a control row vanishing, a control appearing twice, and the `single_sequence` refusal never being exercised. |

**What `C8`/`C9` can and cannot establish** is worth being blunt about, because it is the
difference between a control and a calibration. The run is offline, so the target
alignment is **synthetic** — homologs built from the target's own sequence by point
substitution. §11.1 measured what that costs directly: on a synthetic a3m, `iptm` came out
*lower* with the alignment than without it. So the absolute numbers this run produces are
not threshold material and the report refuses to be read as though they were. What it
establishes is **separation**: whether 040 can tell a real binder from the same backbone
with its residues shuffled. That is the Gate B Phase 2 question, and §8.6 is why it has to
be asked — `iptm` is a number whatever went in.

`C7` is the control that needs a real target MSA, which needs one fetch from
`api.colabfold.com`. That is outward-facing, so it is not something the script does on its
own: `--allow-egress` turns it on, and until someone turns it on §4's decision stays half
measured. The binder is not an input to node 01b, so even then no designed sequence can be
sent — that is structural rather than a promise.

And three standing checks that are not part of either gate — one audit and two probes:

| Piece | Role |
|---|---|
| `audit_ports.py` | Checks every node's declared `inputs`/`outputs` against the `./inputs/<name>` and `./outputs/<name>` paths its code actually mentions. Written after node 02 wrote three files and declared one, which `silva validate` cannot see because it checks that declarations are well formed, not that they match the code. Run-time-built names are **printed, not skipped** — a tool that quietly ignored what it could not parse would be §8.3's case. An undeclared write has to appear in `DELIBERATELY_UNDECLARED` with its reason, so the exception list cannot grow by accident. Takes a second and stages nothing; it is a second opinion, not a replacement for the suite. |
| `probe_msa_clamp.sh` + `report_msa_clamp.py` | The measurement §11.1 rests on: the same a3m under both `msa_mode` values, offline, asking whether `single_sequence` discards it. It does. The probe's own verdict says that **which arm scores higher means nothing** — the homologs are synthetic and the finding is the difference, not its sign. |
| `smoke_phase2.sh` | Walks the complex path once on the host interpreter with no GPU, printing what each node wrote. It decides nothing and asserts nothing; it exists so the fixtures written against the path are written against something that has actually run. It is also the thing that did *not* catch node 02's undeclared outputs, because it copies `n2/outputs/*` wholesale — which is the difference between a smoke test and a suite. |

Running and judging are separate files on purpose. A runner that also judged could only
report pass or fail, and when something failed there was nothing left to look at — which
is how 037's first suite went. Here a verdict can be re-taken, with a changed check,
against a run already on disk:

```
bash tests/run_fixtures.sh                       # stage, run, judge
python3 tests/check_fixtures.py tests/fixture-runs/latest --verbose
python3 tests/check_fixtures.py tests/fixture-runs/latest --no-falsify
bash tests/run_fixtures.sh --only rmsd-identity  # or --section scoring
```

### How a case says what it needs

Every CI entry carries a `harness` field: one space-separated `key=value` grammar for all
four sections, not one per section, because a harness that has to know which section a
case came from is a harness that can stage the wrong thing.

| Token | Meaning |
|---|---|
| `shape=` | a `make_designs.py` shape (default `compact`) |
| `n=` / `len=` | how many designs, and the binder length |
| `metrics=` | `node04` (default) runs **workflow-037's own node 04** to produce the real `design_metrics.json`; `none` withholds it |
| `mutate=` | a named edit applied to a staged file after it is written and before the next node reads it |
| `also=` | a second, independent tree in the same case, for a claim that needs two shapes |
| `stage=` | how far to run: `ingest`, `rerun`, `predict-dry`, `predict-dry-cplx`, `predict-timeout`, `predict-partial`, `full` |
| `mp=` | everything after it goes to `make_predictions.py` verbatim, so it comes last |

The harness never restates what the entry already says: parameters come from `params`, and
whether a case must fail comes from `expect_error`.

`design_metrics.json` is **not synthesised.** It is 037 node 04's output, the whole
`designed_chain` contract rests on what that node actually writes, and a stand-in would
agree with the spec instead of with 037. Negative cases mutate the real file.

### Two stand-in predictors, on PATH

Node 03 resolves `colabfold_batch` from `PATH` and that is deliberately not overridable by
an environment variable, so two cases put a stand-in of that name in the case's own `bin/`
and prepend it. One installs `trap '' TERM` and leaves a worker in the same process group;
the other fails for exactly one design and writes plausible output for the rest.
`nvidia-smi` is stood in the same way, because those cases are about the subprocess
contract rather than about device detection. Both stand-ins print what they are, into the
log the checker reads.

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

| Shape | Construction | True `sc_rmsd_A` | Reported | Why this one |
|---|---|---|---|---|
| `identity` | the design backbone, copied | **0, by construction** | `0.0` | needs no Kabsch to know the answer |
| `rotated` | the design backbone under a rigid rotation and a large translation | **0, by construction** | `0.0005` | needs no Kabsch to know the answer, and it is the keystone test |
| `displaced:<d>` | a displacement field of RMS magnitude exactly *d*, projected orthogonal to every rigid motion, then rigidly moved | **d** | `1.9999` at *d* = 2 | calibrates magnitude |

#### Why `displaced` needs the projection, and how the generator proves it

A displacement with any component along a translation or a rotation is partly absorbed by
the superposition, so the true RMSD falls below *d* by an amount that depends on the
geometry — and then no fixture can state it.

With `P = Q + D`, both centred, the cross-covariance is `H = Qᵀ Q + Dᵀ Q`, and the optimal
rotation is the identity exactly when `H` is symmetric positive definite. `Qᵀ Q` always is.
The antisymmetric part of `Dᵀ Q` vanishes iff `Σᵢ (Qᵢ × Dᵢ) = 0`, which is precisely the
statement that `D` is orthogonal to every infinitesimal rotation. So the generator starts
from an alternating field, removes its mean, removes the least-squares rotation generator,
and only then scales the remainder to RMS *d*.

That is an argument, so it is also checked. `make_predictions.py` runs **its own Kabsch,
which node 04 never sees**, and refuses to write the case unless the construction measures
*d* to 1e-9 **and** the recovered rotation is the identity to 1e-9. If a geometry ever makes
`H` non-positive — *d* large against the structure's own extent — the generator fails loudly
rather than quietly handing node 04 an expected value that is wrong.

#### The reported numbers are not the constructed ones, and that is the format's floor

A PDB `ATOM` record carries three decimal places. Every coordinate is therefore rounded to
0.001 Å on the way out — about 0.0005 Å per axis — which puts a floor under the RMSD of
roughly 1e-4 Å. Measured on the 30-residue case:

| Shape | Constructed | As written to the PDB | Node 04 reports |
|---|---|---|---|
| `identity` | 0 | 2.06e-15 | `0.0` |
| `rotated` | 0 | 0.00050538 | `0.0005` |
| `displaced:2.0` | 2.0 | 1.99988314 | `1.9999` |

`identity` loses nothing because the lattice walk's coordinates are already whole
thousandths. `rotated` does, because rotating them lands between representable values.

Two consequences, both deliberate. `truth.json` records `sc_rmsd_A` **and**
`sc_rmsd_A_as_written`, so a check that wants better than 1e-3 compares against the file
rather than against the intention. And node 04 reports `sc_rmsd_A` to **four** decimal
places — one more than the data carries — because rounding to three would print the
keystone case as `0.001`, sitting exactly on its own tolerance.

None of which blunts the keystone: 0.0005 Å is five orders of magnitude below the
**259.98 Å** that the same two structures give when the superposition is skipped.

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
| `argv-model-type-mono` | `alphafold2_ptm` **and** `--msa-mode single_sequence` on the mono path | §6, §4 |
| `argv-model-type-cplx` | `alphafold2_multimer_v3` on the cplx path, still with no `--host-url` | 029 |
| `timeout-kills-process-group` | a child that ignores SIGTERM is still killed, because the predictor spawns workers that outlive the direct child | `755d4af2` |
| `per-design-failure-continues` | one failed prediction is recorded and the batch continues | §6 |

`--dry-run` is node 03's own flag, copied from `026/02_predict/script.py:37`. It resolves
every argv, writes the manifest and exits without starting the predictor or touching a
device — which is what makes the offline and model-type guarantees assertable on a machine
with no GPU. The flags *are* the guarantee, and they are decided in this node.

**What the node cannot assert about the network, and who does.** It can guarantee
`--host-url` is absent, and records `host_url_passed: false` in the manifest. It cannot
prove from inside its own process that no egress happened. That proof is the harness
running it under `docker run --network none`, which is how the offline claim was measured
in the first place — see [reference/](reference/). So `no_network_reached` is a check on
the *harness*, not on the node, and the node records the flags it would have needed.

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
| `missing-plddt-is-null` | an absent `plddt` becomes `null` and the design fails on "not measured", not on an invented number — and this is the one that **gates** | §8.7 |
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

| id | Phase | Experiment | Must | Status |
|---|---|---|---|---|
| `C12` | 1 | **the positive control**: 037's idealised 60-residue helix carrying its `ALEKIAELKA` designed sequence, mono path | score well | **0.9272 Å, 95.05 — passes** |
| `C3` | 1 | the `C1` backbone with its sequence scrambled | score badly | 18.8308 Å, 44.94 — fails |
| `C4` | 1 | the `C1` backbone with a poly-alanine sequence | score badly | 49.0800 Å, **95.81** — fails on RMSD only |
| `C1` | 1 | 4ZQK: PD-1's own sequence on its own backbone, mono path | *was* "score well"; now a `single_sequence` diagnostic | 16.5096 Å, 49.33 |
| `C2` | 1 | 1BRS: barnase's own sequence on its own backbone, mono path, exercising its label numbering starting at 3 | *was* "score well"; now a `single_sequence` diagnostic | 11.5737 Å, 43.09 |
| `C5` | 1 | real archived 037 designs | record the numbers — this is 040's first honest ranking | blocked: inputs not committed |
| `C6` | 1 | a prediction killed at a short budget | the whole process tree dies and `nvidia-smi` shows the device freed | **met** — 2357 MiB held, killed at a 35 s budget, nothing surviving in a container that outlived the node, device back to 1 MiB |
| `C7` | 2 | **4ZQK complex, three ways**: `single_sequence`, target-only MSA, full MSA — reading `iptm`, `pae_interaction` and **both** chains' pLDDT each time | settle §4 by measurement; the §11.1 gate | arm 1 of 3 in hand, from `C1`/`C2` |
| `C8` | 2 | 1BRS barnase on barstar, native | good interface metrics | **harness built, not yet run.** `run_controls_phase2.sh` |
| `C9` | 2 | the same backbone, binder sequence scrambled | bad interface metrics | **harness built, not yet run.** Same run as `C8`, same target |
| `C10` | 2 | a real `cplx_` input with its `:` stripped | come back one chain and be refused | not run. The CI case `cplx-input-missing-colon` now shows node 03 refusing this before the predictor starts, so what `C10` is left to establish is the half nobody can assert offline: what ColabFold actually returns for a colon-less complex |
| `C11` | 2 | a real target carrying an internal gap, folded whole | both chains' pLDDT recorded; the whole-chain fold not degraded as `8e4772a7`'s trim was | not built |

**`C12` was added after the first control run and it is the reason Phase 1 has a verdict
at all.** `C1`–`C4` between them cannot tell a mis-specified positive control apart from a
broken scorer: when the natives fail, both explanations fit. A designed positive in the
same run decides it. That it was missing from revision 1's plan is recorded as §8.4's
fifth instance, because "use a real experimentally validated complex" is the kind of
choice that passes review on its pedigree without anyone checking it against the
prediction mode two sections earlier.

`C7` is the experiment that earns its GPU time twice over: it is both the Phase 2 gate and
the calibration input for `pae_interaction`, and it is the reason `plddt_target` is a
first-class output rather than a diagnostic afterthought. `C1` and `C2` are now its
baseline arm rather than dead weight.

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

65 experiments, 56 distinct checks, all implemented, all invoked, none vacuous.

| Section | Cases | Phase 1 | Phase 2 | Blocked |
|---|---|---|---|---|
| `ingest` | 14 | 11 | 3 | — |
| `predict` | 6 | 5 | 1 | — |
| `scoring` | 22 | 15 | 7 | — |
| `rejections` | 11 | 9 | 2 | — |
| `controls` (GPU) | 12 | 7 | 5 | 1 |

So Phase 1's gate is 40 CI cases plus 7 GPU controls, and Phase 2 adds 13 CI cases and 5
controls. Nothing is blocked on the harness any more: the six generator shapes four
`ingest` cases were waiting on are built, and the one case still blocked is `C5`, which
needs real archived 037 designs that are not in this repo.

Read these off `fixtures.json` rather than trusting the table — it had drifted by three
while the suite grew, which is the ordinary fate of a count kept by hand:

```
python3 -c 'import json; d = json.load(open("tests/fixtures.json"));
print({k: len(v) for k, v in d.items() if isinstance(v, list)})'
```

### Where it stands

All 53 CI cases run, on both paths. Measured on the host interpreter on 2026-11-21,
after the Phase 2 cases landed:

| | |
|---|---|
| cases | **53 pass, 0 FAIL, 0 gap** of 53 |
| checks | 56 named, 56 implemented, 56 invoked |
| falsification | 54 broken by a corruption, 3 not in the tree and accounted for, **0 vacuous** |
| Gate A on the host | **met** |
| Gate A **in the image** | in flight against this code. The last image run — 47 pass, 0 FAIL, 6 gap of 53, 56/56, 51 falsified, 0 vacuous — tested the commit before Phase 2 and does not speak for it |
| Gate B Phase 1 | **met**, measured twice on the GPU, offline; `C6` met separately |
| Gate B Phase 1, falsified | 6 mutations, 6 refused — and that pass found 2 defects in the judge |
| Gate B Phase 2/3 | not met. `C7`–`C11` have never been run |

**The five complex-path gaps are closed.** Every case in the suite now reaches the thing
it is about, which is the first run where that has been true. The five are listed under
"Known gaps" below with what each was waiting for and the four defects that came out of
making them run — including one that would have degraded every complex prediction to a
monomer, silently, on a real Silva run.

**Gate A is only satisfied in `colabfold:2026_10_08`**, which is the only run that
counts; the host figures above are a dry run and the suite says so in its own verdict. The suite prints `RUNNING IN colabfold:2026_10_08 -- this is the gate.` at the top
of that run, and a dry run says so instead, in its header, in `run.json` and in its
verdict.

The image run differs from the host run in exactly one way, and it is the way it was
built to differ: **one gap rather than none.** It is `sigterm_ignoring_child_killed`,
which reports a gap in a container because a pid recorded inside one means nothing on the
host — the claim is established by the host run and the image run says so rather than
quietly passing (§8.4).

**Gate A in the image has not yet been re-measured against the complex path.** The
53-case image run recorded above is the code as of the commit before Phase 2's cases
landed, and it is not withdrawn: it really happened, against the code it names. It is
simply not a statement about the five cases that were gaps when it ran. Until the image
run against this code lands, Gate A for the workflow as it now stands is unproven, and
the host figures do not substitute for it.

The count grew by six — five while nodes 01–04 were written and one from running the GPU
controls — and every one came from building or measuring rather than from planning:

- `num-models-non-integer`, when the built node showed `threshold-non-numeric` was
  asserting a message for the wrong cast.
- `binder-too-short-for-rmsd`, `sequence-outside-alphabet` and
  `inventory-arrays-misaligned`, which are three guards node 02 needed and nothing else in
  the suite could be seen refusing (§8.13). One of them, `sequence-outside-alphabet`,
  **cannot be reached from any generator shape at all**, because node 01 maps every
  unrecognised residue to `X` and `X` is in the alphabet; it mutates the inventory, and
  its `measured` field says so.
- `missing-plddt-is-null`, the companion to `missing-iptm-is-null` and the one that
  matters more, because `plddt_binder` gates and `iptm` does not. A node defaulting it to
  0.0 fails the design and a node defaulting it to 100.0 passes it; both read as verdicts.
- `timeout-tail-falls-back-to-predictor-log`, from control `C6`, and the only case in the
  suite that exists because a **stand-in** was wrong rather than a node. The real
  predictor's stdout is block-buffered and a kill destroys it; `STANDIN_SIGTERM` flushes
  line by line, so `timeout_message_carries_output` was passing on output the real tool
  does not produce at that moment. There is now a second stand-in that reproduces the
  buffering, and the new case requires node 03's stdout capture to be genuinely zero
  bytes so it cannot drift into testing the ordinary path.

Where the built node's output differs from the constructed truth, a `measured` field says
by how much and why.

## Known gaps

### Closed: the complex path, measured

**These five were gaps in every run until 2026-11-21, and are now checks.** They are
recorded here rather than deleted, because what each one was waiting for, and what it
found when it finally ran, is the only reason to believe the check means anything.

| Case | What it was waiting for | What it found |
|---|---|---|
| `cplx-input-has-one-colon` | a `cplx_` file to count colons in | passes. One colon, halves of 20 and 10 residues |
| `gapped-target` → `fusion_disclosed` | a folded target, so a fusion could be disclosed | passes, and the check itself was wrong: see below |
| `argv-model-type-cplx` → `argv_includes` | a record asking for the multimer model | passes, and found the harness reading one argv per design where the complex path runs two |
| `monomer-on-complex-path` | the single-chain refusal to live somewhere reachable | passes. Node 01b refuses: a monomer declares no target chain, so there is nothing to align |
| `cplx-input-missing-colon` | an input that reaches node 03 missing its colon | passes, and moved node 03's guard from a per-design failure to a hard stop |

Four defects came out of making them run, and none of them was visible from reading the
code:

**Node 02 wrote three files and declared one.** `job.toml` listed `mono_*.fasta` and
`prediction_inputs.json` as its outputs while the node also wrote `cplx_*.fasta` and
`pairedmsa_*.a3m`. `silva validate` reported "6 job(s) ok", the smoke test passed, and
every complex would have degraded silently to monomer-only on a real run, because Silva
copies a node's *declared* outputs and nothing else. The suite caught it because its
harness reads those globs out of `job.toml` rather than carrying its own copy of the port
list — which is the one thing the harness docstring says that choice is for.

**`fusion_disclosed` passed on a key name.** It searched `prediction_inputs.json` for the
substring `fus`, and the payload carries a `fusion` key whether or not anything was fused.
It now counts the numbering gaps node 01 measured and requires each one's boundary
residues to appear in the disclosure. §8.4, again: a check that passes for the wrong
reason is worse than no check.

**The argv checks read one prediction per design.** `_argvs` was keyed on the design name,
and on the complex path node 03 runs each design twice — monomer through `alphafold2_ptm`,
complex through `alphafold2_multimer_v3`. One of the two silently replaced the other, so
whichever survived was the one every argv check then read, and `--model-type` would have
been compared against the wrong path's value depending only on manifest order. Keyed on
`design/mode` now, and `argv_includes` takes an optional `mode` so a flag whose correct
value differs by path can be asserted at all.

**A corrupt input was a per-design failure.** Node 03 recorded the missing colon against
that design and carried on, which is right for a predictor that fell over on one sequence
and wrong here: node 02 writes both files and reads them back before returning, so a
disagreement means something changed the tree after it ran, and nothing in node 03 can say
what else that touched. It is a hard stop now, and the per-design contract — 037's
measured lesson (`755d4af2`) — is untouched for actual prediction failures.

A fifth, smaller one: node 03's warning counted prediction records and called them
designs, so a design whose monomer was fine and whose complex was not read as
"1 of 2 design(s) have no usable prediction: design_0". Both counts are stated now.

### One claim each run cannot make

`timeout-kills-process-group` → `sigterm_ignoring_child_killed` passes on the host and
reports `GAP` **in the image**, which is why the two runs are not interchangeable.

The worker records its own pid and the harness looks for `/proc/<pid>` afterwards. Inside
a container that number belongs to the container's own pid namespace, so looking for it on
the host finds nothing whatever node 03 did — and the check would have passed for exactly
the reason §8.4 calls worse than no check. Probing from inside the container answers a
different question, because the container's exit ends the worker regardless. So the
harness records whether it could observe anything, and the check says so.

The claim is established by the host run. Everything else is established by the image run.
Both are needed and neither subsumes the other; the gate is the image run plus this one
line.

Seven other Phase 2 `scoring` cases pass for a different reason: node 04 slices from the
recorded `chain_order` rather than from an assumption, and reads `iptm` and
`pae_interaction` whenever they are present. Those ran green through the whole of Phase 1,
over predictions `make_predictions.py` built to carry interface metrics. They are
evidence about node 04's arithmetic and not about the complex path, which is why they
were never counted as Phase 2 being partly done.

### Carried forward as blocked

Each is a build task, not an excuse — listed here so none of them becomes a silent
omission.

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
  few CIFs plus their `design_metrics.json`, or regenerate them on the GPU host. This is
  now the most valuable missing case in the suite: `C12` establishes that the scorer works
  on a designed sequence, and `C5` is the only thing that would say how it ranks designs
  that were not built to be easy.

~~**Run by hand, not yet run:** `C6`.~~ **Done**, and it found a node 03 defect plus two
defects in its own judge. See [reference/tier2-phase1-controls.md](reference/tier2-phase1-controls.md).

**Needs a generator shape `make_designs.py` does not have:**

~~Four `ingest` cases carry a `needs` field.~~ **Done**, extended in 037's tree in place so
both workflows' suites get the shapes. No entry carries `needs` any more.

## What building the nodes found

Kept here because a suite that never caught anything is not evidence that the code is
right — it is evidence that nothing was asked of it.

**`rmsd-length-mismatch` caught a real defect in node 04.** The node compared the
prediction against node 02's *declaration*, which agreed, and never against the design it
is supposed to be a refolding of. So `zip()` silently truncated to the shorter array: a
32-residue prediction scored against a 30-residue design came out at **`sc_rmsd_A` = 0.0,
reported as PASSING**. Every other check in the node passed. The fix is a second, separate
assertion against the designed length, and the case's `measured` field records the
sequence of events.

That is the whole argument for the case, and it is worth noticing that the case was
written before node 04 existed and the defect was not one anybody predicted — the fixture
specified an outcome, and the outcome is what caught it.

**Two fixtures were wrong about their own subject** (`threshold-inf`,
`threshold-non-numeric`), and the built node is what showed it. Details above, in the
`rejections` section.

**One guard turned out to be unreachable from the generator** (`sequence-outside-alphabet`),
because node 01 maps every unrecognised residue to `X` and `X` is in the alphabet. The
case mutates the inventory instead, which is a less tidy fixture and a more honest one.

**The reported RMSD is not the constructed RMSD**, and that is the PDB format's three
decimal places rather than an inaccuracy. It changed how node 04 rounds. Details above,
under `make_predictions.py`.

## What running the Tier 1 suite found

**`no-score-file-is-failure` caught a second real defect in node 04.** §8.7 says a metric
that could not be computed is `null` with a recorded reason, because a null with no reason
cannot be told apart from a null nobody noticed. `score_one()` returns early in four cases,
and only one of them filled `null_reasons`. So a design whose prediction came back with no
`rank_001` scores file was reported with a blank `plddt_binder` and **nothing beside it** —
exactly the state the section forbids, arrived at through a different door. All four early
returns now go through one `bail()` that writes the row's cause into every metric, so a
branch added later cannot forget. Three of the four were unreachable from any fixture; the
one that was reachable is what found them.

**The case was itself wrong twice before it could find that.** It asked for a reason naming
`no_score_file`, a snake_case token no node in this repo writes, and the harness was
marking the design failed in the manifest — so node 04 was *reading a verdict* rather than
reaching one, and would have passed while noticing nothing. `make_predictions.py` now
records the prediction as `ok` and simply writes no scores file, which is the only version
of this case node 04 has to solve. §8.16, twice over, on one fixture.

**The first stand-in predictor was unfair, and the unfairness was worth writing down.** It
put its worker behind `setsid`, the worker survived the process-group kill, and the case
failed — correctly, because a process-group kill cannot reach a process that has left the
group. That is a real limit of 037's fix rather than a defect in it, and catching it would
mean walking `/proc` or a cgroup, neither of which 037 did or claimed. The stand-in now
reproduces what 037 actually measured; the limit is stated in `validation_pipeline.md` §6
instead of being hidden behind a stand-in nothing could pass.

**One check came out vacuous on the first falsification pass**, and one corruption was a
no-op. `mono_record_count` could not be broken, because the only corruption that changed a
record count merged several files and the case it runs on has one design; it now has a
corruption that appends a second record to a single file. And `blank-node-log` emptied the
`.log` files beside the tree while the checker reads stdout from `run.json`, so two
log-reading checks were being credited to deleting an artifact instead. Both are the
falsification pass doing its job on the suite rather than on the nodes.

**Four of the five harness strings were wrong when first written**, in ways that made cases
test less than they claimed: two lost the `make_predictions.py` options that created the
failure they were about, one withheld nothing where the case needed the fallback path, and
one ran two nodes where the check needed five. Every one of them showed up as a `FAIL` on
the first full run. The lesson is the one already in §8.16 and it keeps arriving: the
harness is the least-reviewed code in the repo and the most likely thing to be wrong.

## What running the Tier 2 controls found

Tier 1 found defects in nodes and in the harness. Tier 2 found defects in the *plan*,
which is the only tier that can.

**The pipeline is correct on the input class it is for.** `C12`, a designed sequence on an
idealised backbone, went through ingest, input construction, a real offline GPU
prediction, superposition, filtering and the report, and came out at **0.9272 Å and 95.05
pLDDT** — the one design exported, carrying its `UNSCREENED` tag with the thresholds
printed beside the values. Node 04's numbers matched the predictor's own log to the digit
(49.3→49.33, 43.1→43.09, 45→44.94, 95.8→95.81), which independently confirms the pLDDT
slicing and the parsing.

**The positive controls in revision 1's plan were the wrong controls.** Native complexes
predicted at `single_sequence` score 49.33 and 43.09 pLDDT, because §4 chose that mode for
designed binders and a native immunoglobulin domain is the opposite case. The plan asked
for the most defensible-sounding control available and nobody carried it back five
sections to check it against the prediction mode. Recorded as §8.4's fifth instance;
`C1`/`C2` are now diagnostics and the `single_sequence` baseline for §4's Phase 2 arms.

**`plddt_binder` cannot gate alone.** `C4`, 115 alanines, scores **95.81** — above the
real design's 95.05, and past the 80.0 threshold. Only `sc_rmsd_A` rejects it, at 49 Å.
This is §8.6 in a sharper form than the class was written for: the metric computable on
garbage does not look ordinary on garbage, it looks *excellent*, and it is a gating
metric. `plddt_binder` may reject a design, never admit one.

**The Gate B judge could lose its own verdict.** `report_controls.py` raised a `KeyError`
on a provenance field one control did not carry, after the results table and before the
verdict line — so a control run's `VERDICT.txt` is a file that reads like a complete
report and contains no verdict at all. Then the fix had the same defect inverted: a
*deleted* provenance record produced no row and passed, while a *partial* one was
correctly refused. Both found by `falsify_controls.sh`, neither visible by reading the
code, and both now §8.3 instances.

**4ZQK chain B has an 8-residue internal gap.** Probing before choosing, per §8.14, moved
the binder to chain A. Had it not, the positive control would have failed for a reason
unrelated to the score — and given that the natives failed anyway, the gap would have been
a second wrong explanation sitting underneath the first.

## What this suite does not establish

Stating this plainly is part of the suite, for the same reason 037 states its own limits
(`5a2d9b3f`): reading the control set as proof of more than it shows would be the most
expensive mistake available here.

- **It does not establish that a passing design binds.** There is no experimental binding
  data anywhere in this repo. Everything here is in silico, and a design that passes every
  threshold is a better hypothesis, not a result.
- **The positive control is a synthetic idealised helix, which is the easiest designed
  positive there is.** `C12` is a textbook backbone carrying a repeating designed
  sequence. It establishes that the scorer is not broken; it does not establish that the
  scorer discriminates among *real* designs, which is `C5` and is still blocked. The
  native complexes that were supposed to do this job cannot: predicted at
  `single_sequence` they fail for a reason unrelated to the score (finding 2).
- **The negative controls are scrambles, and a scramble is an easy negative.** The hard
  case is a design that looks fine and does not bind, and nothing in this repo can produce
  one without wet-lab data. So the suite can catch a scorer that is broken or merely
  generous; it cannot measure how well a working scorer ranks real candidates.
- **Five models, one seed, and a two-run variance figure that is not a distribution.**
  Repeating the whole Phase 1 control set gave `sc_rmsd_A` within ~0.03 Å and
  `plddt_binder` within ~0.04 — enough to say a 2 Å threshold is not being decided by
  noise, and not enough to be called variance. Nothing here varies the seed or the model
  set. 029 recorded all five AF2 models agreeing on one native complex (iPTM 0.850–0.857);
  nothing establishes that designs behave as tidily, and 037 found its own seed
  reproducibility was empirical rather than guaranteed, since `torch_deterministic` is
  left false upstream (`11bb5149`).
- **A passing `plddt_binder` means very little on its own.** Measured: 115 alanines
  score 95.81, above the real design. `sc_rmsd_A` is doing the discriminating, and any
  future metric added to the gate has to be shown to survive the poly-alanine control
  before it is trusted (§8.6).
- **Nothing here tests 037.** 040 consumes 037's output and `validation_pipeline.md` §10
  requires no changes to it. If 037 changes its contract, this suite notices only through
  the `ingest` section's fixtures.
