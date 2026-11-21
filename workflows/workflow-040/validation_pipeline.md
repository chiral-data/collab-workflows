# workflow-040: closing the validation loop — design spec

Status: proposal. Nothing here is built.

Context: workflow-037 generates binder candidates and stops. This specifies what has to
exist for those candidates to be *ranked* rather than merely listed.

## 0. Revision note

Revision 2, 2026-11-20. Revision 1 was written during the workflow-037 review and lived
outside the repo as `~/wf037-notes/validation_loop_spec.md`.

What changed, and why it matters before anyone starts building:

1. **The `single_sequence` default was wrong for the complex path, and the reason is not
   the one revision 1 assumed.** Workflow-029 measured it. §4 and §7.
2. **Three of revision 1's five open questions are closed** by evidence already in this
   repo rather than by new spikes. §11.
3. **Two of the three changes revision 1 asked for in workflow-037 have already
   landed.** Asking for them again would be exactly the defect this repo keeps fixing —
   a document describing a state that has since changed. §10.
4. **The data contract named a field vocabulary that does not exist.** Node 04 of 037
   writes three values for `designed_chain_source`, not two. §5.
5. **New §8** records the defect classes from 037 and 029 that 040 is positioned to
   repeat, each with the evidence and the rule it implies here. This is the section to
   read before writing a node.
6. **Phasing is now a decision rather than a plan.** Phase 1 is what gets built; Phase 2
   has a settled architecture and a gating spike. §12.

## 1. The problem, stated precisely

Workflow-037 emits designs whose quality nobody has measured. Its metrics are geometric
by declaration — contact counts at an 8 Å Cα–Cα cutoff, radius of gyration, net charge,
composition. Every one of them is computable on a design that folds into nothing and
binds nothing, and on such a design they all look ordinary. 037 knows this and says so:
`fd75d77b` added a broken-model check for exactly this reason, and `5a2d9b3f` put the
limits of its own test suite in writing.

The standard way to tell a real candidate from a plausible-looking one is to refold the
designed sequence with a structure predictor that never saw the design, and ask two
questions:

1. **Does the sequence encode the structure that was designed?** Superimpose the refolded
   binder on the designed binder and take the Cα RMSD. This is *self-consistency*. A
   sequence that does not recover its own backbone is not a candidate.
2. **Does the complex hold together?** Predict target+binder jointly and read the
   interface confidence. This is the metric that correlates with experimental success.

Neither is computable today, for a structural reason: **node 06 exports sequences, and a
sequence alone cannot be compared against a backbone that stayed behind.** The handoff is
lossy by construction.

Note what 037 has already fixed, because it changes what is left to do. Its report no
longer ranks on contact count (`e8326b55`), so there is no longer a misleading ranking to
displace — there is simply no ranking at all. 040's job is to supply the first honest one.

## 2. What is already known, and from where

Revision 1 treated this as greenfield. It is not. Three workflows in this repo have
already paid for knowledge 040 needs, and one of them has built most of node 04.

| Question | Answered by | Answer |
|---|---|---|
| What does `colabfold_batch` emit, and does it carry chain boundaries? | `029/04_fold_and_filter/script.py` | `*_scores_rank_001_*.json` carries `plddt` (flat per-residue array), `pae` (full N×N matrix), `iptm`, `ptm`. It carries **no chain boundaries**. 029 derives the split from a chain length it knows a priori. |
| How is `pae_interaction` computed in practice? | same | Mean over both off-diagonal blocks of the PAE matrix, concatenated. 040's §7 definition matches. |
| Can `single_sequence` mode produce usable interface metrics? | `2361277f` | No. See §4. |
| Does supplying the target's coordinates as a template rescue `single_sequence`? | `92cf05b5` | **No — measured.** This rules out the obvious cheap fix. |
| What happens if the target chain is trimmed to the segment that matters? | `8e4772a7` | It models far worse. A carved-out segment reaches ~42 pLDDT where the whole chain reaches 81, and the interface PAE floor moved 7.6 → 15.4 Å against a threshold of 10, putting that filter out of reach for *every* design: two 8×2 runs scored 2/16 before and 0/16 after. |
| Why does either of those break the interface metrics? | `8e4772a7` | **Interface PAE can be no better than the worse-modelled of the two chains.** This single rule explains both results and drives §4. |
| Is a separate monomer prediction needed for Cα RMSD? | `029/04` | Not necessarily — 029 takes its backbone RMSD from the complex prediction alone. See §11.4. |
| Does a `--dry-run` path for a prediction node already exist? | `026/02_predict/script.py:37` | Yes. Copy that pattern. |
| Why can 026 not predict a complex? | `d3d707ed` | `026/01_validate_fasta` checks every sequence against `ACDEFGHIKLMNPQRSTVWYX`, which has no `:`. The only syntax that encodes a complex is rejected before prediction. |

**029's node 04 is the single most useful artifact for whoever builds 040.** It already
computes `iptm`, `plddt_binder`, `pae_interaction` and a Cα RMSD, one design per
`colabfold_batch` invocation, in a temp directory. Read it before writing node 04 here.

It is prior art, not a dependency. 040 must not import from it or depend on it: it is
welded to 029's pipeline, reading ProteinMPNN sequence files and RFdiffusion backbones,
and it carries at least two defects 040 must not inherit (§8.7, §8.11).

## 3. Decision: a new workflow, not an extension of 037

**Recommendation, unchanged from revision 1: build `workflow-040`.** 038 and 039 are
taken by work in progress — both are untracked in the tree and open as draft PRs #280 and
#279 — so 040 is the next free number.

The deciding argument is that the handoff is only lossy because we chose to export the
wrong artifact. The `*_generated.cif` files already contain the designed binder backbone
*and* the target, in one file. A validation workflow that ingests **the CIFs instead of
the FASTA** has everything it needs. No information is lost across the boundary, so there
is no structural reason to fuse the two workflows.

Given that, separation wins on four counts:

1. **Reusability.** Self-consistency and interface scoring apply to binders from any
   source — RFdiffusion, RFD4, Chroma, hallucination, or hand-built. Welding them into
   037 throws that away.
2. **Different runtime shape.** Generation and prediction have different images, GPU
   profiles and failure modes. Users will want to generate once and re-validate several
   times under different settings. A single linear chain forces re-running generation.
3. **Reviewability.** 037 is six nodes and it took 60-odd commits to make those six
   honest. An eight-node chain with two heavy GPU stages is harder to reason about, and
   one failure in the last stage costs the whole run.
4. **No blast radius into 026.** We neither depend on nor modify someone else's workflow.

**Cost accepted:** validation is opt-in, so someone can still read 037's report without
it. 037 already states that its numbers rank nothing, so this no longer risks a reader
mistaking an ordering for a verdict.

### Rejected alternatives

- **Extend 037 to 8 nodes.** Rejected for the four reasons above. Reconsider only if
  Silva gains a way to compose workflows, which would remove the main objection.
- **Fix 026 to accept complexes.** Rejected as a dependency. 026 has its own users and
  loosening its validator to admit `:` is a behaviour change we should not make on their
  behalf while blocking our own work on it. 037 reached the same conclusion twice, in
  `d3d707ed` and `818954e9`, and the rule it settled on holds here: do not reach into
  another workflow to suit this one.
- **Reuse 029's node 04 directly.** Rejected. It reads 029's own artifacts and carries
  defects listed in §8. Read it; do not depend on it.

## 4. Decision: how the complex prediction gets its MSA

This is the substantive change from revision 1, which specified
`--msa-mode single_sequence` as the default and called it "both the right science for a
de novo binder (it has no homologs, so an MSA is noise) and the setting that keeps the
run offline."

The first half of that is true of the **binder**. It is false of the **target**, and the
target is what breaks.

### The measurement

029 ran a native complex through its node 04's exact flags, changing only `msa_mode`
(`2361277f`):

| Metric | `single_sequence` | `mmseqs2_uniref_env` |
|---|---|---|
| iPTM | 0.100 | 0.860 |
| interface PAE | 27.6 Å | 7.7 Å |
| receptor-side pLDDT | 38.4 | 83.7 |

All five AF2 models agreed (iPTM 0.850–0.857). Cost is roughly 2× wall clock per
prediction.

### The mechanism, which is not interface co-evolution

It would be easy to read that table as co-evolutionary signal across the interface. For a
*designed* binder there is no such signal to recover — the binder has no homologs, so
there is no binder MSA to pair against the target's.

What the table actually shows is the **target being modelled at ~40 pLDDT instead of
~84**, and `8e4772a7`'s rule then does the rest: interface PAE can be no better than the
worse-modelled of the two chains. A badly folded target caps the interface metrics
regardless of the design. `8e4772a7` is the same finding arrived at from the other
direction — trimming the target produced a ~42 pLDDT chain and moved the PAE floor out of
reach of the threshold.

Two consequences:

- **029's 0.100 → 0.860 is an upper bound for a designed complex, not a prediction.**
  It was measured on a native pair where both chains have deep MSAs. Nothing in this repo
  has measured the designed case.
- **The deficit is addressable without aligning the binder.** That is what makes the
  decision below available at all.

### Decision

**Phase 1 uses `single_sequence`, and this is correct rather than a compromise.** Phase 1
predicts the binder alone. A de novo designed binder is the one case where AF2 does well
without an MSA — designed proteins are idealised, which is why the self-consistency
convention exists at all. No MSA, no egress, no server dependency.

**Phase 2 aligns the target only, once per run, and reuses it.** The target chain is
public native protein, so aligning it discloses nothing. The designed binder is paired
against it unaligned, and never leaves the machine. Cost is one alignment per *run*
rather than per design, which matters because batches of designs are the entire point of
040.

This is gated on a spike (§11.1): whether `colabfold_batch` accepts a precomputed
per-chain a3m in which one chain is present only as its bare query sequence. The `.a3m`
input path and its `#`-prefixed chain-length/cardinality header are documented upstream;
the unaligned-partner shape is not confirmed, and it is the whole mechanism. **If the
spike fails, Phase 2 does not proceed on this design** — it falls back to the first
rejected alternative below, and that is a decision to re-make, not a detail to patch.

The alignment is its own node for a reason beyond caching: it makes network egress a
single declared boundary in one `job.toml`, which is what `92cf05b5` and `8a6dcebf` had
to retrofit onto 029.

### Rejected alternatives

- **Public MSA server for the whole complex** (029's current default). Known to work,
  zero mechanism risk, and the fallback if the spike fails. Rejected as the default on
  cost and reproducibility rather than on safety: it is one query per *design* against a
  shared, rate-limited community resource, at ~2× wall clock, and the server's databases
  drift, so scores are not reproducible across months.

  Worth stating plainly, because revision 1 overstated it: sending a designed sequence to
  an MSA server is a **custody and disclosure** problem, not really a biosecurity one. The
  server synthesises nothing, and the real gate is DNA synthesis screening, which is
  downstream and unaffected. What is true is that an unscreened, unpublished design ends
  up held by a third party, which sits awkwardly beside the notice 037 carries on every
  artifact (`7709809b`). That is a reason to prefer not doing it, not a reason to call it
  unsafe. If this fallback is taken, the egress is declared in the node's `job.toml` and
  in the README, as 029 now does.
- **Local MSA databases.** Correct science, zero egress, fully reproducible, and the only
  option that works air-gapped. Rejected on cost: ColabFold's database set runs to
  hundreds of GB, the setup lands on whoever runs the workflow rather than whoever builds
  it, and choosing a thinner subset is itself an unmeasured science question. 037 hit
  "no space left on device" on an 82 GB fetch (`8affe1e8`), so this estate has form.
  Revisit if an air-gapped deployment is ever asked for.
- **Target coordinates as a template, keeping `single_sequence`.** Rejected because it was
  **measured not to work** (`92cf05b5`). Recorded here so nobody re-proposes it.

## 5. Data contract

Everything crossing the boundary, and who owns it.

| Artifact | From | Required | Why |
|---|---|---|---|
| `*_generated.cif` | 037 node 03 | yes | designed binder backbone + target, in one file |
| `design_metrics.json` | 037 node 04 | strongly recommended | carries `designed_chain` and `designed_chain_source` |
| `generation_manifest.json` | 037 node 03 | optional | provenance chaining; carries `target_chains` |

**`design_metrics.json` is recommended for a specific reason.** Its rows carry
`designed_chain` per design. Without it, 040 has to re-derive which chain is the binder,
and the only available heuristic is chain order — which is exactly the assumption that
produced finding #27 in 037, where reversing two chains silently relabelled every metric
(`252ee072`). 040 must **read** the designed chain, not re-derive it.

### Read 037's vocabulary, do not invent one

Revision 1 said `designed_chain_source` is either `metrics` or `assumed-first-chain`.
That is wrong on both counts. 037's `resolve_designed_chain()`
(`04-analyze/script.py:248`) writes exactly three values, and 040 must handle each:

| Value | Meaning | What 040 does |
|---|---|---|
| `manifest-target-chains` | resolved by set difference against the target's recorded chain ids. The good path. | Trust it. |
| `single-chain` | the file holds one chain, so it is the design by construction — a monomer or `motif_scaffold` run. | Accept for self-consistency. **Refuse for the complex path**: there is no target to score an interface against. |
| `assumed-first-chain` | 037 declined to answer and fell back to chain order. | Proceed, say so loudly in the output, and carry the value into the report. |

Inventing a fourth value, or collapsing these three into two, would make 040's own output
unreadable against 037's. Read the string; do not re-map it.

### The fallback is the dangerous case, not the absent file

This is the sharpest lesson available from 037, and it is worth more than the rule it
amends. `ef74646e` found that the `assumed-first-chain` fallback had been in force on
**every real run** of 037, because node 03 recorded `target_chains` by listing every
`label_asym_id` in `target.cif` — which on 4ZQK is the target chain, the chain the contig
dropped, a sodium and two solvent groups. The designed chain was in that list, so the set
difference found no candidate and the fallback fired every time.

**And it happened to be right**, because the design does come out first. So nothing
failed, nothing looked wrong, and the only thing standing between a correct answer and a
silent mis-attribution was a chain order nothing checked.

The rule for 040: a fallback that produces the right answer is still a defect, and the
test that catches it is not a test of the fallback. It is a test that **the good path was
taken** on an ordinary input. §9 requires one.

## 6. Nodes

Five nodes for Phase 1, mirroring 037's shape. Phase 2 inserts one.

### Node numbering

Phase 2's target-alignment node is `01b-align-target`, not a renumbering of 02–05.
Renumbering a shipped workflow invalidates every saved `params.json`, every README
cross-reference and every fixture path, for a cosmetic gain. The letter-suffix form is
already used in this repo (`025/04a_run_vina`, `025/04b_run_gnina`). `01b` also sorts
where it runs: it depends only on node 01's inventory and runs once per run, not per
design.

### 01-ingest-designs

Read the CIFs and the optional metrics. For each design, resolve the binder chain and
the target chain(s), extract both sequences, and record residue numbering for each.

Emit `design_inventory.json`: one entry per design with binder chain id, binder sequence,
binder residue ids **in chain order**, target chain ids and sequences, and
`designed_chain_source` carried through verbatim from §5.

Fails when: no CIFs, a CIF has fewer than two chains for a complex run, or the named
designed chain is absent from the file.

Does **not** fail when `design_metrics.json` is absent — that is the stated-fallback path,
not an error. It fails loudly in the report instead (§8.2).

### 01b-align-target — Phase 2 only

Align the target chain once and cache the result. One target per run, so one alignment
per run.

This is the only node in the workflow permitted network egress, and it declares it in
its `job.toml`. The designed sequences are not inputs to this node, which is the
structural reason they cannot leak from it.

Fails when: the alignment cannot be obtained. A Phase 2 run without a target alignment
must not silently fall back to `single_sequence` and report interface numbers — §4 is the
whole argument for why those numbers would be meaningless.

### 02-build-prediction-inputs

Write, per design:

- `mono_<design>.fasta` — the binder sequence alone, one record. Feeds self-consistency.
  Phase 1 and Phase 2.
- `cplx_<design>.*` — one record, target and binder as two chains. Feeds interface
  scoring. Phase 2 only; the file format depends on the §11.1 spike.

One design per file, always. Never a batch file. The conflation bug in 026 — globbing
every design's scores together and taking the job id from whichever sorts first — is a
direct consequence of batching, and it produces a confident report about a design that
does not exist (`818954e9`).

**Record the chain order in the output, because node 04 depends on it.** The scores JSON
carries no chain boundaries (§2), so node 04 must slice the per-residue `plddt` array by
position. 029 slices `plddt[:binder_length]` because its binder is first; if 040 writes
`TARGET:BINDER` then the binder is the *tail*, and a node 04 that copies 029's slice
verbatim will average the target's pLDDT and call it the binder's. The order is a fact
this node decides, so this node records it and node 04 reads it. Node 04 additionally
asserts that the two recorded lengths sum to the array length (§8.1).

Validate here: residue alphabet, length bounds, chain count. Failing at this node costs
nothing; failing inside the predictor costs GPU minutes — which is why 037 made spec
validation its own CPU node (`37c23bb5`).

Note the inverse of 026's bug: a `:` is **required** in a `cplx_` input and **forbidden**
in a `mono_` one. Check both, rather than checking neither and relying on the predictor.

### 03-predict

Run the predictor once per input file. GPU.

| | Phase 1 (monomer) | Phase 2 (complex) |
|---|---|---|
| `--model-type` | `alphafold2_ptm` | `alphafold2_multimer_v3` |
| `--msa-mode` | `single_sequence` | target-only a3m per §4 |
| `--num-models` | 5 | 5 |
| `--num-recycle` | 3 | 3 |
| `--host-url` | not passed | not passed |

`--model-type` was missing from revision 1 entirely. It is not optional: the multimer
model is what produces `iptm` at all, and running a single chain through it is not the
monomer prediction anyone means. 029 sets `alphafold2_multimer_v3` explicitly.

**Do not pass `--host-url`.** 026 passes it unconditionally, defaulting to
`https://api.colabfold.com` (`026/02_predict/script.py:58`), which is how a workflow
reaches a third-party service without anyone deciding that it should.

Required assertions, each of which exists because its absence produced a silent wrong
answer somewhere in this repo:

- A `cplx_` prediction must come back with the expected number of chains. If it returns
  one, the `:` was dropped somewhere and the result is a monomer masquerading as a
  complex. Fail loudly.
- Weights must be baked into the image and the run must not reach the network.
  §11.1 settles whether that holds.

Per-design timeout, as in 037 node 03, and **by the mechanism 037 arrived at**, not by
`subprocess.run(timeout=...)`: generation spawns workers that go on holding the GPU after
the direct child is signalled, so `755d4af2` moved to `Popen` in its own process group
with a `kill_tree()` that escalates SIGTERM → SIGKILL. A prediction that stalls holds the
device exactly the same way. Budgets differ per path: a 5-model complex prediction is not
a monomer prediction.

Retain the predictor's output so a failure can be diagnosed rather than inferred
(`755d4af2` for the OOM case; `2359ffa4` for the opposite mistake — 029 discarded
ColabFold output unless the design passed, so a failing funnel could not be diagnosed
afterwards).

Partial failure is recorded per design and does not abort the batch. The manifest records
which designs have no prediction, and the report shows them as failed rather than
omitting them (§8.5).

### 04-score

Per design, compute:

- **`sc_rmsd_A`** — Cα RMSD, refolded binder vs designed binder, after Kabsch
  superposition. Correspondence is by *position in the chain*, established in node 01.
  **Both predictors and the generator renumber**: ColabFold renumbers from 1 per chain,
  and RFD4 renumbers the preserved target chain from 1 in its output, so for a target
  whose label numbering does not start at 1 — barnase in 1BRS starts at 3 — the generated
  structure's target numbering is already shifted from what 037 node 02 conditioned on
  (`ef74646e`). Map explicitly, assert the mapping (equal lengths, 1:1), never match on
  residue number.
- **`plddt_binder`** — mean pLDDT over binder residues, sliced by the order node 02
  recorded, with the length assertion of §8.1.
- **`plddt_target`** — the same for the target residues. **Diagnostic only: no threshold,
  no contribution to pass or fail.** This is not optional and it is not padding. 029 added
  it (`2359ffa4`) after its absence caused a wrong diagnosis: the interface filters were
  read as unreachable because the receptor had been predicted blind, and the receptor
  pLDDT needed to test that had to be computed by hand, outside the pipeline. Given §4,
  a badly-modelled target is 040's single most likely cause of bad interface numbers, so
  the number that distinguishes it is a first-class output.
- **`iptm`** — interface pTM from the complex prediction. Phase 2 only.
- **`pae_interaction`** — mean PAE over the off-diagonal blocks between target and binder,
  taken in both directions and averaged. Phase 2 only.

Merge with 037's metrics when present, **keyed by design name, never by position**
(`818954e9`; 026 takes its job id from whichever scores file sorts first).

A metric that could not be computed is `null` with a recorded reason. It is never 0.0 and
never `inf` — see §8.7, which is a defect 029 carries today.

### 05-filter-and-report

Apply thresholds, rank, emit `passing_sequences.fasta`, a metrics table and an HTML
report.

Every metric that gates a decision is a rendered column, not a log line and not an
embedded JSON field (`179003b9`: 037's node 04 warned about an unfolded design in the log
while the report rendered it as an ordinary row). Designs that failed prediction appear
with blank metrics and a stated reason.

Ordering: node 04 sorts, node 05 inherits, so there is one source of order
(`179003b9`, `e8326b55`). Sort **numerically**, not lexicographically, or `design_10`
lands between `design_1` and `design_2` (finding #51).

The report says what its order means. 040 is allowed to rank — that is the point of it —
but only on a calibrated metric, and until §12's Phase 3 runs, `pae_interaction` is
presented as a value with a provisional band rather than a gate (§7).

Renderers are isolated and provenance is drawn first. 037's report could lose its
screening disclaimer entirely when offline, because Plotly comes from a CDN and the
exception aborted the rest of the script, taking the run-info block with it
(`cb418341`).

The biosafety statement from 037 travels into every artifact here. A sequence that has
passed validation is *more* likely to be synthesised, not less, so the screening notice
matters more at this stage than at generation. Two mechanical constraints on the tag,
both learned the hard way in `7709809b`: no leading `;` comment banner, because Biopython
1.87 raises on comment lines; and no colon in the tag, because a colon in a record id or
description is read as a chain-break marker. The second is sharper here than in 037,
since 040 uses colons deliberately in its `cplx_` inputs.

## 7. Metrics, thresholds and their provenance

Thresholds are **defaults to calibrate**, not truth. Each is recorded in the output
alongside the value so a reader can re-filter without re-running.

| Metric | Default | Source | Caveat |
|---|---|---|---|
| `sc_rmsd_A` | < 2.0 | common self-consistency convention in the design literature | scale-dependent; a 30-mer and a 150-mer should not share a cutoff |
| `plddt_binder` | > 80 | standard AF2 confidence banding | measures fold confidence, not binding |
| `pae_interaction` | < 10 | Bennett et al. 2023, the strongest in-silico discriminator they report | **see below** |
| `iptm` | > 0.5 | AF2-multimer convention | weaker than `pae_interaction`; secondary |
| `plddt_target` | none | — | diagnostic only; never gates |

029 runs the same four thresholds with `min_iptm = 0.6`, `min_plddt_binder = 80.0`,
`max_pae_interaction = 10.0`. Where 040 chooses a different number it should say why, so
two workflows in one repo do not quietly disagree about what a good design is.

**The `pae_interaction` caveat is important and must be in the README, not just here.**
The published cutoff comes from protocols using *AF2 initial guess* — the predictor is
given the designed complex as a starting point, which sharply improves discrimination.
Vanilla ColabFold does not do this. The number will therefore not transfer unchanged.

And it is now clear that it will not transfer for a second, independent reason: §4's
measurement shows the metric's floor is set by how well the **target** was modelled. Any
cutoff carried over from the literature assumes a well-modelled target. `8e4772a7` is the
cautionary case — a target modelled at ~42 pLDDT moved the floor from 7.6 to 15.4 Å and
made a threshold of 10 unreachable for every design in the batch, which read as "all
designs are bad" rather than as "the target was never modelled". Until Phase 3's
calibration exists, present `pae_interaction` with a provisional band and always beside
`plddt_target`.

**Do not ship a threshold that has not been measured against something.** 037's
compactness check was first drafted at 1.5 and would have fired on three real monomer
designs sitting at 1.44, 1.46 and 1.49; running it across the archive before shipping is
what caught that (`fd75d77b`). The same discipline applies to all four numbers above.

**Narrow every claim to what was measured.** 037's compactness ratio was presented as a
folded-or-not verdict until `7968abeb` measured it against real structures and found it
flags a textbook α-helix at every length in the default binder range — spectrin at 2.27,
tropomyosin at 6.35, where a straight chain is 6.3. The threshold survived; the claim was
narrowed to "not compact for their length". For 040: `sc_rmsd_A` under a cutoff means the
sequence recovers its own backbone. It does not mean the design binds, and the report must
not imply it.

## 8. Defect classes from 037 and 029 that 040 must not repeat

037 took roughly sixty fix commits to make six nodes honest. Most of those fixes are
instances of a dozen or so classes, and 040 is positioned to repeat nearly all of them.
This section is the pre-build checklist.

### 8.1 Do not re-derive what an upstream node already decided

The repo's `CLAUDE.md` names this as the recurring defect class, and 037 fixed it at least
four times (`252ee072`, `ef74646e`, `179003b9`, `f8262f20`). Every instance looked like an
ordinary answer.

040's exposures, in order of likelihood:

- **The designed chain.** Read `designed_chain` (§5). Never chain order.
- **The chain order in the prediction input.** Node 02 decides it; node 04 reads it. The
  scores JSON has no chain boundaries, so a wrong slice silently averages the wrong
  chain's pLDDT. Assert `len(plddt) == n_target + n_binder` before slicing, and fail if
  it does not — that assertion is the only thing standing between a position assumption
  and a plausible wrong number.
- **Residue correspondence for RMSD.** Established in node 01, asserted in node 04.
- **Pass/fail.** Decided once, in node 04, and carried. Not recomputed in node 05 from
  the thresholds. `179003b9` is explicit about why: two nodes applying the same rule
  eventually disagree about the same design.

### 8.2 A fallback that gives the right answer is still a defect

§5 tells the story: 037's first-chain fallback was in force on every real run for weeks,
and was right every time. Carry `designed_chain_source` into the report, not just the
manifest, and test that the good path is taken on an ordinary input (§9).

### 8.3 Exit codes see nothing

All ten of 037's fixture entries exited 0 on nodes 01 and 02 before its suite existed.
The first run of real checks found twelve defects, **none of which moved an exit code**,
because the failure mode for a bad input shape is silence (`13153c89`).

040's equivalents: a prediction that returns a monomer where a complex was asked for; an
RMSD computed against a mis-mapped residue range; a pLDDT slice off by the target's
length. All exit 0. All produce a number. Assert against a measurement of the input, with
a parser that is not the one under test — 037's `probe_entry.py` parses `_atom_site` by
hand for exactly this reason.

### 8.4 A check that passes for the wrong reason is worse than no check

It is also evidence, which is what makes it worse. Four instances in 037:

- A log scan for component `NA` matched inside `NAG` and inside the string `nan`
  (`0e7c4881`).
- Checks matched prose in `report.html`, where node 05 explains its columns — so every
  column label appears in the file whether the column was built or not, and the checks
  passed either way (`bf84359b`). They read the report's own `DATA` payload now.
- A check matched the exact wording of a warning, so rewording it broke the check while
  the behaviour was unchanged (`7968abeb`).
- A check asked for a message that was itself wrong, and so would have failed the
  corrected node (`bf84359b`).

Rules for 040: assert on structured output and numbers, never on rendered prose; and
**verify every check is non-vacuous** by suppressing the behaviour and confirming the
check fails. 037 did this explicitly and recorded it in the commit message.

### 8.5 A design that failed must appear as failed, not vanish

026 drops what it cannot score. A report that omits failures reads as a batch where
everything worked. Node 05 renders failed designs with blank metrics and a stated reason.

### 8.6 Metrics computable on garbage look ordinary

This is 037's reason for existing in the shape it has, and 040 inherits the trap rather
than escaping it. `sc_rmsd_A` is computable on a prediction that converged on nothing;
`iptm` is a number whatever went in. The only defence is a negative control that must
score badly (§9).

### 8.7 Never default a missing metric to a number

029 does this today, in two places: `scores.get('iptm', 0.0)` and `pae_interaction =
float('inf')` when the matrix is absent. A monomer prediction has no `iptm`, so 0.0 reads
as a measured bad interface; `inf` fails the threshold, so a missing measurement is
rendered as a verdict. Both are indistinguishable from a real result downstream.

040: `null` plus a recorded reason. Never a sentinel that survives into a comparison.

### 8.8 Refuse non-finite parameters by name

`numeric_param()` in 037 accepted `nan` and `inf`. Neither is caught by a range check,
because every comparison against `nan` is false: it passed "must be positive" and passed
"outside the expected range" without a warning, then made every distance test false, so
the node reported no contacts on a structure that has plenty (`55edd8fc`).

040 has four float thresholds, all of them compared against measurements. Reuse 037's
`int_param()` / `bool_param()` / `numeric_param()` shape, including the detail that every
bad parameter must fail the same way — `975ab379` exists because 037 added two clean
parameters and left two raw `int()` casts ten lines above them, so a typo read as the
workflow crashing.

### 8.9 Declare each shared parameter once

029's node 05 declared `receptor_chain = "B"` while the workflow declared `"A"`, and
Silva's TUI and headless runner differ in whether a node's `job.toml` defaults are
synthesised for overlay. In the TUI, node 05 inherited the workflow's `"A"`,
`PARAM_BINDER_CHAIN` was never set, and **every design was scored against itself** —
`prodigy --selection A A` (`72e3a492`). The fix was to recognise that neither value was
ever configuration and make them module constants.

040: declare shared parameters at workflow level only. Anything not genuinely
user-configurable is a module constant, not a parameter with a default in two places.

### 8.10 Declare resources, and never have two dependencies write the same filename

`[resources] gpu = "required"` on node 03, `"none"` elsewhere, `schema_version = "1.1"`
(`7eee47a8`). Note the trade-off `97a16f7c` records: under silva 0.5.0, whose schema
cannot see `[resources]`, `use_gpu = true` is the only thing that passes the device
through, so the GPU node carries both keys and accepts one cosmetic note under 0.6.0.

And the contract bug an end-to-end run passed only by luck: 037's node 03 depended on
both 01 and 02, and both wrote `target.cif`. Silva copies every dependency's outputs into
the dependent's inputs, so which file survived was down to iteration order — harmless for
three of four tasks and fatal for the fourth (`7fba161c`). Declare the minimum edges, and
never let two dependencies write the same name.

### 8.11 Guard before rewriting, and make it idempotent

029's `swap_binder_to_chain_a` flipped A and B unconditionally and was not idempotent, so
re-running the node over a non-empty output directory swapped already-swapped files back
and put the receptor on chain A. This is not hypothetical: it fired during a GPU smoke
test, on a directory left behind by an earlier run (`4d1a523e`). The guard now asserts the
chain letters and lengths *before* the file is rewritten, so a refusal leaves the input as
it was.

040 rewrites nothing upstream, but node 02 and node 03 both write per-design files into
directories that may hold output from an earlier run. Re-running must be safe or refused.

### 8.12 Do not trim a gapped target, and disclose the fusion

`8e4772a7` in full: trimming a gapped receptor to the single segment carrying the epitope
was argued for on the grounds that ColabFold cannot express a gap inside a chain, so
folding the fragments concatenated presents two lobes as bonded. The argument was never
measured against the alternative, and it was wrong — ~42 pLDDT trimmed versus 81 whole.

040 folds the target chain whole and **warns that each internal break is fused rather than
cutting around it**. The fused junction is a real defect, kept knowingly, and the README
says so. 037's targets come from RCSB, so disordered gaps are ordinary, not exotic.

### 8.13 A guard that cannot be shown refusing is not a guard

And one that refuses everything is not one either. 037's suite has a `rejections` section
asserting both the exit status and that the refusal gives the stated reason, plus one
entry proving a legitimate override is still accepted (`4b1a611a`, `b6b741e9`).

### 8.14 A bundled example must be what it says it is

029's README described its benchmark as IL-7Rα and cited hotspots to Watson et al. 2023.
PDB 2B5I is the quaternary IL-2 receptor complex; chain A is interleukin-2, and the string
"interleukin-7" does not appear in the file. The documented hotspots were residue numbers
from a different protein, and they "passed validation because they exist, not because they
meant anything" — one of them was not in the deposited coordinates at all (`11ddf044`).

040's positive controls carry more weight than a demo input: they are what calibrates the
thresholds. Verify by measurement which chain is which and what the binder actually is,
and record the measurement beside the fixture, as 037's `fixtures.json` does with its
`evidence` strings.

### 8.15 Documentation is part of the change

The repo's `CLAUDE.md` is explicit, and 037 has three commits whose entire content is
catching up (`b1456f2c`, `4757ebfe`, `b0baa35e`). `b1456f2c` is the instructive one: four
places in one README described the node 06 handoff, three correctly, and the wrong one was
the paragraph a reader meets first. Following it hands 026 every design at once and gets a
conflated report — the exact failure two earlier commits existed to prevent.

For every behaviour change in 040: the node README section, `tests/README.md` and its
known-gaps list, `tests/fixtures.json`, the node's `job.toml` inputs/outputs/hints, the
module docstring, and `workflows/CHANGELOG.md` + root `README.md` + `apps/applications.json`
when the workflow is registered. Same commit.

### 8.16 When a test fails, suspect the harness first

Building 037's downstream fixtures, the target was placed at a fixed 14 Å, which at
design length 30 left a 2.6 Å gap — so all three designs came out flagged as clashing and
the case appeared to show node 05 failing when node 05 was correct (`669f4ae9`). The
separation is now measured from the design's own extent.

040's synthetic predictions have the same shape of risk: a fixture that "proves" the
scorer wrong is more likely to be a fixture with the wrong geometry.

## 9. Test plan

**The suite itself is specified in `tests/README.md` and `tests/fixtures.json`, and it is
the merge gate.** 58 experiments across two tiers — a CI tier with no GPU, and a GPU tier
run by hand with its numbers recorded in the PR. This section is the summary; those files
are authoritative, and a disagreement between them is a bug in this one.

Mirrors 037's suite: real inputs, measured properties, assertions against the input
rather than against an exit code. Plus the thing 037's suite could never have.

### Positive control — the piece that has been missing

Take a **native, experimentally validated** complex from the PDB and feed it through 040
as though it were a design: the native binder chain's own sequence, against its own
backbone.

- 4ZQK (PD-1/PD-L1) — already 037's worked example, and already a 037 fixture.
- 1BRS (barnase–barstar) — a textbook high-affinity complex, also already a 037 fixture.
  Note barnase's label numbering starts at 3, which makes it the right case for the
  renumbering assertion in §6/node 04.

A native binder on its own backbone **must** score well: `sc_rmsd_A` low, `plddt_binder`
high, and in Phase 2 `pae_interaction` good. **If 040 scores a real complex badly, 040 is
wrong.** This is the first test that can falsify the scoring itself rather than the
plumbing, and it is what calibrates §7's thresholds.

Verify each entry by measurement before trusting it, per §8.14, and record the
measurement. 1BRS holds three copies of the complex across six chains — 037 found that
five of the six vanished silently from its own pipeline (`53d1d460`) — so "chain A" is not
a sufficient description of what is being tested.

### Negative control

The same backbone with a scrambled or poly-alanine sequence. Must score badly. A filter
that passes everything is not a filter, and unlike the positive control this one catches a
scorer that is simply generous (§8.6).

### Diagnostic control — new, and specific to §4

Run one positive control through the complex path **three ways**: `single_sequence`,
target-only MSA, and full MSA. Read `iptm`, `pae_interaction` and **both** chains' pLDDT
each time.

This is the experiment that settles §4 empirically rather than by argument. It is one
structure and three predictions. It tells you whether the target-only alignment recovers
the full-MSA numbers, and `2359ffa4` exists precisely because the per-chain pLDDT needed
to interpret such a result was missing once before.

### Plumbing, no GPU

Reuse `workflow-037/tests/make_designs.py`. Its serpentine walk through a cubic lattice
puts every consecutive Cα exactly 3.8 Å from the last while staying globular, and it
already builds folded, extended, broken-backbone, interpenetrating and helical shapes —
so a synthetic prediction with a *known* RMSD against a known backbone is available
without a GPU.

Node 03 needs a `--dry-run` path so nodes 01, 02, 04 and 05 run in CI without one. Copy
026's (`026/02_predict/script.py:37`).

### Specific regressions to pin

Each of these is a defect that has actually occurred in this repo:

- A complex prediction returning one chain must fail, not score. (026's bug.)
- Metrics must never be matched across designs by position. (026's conflation bug.)
- A pLDDT slice must be rejected when the array length does not match the recorded chain
  lengths. (§8.1; 029 slices binder-first and 040 may not be.)
- Absent `design_metrics.json` must produce a stated fallback, not a silent guess. (#27.)
- `designed_chain_source == manifest-target-chains` on an ordinary input — i.e. the good
  path was taken, not merely available. (§8.2; this is the test `ef74646e` did not have.)
- `designed_chain_source == single-chain` must be refused on the complex path, not scored
  against an absent target. (§5.)
- A design whose prediction failed must appear in the report as failed, not vanish.
- A missing metric must be `null`, never 0.0 or `inf`. (§8.7.)
- `nan` and `inf` must be refused for every threshold parameter. (§8.8.)
- Re-running node 02 or 03 over a populated output directory must be safe or refused.
  (§8.11.)
- Numeric, not lexicographic, design ordering at twelve or more designs. (#51.)

## 10. Changes required in workflow-037

**Revision 1 asked for three. Two have since landed.** Verify before touching anything —
re-asking for a completed change is the same defect class as §8.15, pointed the other way.

1. ~~Node 05 must stop implying a ranking it cannot support.~~ **Done**, in `e8326b55`.
   Designs are listed by name, clashing ones sunk, geometry-flagged ones sunk into a
   separate bucket (`bf84359b`), and the report states that the order is not a ranking
   (`README.md:280`). Nothing left to do.
2. ~~Document the handoff as the CIFs, not the FASTA.~~ **Done**, in `252ee072` and
   `b1456f2c`. The README now says in four places that a validation step reads
   `*_generated.cif` and `design_metrics.json` rather than the exported FASTA
   (`README.md:58`, `:92`, `:259`, `:725`). The FASTA export stays, as it has other
   consumers.
3. **No changes to nodes 01–04.** Still true. The contract 040 needs is already there, and
   §5's three-value `designed_chain_source` vocabulary is what 037 already writes.

So 040 requires **no changes to 037 at all**. If a gap turns up during the build, it goes
in this section with its evidence rather than being fixed opportunistically across a
workflow boundary.

One thing worth adding to 037 only if Phase 2 proceeds: nothing in 037 records whether a
target chain carries internal gaps, and §8.12 makes that fact load-bearing for the complex
prediction. 040 node 01 can measure it from the CIF itself, so this is a preference, not a
dependency.

## 11. Open questions

Revision 1 listed five. Three are now closed by evidence already in the repo.

### 11.1 Open, and gating: does `colabfold_batch` accept a partially-aligned complex a3m?

Specifically: an a3m whose `#` header declares two chains, where the target carries a deep
alignment and the binder is present only as its bare query sequence. This is the entire
mechanism of §4's decision. The `.a3m` input path and the header format are documented;
this shape is not confirmed.

**Spike: ~1 hour.** If it fails, §4's first rejected alternative becomes the design, and
that is a decision to re-make rather than a detail to patch.

Folded into the same spike, since it needs the same image: **does the prediction run fully
offline with baked weights?** Weights must be in the image, and the node must not pass
`--host-url`. Revision 1's question #1, still open, now cheaper to answer alongside this
one.

### 11.2 Closed: what does `colabfold_batch` emit for PAE, and does it carry chain boundaries?

`*_scores_rank_001_*.json` carries `plddt` (flat per-residue array), `pae` (full N×N
matrix), `iptm` and `ptm`. It carries **no chain boundaries**, which is why §6 makes the
chain order an explicit output of node 02 and an assertion in node 04. Verified in
`029/04_fold_and_filter/script.py:203-223`.

Residual, and small: confirm the key names for the monomer model (`alphafold2_ptm`), which
has no `iptm`. §8.7 is the rule that keeps a missing key from becoming a number.

### 11.3 Closed, and reframed: how much discrimination is lost without AF2 initial guess?

The mechanism is now known. The dominant term is not the initial guess but **how well the
target is modelled** (§4), and a badly-modelled target caps the interface metrics
regardless (`8e4772a7`). Supplying the target's coordinates as a template was measured not
to rescue a `single_sequence` run (`92cf05b5`), which rules out the cheap fix.

What remains is a measurement, not a question: §9's diagnostic control. If native
complexes do not separate from scrambled ones with a well-modelled target, vanilla
ColabFold is insufficient and initial-guess support becomes a prerequisite rather than a
refinement.

### 11.4 Partly closed: is a monomer prediction needed at all?

029 takes its Cα RMSD from the complex prediction alone, with no separate monomer run, so
this is possible in principle.

It does not change Phase 1, which has no complex path and therefore needs the monomer run
regardless. It becomes a Phase 2 question with a concrete test: compare `sc_rmsd_A`
computed from the monomer prediction against the same quantity read off the binder chain
of the complex prediction, on the positive controls. If they agree, drop the monomer run
and halve the GPU cost. Expect them not to agree exactly — the monomer number is the
cleaner self-consistency signal, since the complex prediction has the target to lean on.

### 11.5 Open: GPU budget per design

029 measured roughly 2× wall clock for MSA mode versus single-sequence, but not absolute
cost per design at 5 models. 037 measured 35–50 s per replicate for generation, which is
the wrong unit here.

Measure before promising a batch size. It sets the per-design timeout of §6/node 03, and
037's experience says the default and the documented production setting must be checked
against each other — `b0baa35e` found 037's own documented production advice of "100+
replicates" landing on its own 3600 s default, killing the job hours in.

## 12. Phasing

**Phase 1 — self-consistency only. This is what gets built.**

Nodes 01, 02 (mono only), 03 (`alphafold2_ptm`, `single_sequence`, offline), 04
(`sc_rmsd_A`, `plddt_binder`, `plddt_target`), 05.

No complex prediction, no `iptm`, no `pae_interaction`, no MSA, no network. `single_sequence`
is the correct setting here rather than a compromise (§4). This alone turns 037's output
from unlisted into filterable, and it sidesteps every open question except §11.1's offline
half and §11.5.

Carry `plddt_target` from the start even though Phase 1 has no interface metric for it to
explain. It costs nothing, and it is the number that makes a Phase 2 result interpretable.

**Phase 2 — interface scoring. Architecture settled, gated on a spike.**

Adds `01b-align-target`, the `cplx_` path through nodes 02–05, and `iptm` +
`pae_interaction`. Target-only MSA per §4, which is a decision, not an option left open.

Gate: §11.1. If the a3m shape is not accepted, stop and re-make the §4 decision rather
than quietly falling back to the public server.

**Phase 3 — calibration.**

Run §9's controls, including the three-way diagnostic. Set thresholds from measurement,
and write what was measured into the README beside each number. Until this phase lands,
`pae_interaction` is a provisional band and not a gate (§7).

Phase 1 is the smallest thing that is genuinely useful. Ship it before starting Phase 2.
