# Tier 2, Phase 2 controls — measured

Real ColabFold multimer predictions of a native complex and the same complex with the
binder's residues shuffled, through the whole workflow including the interface path, on the
GPU, with no network. The numbers below are what came out, and they say the control was
built on an assumption that is false.

Reproduced by `bash tests/run_controls_phase2.sh`. Judged by
`tests/report_controls_phase2.py`. **`bash tests/falsify_controls_phase2.sh --synthetic`
shows that judge refusing**, which is the only reason to believe it when it agrees.

**Verdict: Gate B Phase 2 NOT MET** — `iptm` does not separate. The positive scores
**0.0900** and the negative **0.1700**, a margin of **−0.0800**: the shuffle scores
*better* than the native complex.

**The cause is not 040's scoring. It is this control's target alignment**, and the workflow
said so itself before the judge did. See finding 1.

## Provenance

| | |
|---|---|
| Host | `amd-c12-1`, NVIDIA GeForce RTX 3060 Ti 8 GB |
| Image | `colabfold:2026_10_08` (built from `workflow-026/Dockerfile`) |
| Run | `~/wf040-controls-phase2/20261122-065331`, 2026-11-22 |
| Monomer arm | `--model-type alphafold2_ptm --msa-mode single_sequence --num-models 5 --num-recycle 3` |
| Complex arm | `--model-type alphafold2_multimer_v3 --msa-mode mmseqs2_uniref_env --num-models 5 --num-recycle 6` |
| Network | `docker run --network none`, `host_url_passed=False`, `network_used=False`, every prediction offline |
| Target MSA | 9 sequences from `file`, sha256 `f9c719180448b2dc` — **built offline by `tests/make_target_msa.py`** |
| `binder_aligned` | `False` — the designed chain is given no alignment, which is §4's rule |
| GPU | 1 MiB used before, 1 MiB after — the device was left free |
| Wall clock | monomer 70.2 / 70.6 s, complex 260.1 / 220.6 s |
| Disk | 18G free before, 18G after, against an 8G floor |
| Entry | `1BRS`. **Binder** = chain A, barnase, 108 res resolved, ids 3..110, no gaps. **Target** = chain D, barstar, 87 res resolved, ids 1..89, gap 63→66. Both read off the run's own probe, not from memory |
| Chain relabelling | `make_controls.py` writes the target as chain A and the design as chain B, so the *output* chain A is barstar. The probe's chain letters and the prediction's are not the same letters — see the correction note below |
| Nodes run | 037 node 04, then 040 nodes 01, 01b, 02, 03, 04, 05 |

All four predictions completed: `status=ok` for both designs in both modes. Nothing here is
a crash being read as a result.

## What came out

| | role | `iptm` | `pae_interaction` | `plddt_target` | `plddt_binder_in_complex` |
|---|---|---|---|---|---|
| `C8` | positive, native | 0.0900 | 25.7820 | **38.2800** | 34.9900 |
| `C9` | negative, scrambled | **0.1700** | 23.6980 | **39.1900** | 30.0500 |

Separation, worst positive against best negative:

| metric | margin | gates? |
|---|---|---|
| `iptm` | **−0.0800** | **yes** |
| `pae_interaction` | −2.0840 | no |
| `plddt_target` | −0.9100 | no |
| `plddt_binder_in_complex` | +4.9400 | no |

Only `iptm` gates. `pae_interaction` ships as a provisional band in §7 that no control has
calibrated, and a threshold nobody has measured must not be the thing that fails a run.

The monomer arm ran in the same job and also failed, as expected for native chains at
`single_sequence`: `C8` `sc_rmsd_A` 11.5527 at 43.12 `plddt_binder`, `C9` 16.9880 at 49.19.
Those reproduce Phase 1's `C2` (11.57 Å, 43.09) to within 0.02 Å and 0.03 pLDDT, which is
the run-to-run noise Phase 1 measured. The monomer path is behaving identically; nothing
about this run is anomalous apart from the interface numbers.

## Finding 1: the target was never modelled, so no interface number means anything

**Node 04 diagnosed this itself, in the run, before any judge was consulted:**

> 2 design(s) sit against a target modelled below 70.0 pLDDT: C8 (38.28), C9 (39.19). No
> interface number computed against a target that was not modelled means anything, and
> given §4 this is the most likely cause of bad Phase 2 figures.

38.28 is not merely low. **§4's measured table records a receptor-side pLDDT of 38.4 at
`msa_mode = single_sequence`** — and this run passed `mmseqs2_uniref_env` with a
precomputed 9-row alignment. Getting the single-sequence number while supplying an
alignment means the alignment carried no information.

Measured directly, it carries none:

| | |
|---|---|
| rows in `target_msa.a3m` | 9 (one query + 8 homologs) |
| query length | 87 |
| positions each homolog differs at | **1 of 87**, every one |
| identity to query | **98.85%**, every row |

`make_target_msa.py` builds each homolog by one point substitution, which is exactly what
it says it does. Eight sequences that are 98.85% identical to the query contain no
covariation, and covariation is the entire mechanism by which an MSA improves a
prediction. So AF2 had the query and nothing else, and the target folded as if at
`single_sequence` — which is what the 38.28 says.

And per `8e4772a7`, **interface PAE can be no better than the worse-modelled chain**. With
the target at ~38 pLDDT in *both* arms, the interface metrics in both arms are measuring
the same floor. The 0.09 against 0.17 is not the native losing to the shuffle; it is two
numbers drawn from noise, and the sign of the difference is arbitrary.

## Finding 2: the refuted assumption, and where it was written down

`run_controls_phase2.sh` and `report_controls_phase2.py` were both built on the premise
that **a depressed baseline still permits separation** — that if the native and the shuffle
are modelled under the same handicap, the native would still come out ahead by enough to
read. That premise was never measured. It has now been measured, and it is false.

This is §8.1's failure mode in the mirror: nothing was re-derived, something *unknown was
asserted* and a gate was built on it.

**And the contradicting measurement was already in the spec, in the very section that
specified this control.** §11.1's table, same a3m, same seed, both arms offline:

| `msa_mode` | `iptm` |
|---|---|
| `single_sequence` | 0.270 |
| `mmseqs2_uniref_env` with the synthetic a3m | **0.120** |

§11.1 even wrote the conclusion out — *"a tiny, synthetic, heavily gapped MSA can be worse
than none"* — and then declined to draw it, on the narrow and correct ground that the
finding there was the *difference* between the arms rather than its sign. Narrowly right.
It also left a measured warning sitting unused two sections above a control built on the
opposite being true.

§4 held the other half: the 38.4 receptor-side pLDDT at `single_sequence` that this run
reproduced as 38.28. So this was not a gap in the evidence. It was a gap in reading
evidence already collected, which is the harder kind to notice, because the premise felt
like background rather than a claim and so was never looked up.

It is the same shape as Phase 1's finding 2, where native positives were run at
`single_sequence` and §4 had said in advance why they could not work. Phase 1 recorded that
lesson about the *binder*. The identical mistake was then made about the *target*, one
phase later, in a script written after that finding was filed. Recording it here because a
lesson that only applies to the example it was learned on is not a lesson.

**What this means for the verdict.** Gate B Phase 2 is NOT MET, and it must stay NOT MET —
the gate is honest and the judge is right to refuse. But NOT MET here is a statement about
the *control*, not about 040's interface scoring, and the distinction matters because the
two have different remedies:

- If 040's scoring were wrong, 040 would not merge.
- The control's alignment being uninformative is fixed by giving it a real alignment.

Phase 1 had `C12` to tell those two explanations apart. **Phase 2 has nothing equivalent,
and that is the real gap**: there is no Phase 2 positive control whose target is modelled
well enough for its interface numbers to mean anything. Until there is, Phase 2 cannot
distinguish a broken scorer from a broken control, and no amount of re-running it offline
will change that.

## Pre-registered: what the real-MSA run must produce

Written **before** that run, and deliberately so. The whole reason Phase 2's first attempt
cost a day is that its central premise was never written down where it could be checked
(§8.17). So the next run's expectations go here first, with their provenance, and the
result is read against them rather than explained after the fact.

Every number below is either **measured in this repo** with a commit to point at, or a
**stated convention**, or **explicitly unknown**. Nothing is a guess dressed as a target.

| | prediction | where it comes from | status |
|---|---|---|---|
| `plddt_target` | **rises from 38.28 to > 70** | §4's table: 38.4 → 83.7 changing only `msa_mode` on a native complex (029, `2361277f`). 70.0 is node 04's own warning threshold | **measured, and this is the primary claim** |
| `iptm`, `C8` | **> `C9`'s, by a margin that is not noise** | §8.6: the gate is separation, not an absolute | **measured requirement** |
| `iptm`, `C8` absolute | somewhere in **0.100 – 0.860**, and nearer the bottom is not a failure | 029 got 0.860 on a native pair where **both** chains had deep MSAs. Here the binder is given none, by §4's design (`binder_aligned=False`) | **bounded, not predicted** |
| `pae_interaction`, `C8` | between **27.6 Å and 7.7 Å** | same table, same caveat | **bounded, not predicted** |
| `iptm` ≥ 0.5 as "confident interface" | **not used as a gate** | AF2-multimer convention, recorded in §7's threshold table | **convention, uncalibrated** |
| `C9` | stays low | it is the same backbone with the residues permuted, and 037's node 04 reports it as **identical** to `C8` on every column it has (§8.6's measured table) | **measured requirement** |

### The primary claim, and why it is the one that matters

**`plddt_target` must come off the floor.** That is mechanical, not hopeful: §4 measured a
native receptor going 38.4 → 83.7 on `msa_mode` alone, and this run's 38.28 is that same
floor. If a real alignment does not move it past 70, then the fetch did not do what it was
for, and **no interface number from that run is interpretable either** — `8e4772a7` again,
the interface can be no better than the worse-modelled chain.

So the reading order is fixed in advance:

1. **`plddt_target` > 70?** If no, stop. The run says nothing about separation, exactly as
   this one did, and the conclusion is about the alignment, not the scorer.
2. **Only then, does `iptm` separate?** If yes, Gate B Phase 2 is met. If no, *that* is the
   first evidence this repo has ever had that 040's interface scoring does not work, and it
   would be a finding about the workflow rather than the control.

Step 1 is what this file's first version lacked the discipline to require, and it is the
reason the 0.09-against-0.17 result was briefly readable as "040 cannot score an
interface". It never was.

### What would falsify the diagnosis in finding 1

Finding 1 claims the uninformative alignment is why Phase 2 failed. That claim is wrong if
a real alignment lifts `plddt_target` past 70 **and** `iptm` still fails to separate. That
is a real possibility and it is written here so it cannot be quietly reinterpreted later:
the diagnosis is a hypothesis with a stated refutation condition, not a conclusion.

## Run 2: the real alignment. The primary claim held; the control still cannot tell

Run `~/wf040-controls-phase2/20261122-115101`, 2026-11-22, with `--allow-egress`. Node 01b
submitted the 87-residue barstar target to `api.colabfold.com` and received **3600
sequences (3599 homologs)**, sha256 `479c3311e183d652`, `source: server`,
`network_used: true` — against run 1's 9 synthetic rows. Egress was scoped to node 01b;
every other node stayed on `--network=none`.

| | `iptm` | `pae_interaction` | `plddt_target` | `plddt_binder_in_complex` |
|---|---|---|---|---|
| `C8` positive, native | 0.1900 | 20.8920 | **85.4200** | **38.5500** |
| `C9` negative, scrambled | **0.2500** | 18.7060 | **85.9800** | **43.2600** |

Monomer arm, unchanged by the target alignment and reported for reproducibility:

| | `sc_rmsd_A` | `plddt_binder` | run 1 | Phase 1 |
|---|---|---|---|---|
| `C8` | 11.3239 | 43.47 | 11.5527 / 43.12 | `C2` 11.5737 / 43.09 |
| `C9` | 16.9864 | 49.22 | 16.9880 / 49.19 | — |

### Against the pre-registration, in the order it fixed

**1. `plddt_target` > 70? YES, and almost exactly as predicted.** 38.28 → **85.42**. §4's
table said 38.4 → 83.7 changing only `msa_mode` on a native complex. The real alignment
did precisely what §4 said it would, on a different entry, two years of commits later.
Finding 1's narrow diagnosis — that the 9-row synthetic a3m floored the target — is
**confirmed by measurement**.

**2. Does `iptm` separate? NO.** 0.1900 against 0.2500, margin −0.0600. The scramble again
scores higher than the native, as it did in run 1 (0.0900 against 0.1700). Two independent
runs, different target alignments, same sign.

### My refutation condition triggered, and my pre-registration was too weak

This file's pre-registration said, in as many words:

> Finding 1 claims the uninformative alignment is why Phase 2 failed. That claim is wrong
> if a real alignment lifts `plddt_target` past 70 **and** `iptm` still fails to separate.

Both halves happened. By the letter of what was written before the run, finding 1's
diagnosis is refuted and this is the first evidence against 040's interface scoring.

**That condition was defective, and the defect is mine.** `8e4772a7`'s rule — the
interface can be no better than the worse-modelled chain — is **symmetric over the two
chains**, and the pre-registration applied it to the target only. It required one chain to
be modelled and said nothing about the other. So it could not distinguish "the scorer
failed" from "the other chain is on the floor", which is the single distinction the whole
exercise existed to make.

And the other chain is on the floor. **`plddt_binder_in_complex` is 38.55 and 43.26.**
That is not a new number: Phase 1 measured native barnase at `single_sequence` as **43.09
`plddt_binder`**, and this run's own monomer arm gives **43.47**. The binder is sitting at
its own no-alignment floor, in isolation and in complex, to within a pLDDT point.

So: writing the premise down was right and insufficient. **A premise can be written down
and still be wrong in the same way the unwritten one was** — §8.17 recurring inside the
correction for §8.17. What would have caught it is not more discipline about writing
things down; it is deriving the condition from the rule (`8e4772a7` is symmetric, so the
condition must be) rather than from the symptom in front of me (run 1's floored target).

### Why the binder is floored, and why that is not a bug

§4 gives the binder **no alignment at all**, by design: *a de novo binder has no homologs,
and designed proteins are idealised*. `binder_aligned=False` in every run above. That is
correct for the input class this workflow exists for.

`C8`'s binder is **native barnase**. A native immunoglobulin-like domain cannot be folded
from its sequence alone — which is Phase 1's finding 2, verbatim, about `C1`/`C2`. So
`C8`/`C9` are mis-specified in the exact mirror image of the controls Phase 1 already
retired: Phase 1 used a native **binder** as a positive and had to replace it with a
designed one; Phase 2 then fixed its native **target** with an alignment and left the
native binder in place.

**Three times now, the same error about a different chain.** §8.17's two practical rules
were both followed this round — the premise was written down, with provenance — and the
error still landed, because the rules were about *recording* assumptions rather than about
*deriving* them from what is already measured.

### What this does and does not establish

- **`plddt_target` is fixed, and the egress path works end to end.** 3600 real sequences,
  the target modelled at 85, node 01b's guard shown refusing, the alignment's digest and
  provenance recorded. That part of Phase 2 is done.
- **040's interface scoring is still unvalidated — not disproven.** There has still never
  been a Phase 2 control with **both** chains modelled. Nothing here licenses a statement
  in either direction, and the judge now refuses to make one.
- **The negative outscoring the positive is reproducible** across two runs and both
  alignments. On current evidence that is a property of *native barnase given no
  alignment* versus *a scramble of it*, not a property of 040. Phase 1's `C4` is the
  precedent: 115 alanines scored 95.81, above the real design. A sequence encoding nothing
  can be modelled confidently (§8.6).

### The judge now distinguishes the two outcomes

"NOT MET" was doing two jobs: *the scorer failed* and *this control cannot tell*. Run 1's
verdict was briefly readable as the first when it was always the second, and that reading
is what cost a day.

`report_controls_phase2.py` now applies the modelled-chain floor to **both** chains of
**every** control before it reports separation at all, and returns a distinct
**NOT INTERPRETABLE** verdict with exit code 2 when a chain is below it. Re-judging run 2
gives:

> Gate B Phase 2: NOT INTERPRETABLE — this control cannot tell
>   - `C8` `plddt_binder_in_complex` is 38.55, below the 70.0 floor, so that chain was not
>     modelled and no interface number computed against it means anything (`8e4772a7`).
>   - `C9` `plddt_binder_in_complex` is 43.26, below the 70.0 floor …

Separation is not printed at all in that state, so there is no margin lying around to be
quoted out of context. An eighth falsifier mutation (`binder-chain-not-modelled`) floors
one chain and requires this verdict rather than MET **or** NOT MET.

**That change introduced a defect of its own, caught by running it.** The new verdict
printed in an output shape of its own, and `falsify_controls_phase2.sh` quotes the first
`NOT MET:`-or-`  - ` line as the reason a mutation was caught — so the eighth mutation was
recorded as refused with an **empty reason**, a pass that looks like a pass and says
nothing. The same thing had already happened once in this script, for the early-exit
versus bulleted shapes. Two output formats were one too many; a third was the point to fix
the writer instead of adding a third grep to the reader, so every refusal this judge makes
now prints as a `  - ` bullet.

### What Phase 2's positive control actually requires

Stated as a specification, because this is the third control to fail on it:

**Both chains must be modellable under the flags §4 mandates** — the target with an
alignment, the binder with **none**. That means the binder must be a **designed** sequence,
not a native one. It is Phase 1's `C12` move, one phase up: `C12` recovered its own
backbone at 0.93 Å / 95.05 pLDDT from sequence alone precisely because it was designed.

And for a *positive*, the designed binder must be one that is **known to bind** its target,
with a solved complex to measure against. That is a stronger requirement than any control
in this repo has met:

- `C12` is designed but has **no interface at all**.
- `C5`'s 037 designs are designed, but nothing establishes that any of them binds, so they
  cannot serve as a positive — only as unknowns.
- `C8`/`C9` have a real interface but a native binder.

So the next step is identifying a deposited complex whose binder chain is an
experimentally validated de novo design and whose target is native. That is a literature
and PDB question, not a code question, and it is the first thing Phase 2 needs. Until such
an entry is in hand and probed the way `make_controls.py --list` probes every other entry,
Phase 2 has no positive control and Gate B Phase 2 cannot be met by any amount of GPU time.

## Finding 3: the egress path, now exercised

**`C8` and `C9` needed a real target MSA, which needed one fetch from
`api.colabfold.com`.** That fetch has been made, with `--allow-egress`, and run 2 above is
its result.

What crossed the boundary: **one POST carrying the 87-residue barstar target sequence**,
public deposited protein from 1BRS. `--allow-egress` switches `PARAM_MSA_SOURCE` to
`server` and opens `--network=bridge` **for node 01b only** — every other node in the run
stayed on `--network=none`, and the script resets it immediately after. The flag exists so
that reaching the network is a decision somebody made rather than something a test script
did on its own.

Recorded provenance: 3600 sequences, `source: server`, `network_used: true`, sha256
`479c3311e183d652`, `host_url_passed: false`. The alignment's digest is in
`target_alignment.json` so the run can be told apart from an offline one without reading
the log.

**And the boundary was checked before it was used.** The claim that a designed binder could
never cross it was false — three files said the designed sequences were not inputs to node
01b, and `design_inventory.json` carries `binder_sequence` for every design and is 01b's
only input. It is now enforced by `assert_no_designed_sequence()` immediately before the
POST, with the fixture case `egress-refuses-a-designed-sequence` showing it refuse. See
that commit; it is the one piece of this work that was about disclosure rather than about a
number.

**Separation is not measurable offline with a synthetic alignment.** That sentence is run
1's result and it stands. An offline Phase 2 control is not a cheaper version of the real
one; it is a measurement of the floor. Run 2 then showed that a real alignment is
*necessary and not sufficient*.

## Finding 4: a real internal gap, on real deposited coordinates

Unplanned, and worth having. **Barstar — 1BRS chain D, relabelled to chain A on the way
in — is missing 2 residues between 63 and 66**, so node 01's §8.12 fusion warning fired on
deposited coordinates rather than on a fixture, and node 02 disclosed it in
`prediction_inputs.json`:

> chain A is missing 2 residue(s) between 63 and 66, and the two sides are fused

The judge printed it in the provenance block. `C11`'s shape — a gapped target disclosed
through to the report — is therefore present in a real GPU run, not only in the Tier 1
tree. The disclosure is what §8.12 requires: the gap is not silently trimmed, and the
fusion is declared where a reader of the inputs will see it.

### A correction to this file's first version

The provenance row originally read "chain A (barnase, 87 res) on chain D (barstar, 108
res)". **The two lengths were swapped**: barnase is the 108-residue chain and barstar the
87-residue one, as the run's own `probe-1BRS.txt` says and as `--binder A --target D` in
`run_controls_phase2.sh` implies. Nothing downstream depended on it — no measured number
changes — but it is corrected here rather than silently, because a provenance table whose
numbers are decorative is worse than no provenance table.

The confusion has a cause worth naming: **there are two sets of chain letters.** The PDB's
(barnase A, barstar D) and the workflow's, because `make_controls.py` relabels the target
to A and the design to B. So "chain A" means barnase in the probe and barstar in
`prediction_inputs.json`. Anything quoting a chain letter in this workflow has to say
which set it is using.

## `C11` is MET, from run 2's data, with no further compute

`C11` asks for a real target carrying an internal gap, folded whole, with both chains'
pLDDT recorded and the whole-chain fold **not** degraded the way `8e4772a7`'s trim was. Run
2 satisfies every clause, because barstar turned out to carry a genuine gap. Nothing was
run for this; it is a claim read off a run that already happened.

| clause | measured |
|---|---|
| a real target with an internal gap | 1BRS barstar, **2 residues missing between 63 and 66**, ids 1..89 over 87 resolved residues, on deposited coordinates |
| node 01 detects it | `design_inventory.json` records `{"after": 63, "before": 66, "missing": 2}` |
| node 01 warns | the §8.12 warning fired, citing `8e4772a7` and the measured cost of trimming |
| folded whole, not trimmed | one 87-residue query; no segment was dropped |
| the fusion is disclosed | `prediction_inputs.json` → `target.fusion`, and carried into `target_alignment.json` |
| both chains' pLDDT recorded | target **85.42**, binder **38.55** |
| the whole-chain fold not degraded | **85.42** against `8e4772a7`'s trimmed ~42 |

The disclosure propagates through three files rather than being computed once and dropped,
which is what §8.12 asks for: the gap is not silently repaired, and a reader of the
prediction inputs is told the junction is a fusion.

**An incidental corroboration worth noting.** Node 01's warning quotes `8e4772a7`'s figure
for the undamaged case as "the whole chain reaches 81". This run measured **85.42** on a
different entry, a different organism and a different alignment. Those are not the same
experiment and the agreement is not a replication — but a number quoted from another
workflow's commit has now been found in the right neighbourhood by 040's own measurement,
which is more than it had before.

### What `C11` does not establish

- **It is not a controlled trim-versus-fuse comparison on this entry.** Only the fused arm
  was run. The controlled comparison is `8e4772a7`'s and already exists; re-deriving a
  settled decision is what §8.1 warns against, and the trimmed arm would also need the
  3600-row alignment sliced to the retained columns, which is a new opportunity to be
  wrong about something already known.
- **The gap is small.** Two residues, and the two flanks are close in space. A 20-residue
  disordered loop fused across is a different proposition, and nothing here speaks to it.
  `8e4772a7`'s case was the larger one, which is the reason to keep citing it.
- **85.42 is the gapped chain modelled well, not the gap modelled well.** The two residues
  are absent from the deposited coordinates and remain absent; pLDDT says the rest of the
  chain is confident, not that the junction is right. The disclosure exists precisely
  because that distinction cannot be read off the number.

## `C10` is MET, and it corrected its own rationale

`bash tests/run_c10.sh <phase2-run>`, judged by `tests/report_c10.py`. Two arms differing
by **exactly one byte**: the real `cplx_C8.fasta`, and the same file with its `:` deleted
by `tr -d` rather than retyped, so the arms cannot drift in any other character. Both go
to `colabfold_batch` **directly**.

**Bypassing node 03 is the right thing to do here and only here.** Node 03 refuses a
colon-less complex, so what the predictor does with one cannot be observed through the
workflow. The subject is the predictor's behaviour, not the workflow's, and the workflow is
built so this input never reaches it.

`alphafold2_multimer_v3`, `--msa-mode single_sequence`, no a3m, 1 model, 3 recycles. The
clamp §11.1 found matters only when an alignment is handed in and then discarded; nothing
is handed in here. The chain count does not depend on MSA depth.

| arm | input | chains returned | `iptm` | `ptm` | complained? |
|---|---|---|---|---|---|
| A | `target:binder`, 1 colon | **2** — 87 and 108 | 0.08 | 0.20 | no |
| B | colon deleted, 195 residues | **1** — 195 | **none** | 0.20 | **no** |

### The measurement

**The fusion is real and it is silent.** The two sequences become a single 195-residue
chain, a chimera that is not the complex anyone asked for, and `colabfold_batch` exits zero
and says nothing about it. A pipeline without node 03's guard would carry that structure
forward with no signal at all that the input had been degraded.

**And no `iptm` comes out.** ColabFold omits it for a single-chain prediction, which is the
sane behaviour: there is no interface to score. `ptm` is still reported at 0.20.

### Which refutes half of the rationale I wrote for this control

`run_c10.sh`'s header originally justified the control by saying the hazard was a fused
chain **still producing an `iptm`** — a plausible interface number from an input declaring
no interface, §8.6's trap in its sharpest form. The fusion half is correct. The `iptm` half
was **asserted, not measured**, and measurement refuted it.

That is §8.17's error again, this time in the rationale of the control sent to settle the
question — and the third rule earns its keep: the claim was not derived from anything, it
was carried over from the shape of a different finding (§8.6's "metrics computable on
garbage look ordinary", which is true of `sc_rmsd_A` and `plddt_binder` and happens not to
be true of `iptm` on a single chain).

**Node 03's own message never made that claim**, which is worth stating plainly: it says
only that the colon count is wrong and that node 02 verified the file on write, so the file
has changed since. That remains exactly right, and the guard remains necessary — for the
silent fusion, not for a fabricated interface number. Both new files are corrected in
place.

### What `C10` does not establish

- **One entry, one model, one seed.** Whether ColabFold's `iptm` omission is a stable
  contract or an artifact of this version is not established, and node 04 should keep
  nulling an absent metric rather than trusting the omission (§8.7).
- `iptm` 0.08 on arm A is *not* a finding about the colon. It is `single_sequence` with no
  alignment, consistent with every other unaligned number in this file. The control was
  about chain count.

## `C13`/`C14`: the designed-binder control, pre-registered before it is built

Phase 2's three failures so far were all one error about a different chain: a chain the
flags could not model, and no requirement that it be modelled. `C13`/`C14` are the attempt
to fix the control rather than the workflow, and **this section is written before the
scripts exist**, let alone before anything runs.

### The entry

`7XYQ`, confirmed by probing the file rather than from the paper or from memory
(sha256 `883c9c5542b2f65b`):

| chain | role | residues | ids | gaps | provenance |
|---|---|---|---|---|---|
| A | **target** | 203 | 1..207 | **165→170, 4 missing** | human PD-L1 (`CD274`), native, so it has homologs |
| B | **binder** | 124 | 1..124 | none | `DBL1_03`, deposited as a **synthetic construct** |

This is the first control in the repo where **both** chains can be modelled under the flags
§4 mandates: the target from an alignment, the binder from **nothing**, because a designed
sequence is supposed to fold unaided. That is the whole point, and it is exactly the move
`C12` made for Phase 1.

The target's 4-residue gap will be fused and disclosed. That is §8.12's path and `C11` just
measured it as benign at this scale, so it is a known cost rather than a confound.

### Pre-registered, with provenance for every figure

| | prediction | where it comes from | status |
|---|---|---|---|
| `plddt_target` | **> 70**, expect ~85 | barstar measured **85.42** with a real server MSA; §4's table says 83.7 | measured precedent |
| `plddt_binder_in_complex` | **> 70** | `C12` measured **95.05** for a designed sequence folded from sequence alone | **the new claim, and this control's entire premise** |
| **both** of the above | required before separation is reported at all | `8e4772a7` is **symmetric over the two chains** | **derived from the rule, not from the last failure** |
| `iptm`, `C13` vs `C14` | `C13` > `C14`, by a margin that is not noise | §8.6: the gate is separation, not an absolute | measured requirement |
| `iptm`, `C13` absolute | **not predicted** | 029's 0.860 was a native pair with deep MSAs on *both* chains; here the binder has none by design | explicitly unknown |
| target fusion | the 165→170 gap fused and disclosed through to the report | §8.12, `C11` | required behaviour |

The third row is the one that matters and the one I got wrong last time. `8e4772a7`
quantifies over **both** chains, so the condition does too. Writing the premise down was
never the missing piece; deriving it from the governing finding was (§8.17's third rule).

### Reading order, fixed in advance

1. **Are both chains above 70?** If not, `NOT INTERPRETABLE`, and the judge already
   enforces this and names the floored chain.
2. **If the BINDER is the one below 70**, that refutes this control's premise — that a
   deposited de novo binder folds unaided — and it is a finding about `DBL1_03` or about
   designed binders at 124 residues, **not** about 040. 124 residues is large for a
   minibinder and `C12` was 60 idealised ones, so this is a real possibility and it is
   named here so it cannot be reinterpreted later.
3. **Only if both clear 70**: does `iptm` separate? If yes, **Gate B Phase 2 is MET**. If
   no, that is the **first interpretable evidence against 040's interface scoring** this
   repo has ever had, and it would be a finding about the workflow.

### What would refute the standing diagnosis

The current diagnosis is that `C8`/`C9` failed because their binder was native and
unalignable. **That is wrong if both chains clear 70 here and `iptm` still fails to
separate.** Stated now, symmetric this time, with no escape hatch: if that happens, the
problem is 040's and not the control's.

### What would invalidate `C13`/`C14` themselves

- **If `DBL1_03` needs an alignment to fold.** The deposition says synthetic construct, but
  a binder built by grafting a motif onto a natural scaffold can carry a natural sequence,
  and a natural sequence needs an MSA like any other. Step 2 above is what detects this.
- **If the fused target models badly.** `C11` measured a 2-residue fusion as benign; this
  gap is 4 and in a different fold. If `plddt_target` comes in low, the fusion is the first
  suspect and `8e4772a7`'s larger case is the precedent to read.
- **`C14`'s scramble could fold well.** Phase 1's `C4` scored 95.81 on 115 alanines, so a
  well-modelled negative is expected and is not a defect — what matters is `iptm`, not the
  negative's pLDDT.

## Is the verdict worth anything?

Seven mutations, every one refused, and the honest tree checked first:

| mutation | refused with |
|---|---|
| positive `iptm` → 0.05 | `iptm` does not separate, margin −0.0500 |
| negative `iptm` ← positive's | `iptm` does not separate, margin +0.0000 |
| positive and negative swapped | `iptm` does not separate, margin −0.8000 |
| positive `iptm` removed | missing on a control and it is a gating metric (§8.7) |
| the negative's row deleted | `['C9']` produced no row in `design_scores.json` |
| the positive's row duplicated | design `'C8'` appears twice |
| `n3-refuse.log` deleted | the `single_sequence` refusal was never exercised (§8.13) |

### These seven were graded against a fabricated tree, and that is a real limitation

`falsify_controls_phase2.sh` checks the honest run passes *first* (§8.16), because if it
does not, all seven "refusals" are refusals of the pre-existing failure and prove nothing.
Run against the measured run above, it therefore **aborts**, correctly:

> ABORTING: the unmutated run does not pass the judge, so nothing below would mean
> anything. Fix that first (§8.16).

So the judge was graded against `--synthetic`: a minimal tree with typed, round,
deliberately unrealistic numbers (`iptm` 0.90 against 0.10) and a `_FABRICATED` key in the
file saying what it is. That mode was added *because* of this abort, and it earns its keep
— §8.13 says a guard that cannot be shown refusing is not a guard, and without it this
judge could not be shown refusing at all while Phase 2 fails.

**What it does not establish.** A judge shown refusing seven mutations of a fabricated tree
is a judge whose comparison logic works. It is not a judge shown refusing a *real* run that
should have failed, because no real passing run exists to mutate. Phase 1's judge has the
stronger record: its six mutations were applied to a run that genuinely passed. This one's
record is weaker, and it stays weaker until `C8`/`C9` run with a real alignment.

## The guard that was shown refusing

Before spending any GPU time, `run_controls_phase2.sh` runs an `n3-refuse` step requiring
node 03 to reject `PARAM_MSA_MODE_CPLX=single_sequence` on the complex path:

| | measured |
|---|---|
| names `single_sequence` in the refusal | yes |
| non-zero exit | yes, `exit=1` |

That guard exists because of the `max_seq` clamp (`colabfold/batch.py:137-142`): in
`single_sequence` mode ColabFold sets `max_seq = 1`, so a precomputed a3m is parsed and
then effectively discarded — while `config.json`, the re-saved a3m and `coverage.png` all
record that the alignment *was* used (`batch.py:284-287`, written from the parsed alignment
before the clamp). It is the §8.4 failure in library form: every artifact says the MSA was
used, and it was not. Node 03 refuses the combination rather than producing a run whose own
provenance lies about it.

The judge requires that log to exist. Deleting it is the seventh mutation above.

## What this does not say

- **It does not say 040's interface scoring is wrong.** It says this control cannot tell.
  The distinction is finding 2, and conflating the two would be the §8.2 error: a control
  that fails for a reason unrelated to the thing under test is a defect in the control.
- **It does not say 040's interface scoring is right either.** That is the point of
  finding 2's gap. Phase 2 currently has no positive control capable of returning a
  meaningful interface number, so the interface path is *plumbed and guarded* but
  **unvalidated as a score**.
- `C9` is a composition-preserving shuffle — an easy negative, as Phase 1's `C3` was. The
  hard negative is a design that looks right and does not bind, and nothing in this repo
  can produce one.
- Both arms used 5 models and a fixed seed handling, one run each. There is no variance
  estimate for the complex path at all; Phase 1's ~0.03 Å / ~0.04 pLDDT figure is for the
  monomer path and does not transfer.
- Nothing here is experimental validation. These are predictions of a deposited complex,
  and `1BRS` is in AF2's training data.
