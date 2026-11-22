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

## Finding 3: what has to happen next, and why it needs the user

**`C8` and `C9` need a real target MSA, which needs one fetch from `api.colabfold.com`.**

That is an outward-facing network call from the GPU host, so `run_controls_phase2.sh` will
not make it without `--allow-egress`, which switches `PARAM_MSA_SOURCE` to `server` and
opens `--network=bridge` for that one node only. The flag exists so that reaching the
network is a decision somebody made, not a thing a test script did.

This is the user's call and has not been made. It is the same fetch `C7` was already
waiting on, so `C7`, `C8` and `C9` now unblock together.

**Separation is not measurable offline with a synthetic alignment.** That sentence is the
result of this run. An offline Phase 2 control is not a cheaper version of the real one; it
is a measurement of the floor.

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
