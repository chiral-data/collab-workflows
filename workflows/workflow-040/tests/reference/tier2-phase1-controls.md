# Tier 2, Phase 1 controls — measured

Real ColabFold predictions of five controls through the whole workflow, on the GPU, with
no network. The numbers below are what came out, including the ones that say the control
set was wrong.

Reproduced by `bash tests/run_controls.sh`. Judged by `tests/report_controls.py`, which
prints this table and the separation margin. **`bash tests/falsify_controls.sh <run>`
shows that judge refusing**, which is the only reason to believe it when it agrees.

**Verdict: Gate B Phase 1 MET** — `C12` passes, `C3` and `C4` fail, and `sc_rmsd_A`
separates them by 17.90 Å. `plddt_binder` does not separate, by −0.76, and that is a
recorded finding rather than a failure; see finding 3.

## Provenance

| | |
|---|---|
| Host | `amd-c12-1`, NVIDIA GeForce RTX 3060 Ti 8 GB |
| Image | `colabfold:2026_10_08` (built from `workflow-026/Dockerfile`) |
| Settings | `--model-type alphafold2_ptm --msa-mode single_sequence --num-models 5 --num-recycle 3` |
| Network | `docker run --network none`, `host_url_passed=false`, every prediction offline |
| GPU | 1 MiB used before, 1 MiB after — the device was left free |
| Wall clock | 45–67 s per design |
| Runs | two independent runs, 2026-11-21 (`wf040-controls-run`, `...-run2`) |
| Entries | `4ZQK` sha256 `f5012ba58d76099b…`, `1BRS` sha256 `1c98dcc3eac7ccb8…`, fetched from RCSB and re-fetched independently on the host to the same digests |

Chains were chosen from a probe of each entry, not from memory — `tests/README.md`,
"Choosing a fixture". That mattered: **4ZQK chain B carries an 8-residue internal gap**
(52→61), so it would have been a positive control that failed for a reason unrelated to
the score. Chain A, 115 residues, no gaps, is the binder. In 1BRS, chain A is barnase and
**starts at residue 3**, which is the renumbering hazard (`ef74646e`) in the control set
for real.

## What came out

| | role | len | `sc_rmsd_A` | `plddt_binder` | `ptm` | passes |
|---|---|---|---|---|---|---|
| `C12` | positive, **designed** | 60 | **0.9272** | **95.05** | 0.86 | **yes** |
| `C2` | positive, native (barnase) | 108 | 11.5737 | 43.09 | 0.36 | no |
| `C1` | positive, native (PD-1) | 115 | 16.5096 | 49.33 | 0.40 | no |
| `C3` | negative, scrambled | 115 | 18.8308 | 44.94 | 0.30 | no |
| `C4` | negative, poly-alanine | 115 | 49.0800 | **95.81** | 0.53 | no |

Thresholds: `max_sc_rmsd_a = 2.0`, `min_plddt_binder = 80.0`, both uncalibrated.

Separation, worst positive against best negative:

| metric | margin | worst positive | best negative |
|---|---|---|---|
| `sc_rmsd_A` | +2.3212 | `C1` 16.5096 | `C3` 18.8308 |
| `plddt_binder` | **−52.72** | `C2` 43.09 | `C4` 95.81 |

## Reproducibility

Two independent runs, same settings, same seed handling:

| | run 1 | run 2 | Δ |
|---|---|---|---|
| `C1` `sc_rmsd_A` | 16.5384 | 16.5096 | 0.029 |
| `C2` `sc_rmsd_A` | 11.6083 | 11.5737 | 0.035 |
| `C3` `sc_rmsd_A` | 18.8214 | 18.8308 | 0.009 |
| `C4` `sc_rmsd_A` | 49.0691 | 49.0800 | 0.011 |
| `C1` `plddt_binder` | 49.29 | 49.33 | 0.04 |

Run-to-run variation is ~0.03 Å and ~0.04 pLDDT. That is the first variance estimate
either workflow has had, and it is two runs, not a distribution — it says the pipeline is
deterministic enough that a 2 Å threshold is not being decided by noise. It says nothing
about variance across seeds or models, which `tests/README.md` still lists as unestablished.

## Three findings

### 1. The pipeline is right, measured end to end on the input class it is for

`C12` is a textbook alpha helix (CA radius 2.3 Å, rise 1.5 Å, 100°/residue, cross-checked
against 2ZTA) carrying a designed `ALEKIAELKA` sequence. Refolded from sequence alone,
with no alignment, it recovers its own backbone to **0.93 Å at 95.05 pLDDT** and is the
one design exported. Ingest, input construction, a real offline GPU prediction,
superposition, filtering and the report all produce the right answer on a designed,
idealised binder — which is what this workflow exists to score.

### 2. The native positive controls are the wrong controls, and §4 said so in advance

`C1` 49.33 pLDDT, `C2` 43.09. Neither native protein folds from its sequence alone.

This is not a surprise so much as an unheeded prediction: §4 justifies
`msa_mode = single_sequence` on the grounds that *a de novo binder has no homologs, and
designed proteins are idealised*. A native immunoglobulin domain is the opposite — it
depends on the alignment it is not being given. So a native positive control run at
`single_sequence` measures AF2's single-sequence weakness on natural proteins, not 040's
scoring, and `C1`/`C2` cannot do the job §9 assigned them.

`C12` is what separates the two explanations. Without it, "the positives failed" is
equally consistent with a broken scorer; with it, the scorer is demonstrably fine and the
control set was mis-specified. Neither answer was available from C1–C4 alone.

`C1` and `C2` are kept, as **diagnostics** rather than positives: they are now the measured
record of what `single_sequence` does to a native chain, which is a number this repo did
not have and which bears directly on the Phase 2 MSA decision in §4.

### 3. `plddt_binder` is unsound as a gate, and the poly-alanine control proves it

**`C4` is 115 alanines and scores 95.81 pLDDT — the highest in the set, a hair above the
real design's 95.05, and comfortably past the 80.0 threshold.**

It passes `min_plddt_binder` while sitting 49 Å from the backbone it is supposed to
reproduce. Only `sc_rmsd_A` catches it, and it catches it by a mile.

The reason is not a bug: poly-alanine makes an excellent helix, and AF2 is correctly
confident that it does. Confidence in *a* fold is not evidence of *the designed* fold. §8.6
says metrics computable on garbage look ordinary; this is the sharper version — on garbage
this metric looks **excellent**, and it is a gating metric.

So `plddt_binder` is necessary-and-nowhere-near-sufficient. It can reject a design the
predictor is unsure about; it cannot admit one. `sc_rmsd_A` is the metric that
discriminates, and the `plddt_binder` ≥ 80 threshold must be read knowing that a sequence
encoding nothing clears it at 95.8.

## Is the verdict worth anything?

A judge that only ever says MET is a judge that is not reading. `falsify_controls.sh`
mutates one number or flag in a copy of this run and requires a NOT MET and a non-zero
exit each time. All six are refused, and the honest run still passes:

| mutation | refused with |
|---|---|
| `C12` forced to fail | wrong verdict for `['C12']` |
| `C4` forced to pass | wrong verdict for `['C4']` |
| `C3` `sc_rmsd_A` → 0.5 | no separation on `['sc_rmsd_A']` |
| `C12` `sc_rmsd_A` → 99.0 | no separation on `['sc_rmsd_A']` |
| `C12` `sc_rmsd_A` → `null` | no separation on `['sc_rmsd_A']` |
| `C12_provenance.json` deleted | gating control ran without recorded provenance |

The script checks the unmutated run passes first, because if it did not, all six
"refusals" would be the pre-existing failure and would prove nothing (§8.16).

### Two defects in the judge, both found this way

**The verdict could be lost without anything looking wrong.** `report_controls.py`
indexed its provenance fields directly, and `C12`'s record — the one written by hand
rather than measured by `make_controls.py` — did not carry them. The `KeyError` landed
*after* the table and *before* the verdict, so the second control run's `VERDICT.txt` is a
file that reads like a complete report and contains **no verdict at all**. The numbers in
it were right; the conclusion was simply absent, and the only signal was an exit code. §8.3
in its purest form.

Both ends are fixed. `make_c12_provenance.py` now measures `C12`'s record off the
generated file with the same reader and the same field names as the other four, and the
printer reports an absent field instead of dying on it.

**And an absent record passed while a partial one refused.** The first version of the
provenance gate iterated over the records that *existed*, so deleting
`C12_provenance.json` outright produced no row, no note and a clean MET — while a
*partially* filled record was correctly refused. Exactly inverted: the worse input was the
one that got through. The controls expected to have provenance are now enumerated, so an
absence produces a row. This is the sixth line of the table above, and it is in the table
because it was not caught by reading the code.

## What this does not say

- Nothing here is experimental. `C12` recovering its backbone says its sequence encodes its
  structure; it says nothing about binding, and `C12` has no interface at all.
- `C12` is synthetic — an idealised helix with a repeating designed sequence. It is the
  easiest possible designed positive. It establishes that the scorer works; it does not
  establish how the scorer ranks real 037 designs. That is `C5`, which is still blocked on
  inputs.
- The negatives are a scramble and a poly-alanine. Both are easy negatives. The hard case
  is a design that looks fine and does not bind, and nothing in this repo can produce one.
- `C3`, the composition-preserving scramble, lands 18.83 Å against the native's 16.51 Å.
  Those 2.3 Å are not a usable separation; both are simply "did not refold". When the
  positive does not refold either, RMSD cannot discriminate — finding 3 of §9's margin is
  downstream of finding 2, and the margin becomes meaningful only once the positives are
  designs.
