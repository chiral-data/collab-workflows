# Tier 2, Phase 1 controls — measured

Real ColabFold predictions of five controls through the whole workflow, on the GPU, with
no network. The numbers below are what came out, including the ones that say the control
set was wrong.

Reproduced by `bash tests/run_controls.sh`. Judged by `tests/report_controls.py`, which
prints this table and the separation margin. **`bash tests/falsify_controls.sh <run>`
shows that judge refusing**, which is the only reason to believe it when it agrees.

**Verdict: Gate B Phase 1 MET** — `C12` passes, `C3` and `C4` fail, and `sc_rmsd_A`
separates them by 17.90 Å. `plddt_binder` does not separate, by −0.76, and that is a
recorded finding rather than a failure; see finding 3. `C6`, the killed prediction, is
also MET. Only `C5` is outstanding, and it is blocked on inputs.

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
| Nodes run | 037 node 04, then 040 nodes 01, 02, 03, 04, 05 — **not 01b**, which did not exist yet |

**Both runs predate node `01b-align-target`, and the numbers are unaffected.** With
`score_interface` false, node 02 never calls `load_target` and never reads anything 01b
writes, so 01b contributes nothing to any value in this file. `run_controls.sh` now runs
it anyway, so the script matches the DAG and the no-op branch gets exercised on the GPU
host — but that is DAG fidelity, not a correction, and re-running these controls would
produce the same table for the same reason. Said here rather than left to be worked out,
because "the evidence was produced by an earlier version of the script" is exactly the
sort of thing that quietly stops being true.

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

## C6 — a real prediction killed at a short budget

The last Phase 1 control, and the one that tests node 03's `kill_tree()` against a real
`colabfold_batch` holding GPU memory rather than against a shell stand-in. Run by
`bash tests/run_c6.sh <control-run> <out> 35`, judged by `tests/report_c6.py`.

**Verdict: MET**, all five claims, 2026-11-21.

| claim | measured |
|---|---|
| stopped at the budget, did not complete | 177 s over 5 designs at a 35 s budget; a completed prediction takes 45–67 s **each** |
| the predictor was really running | **85** `colabfold` process lines sampled inside the container; GPU peaked at **2357 MiB** against a 1 MiB baseline |
| nothing of the predictor survived | a real process table, no predictor process, in a container that **outlived the node** |
| the GPU was freed | 1 MiB before, 1 MiB after |
| every design recorded failed, with the timeout and budget named | 5 of 5, each reason naming the 35 s budget |

### Why the container has to outlive the node

The obvious version of this control is `docker run` a node that times out. It is worthless:
a container's exit kills every process in it and frees the device no matter how badly
`kill_tree()` is written, so the check would pass for the wrong reason (§8.4). The same
problem is why the Tier 1 suite reports `GAP` rather than `PASS` for the worker-survival
claim in the image — a pid recorded inside a container means nothing on the host.

So `run_c6.sh` starts the container first with `sleep infinity`, runs node 03 through
`docker exec`, and *then* asks what is still alive inside it. And because "nothing
survived" is vacuous if nothing was ever there, a sampler records what was running
mid-prediction. That sampler is why claim 2 exists and why it is checked before claim 3.

### Two defects in this control, before it was allowed to pass

**The first version passed claim 3 for exactly the wrong reason.** It sampled with
`docker exec ps`, and the image has no `ps`, so the after-listing contained
`OCI runtime exec failed: ... "ps": executable file not found` — no process lines, which
the judge read as "nothing survived" and passed. §8.4 again, and this time inside the
checking code. The fix is `docker top`, which runs the *host's* `ps` against the
container's pids and prints a `PID` header the judge can require as proof it really
listed. `report_c6.py` now refuses a listing it cannot show to be a listing: a listing
that could not be taken is not an empty listing.

**And at a 20 s budget the recorded reason was empty**, which turned out to be a real
defect in node 03 rather than a cosmetic one. See below.

## The defect C6 found in node 03

`colabfold_batch`'s stdout is a pipe, not a tty, so it is block-buffered. A process killed
before the buffer fills takes every line with it:

| | `<design>_predict.log` (node 03's stdout capture) | `<design>_log.txt` (colabfold's own) |
|---|---|---|
| completed prediction | 3042 bytes | written |
| killed at 35 s | **0 bytes** | **written, 5 progress lines** |

So node 03's timeout reason read `Last lines: ` with nothing after it. The one thing a
reader needs at that moment was the one thing lost — and it was not actually lost, because
ColabFold logs through the `logging` module to its own file, which node 03's `collect()`
has *already copied into `./outputs`* before the reason is built.

Fixed: `predictor_tail()` falls back to that file and **names it as the source**, because
"the predictor said nothing" and "the predictor's output was buffered away" are different
facts and a reason that conflates them is worse than one that admits the gap (§8.7). The
same blind spot applied to the OOM scan, which only ever saw stdout — an allocator OOM
followed by a timeout was invisible, since the allocator warns, the job keeps the GPU
without progressing, the budget kills it, and every line saying so was in the lost buffer.
The fallback file is scanned for the markers too.

**Why Tier 1 could not have found this.** `STANDIN_SIGTERM` prints to stdout and flushes
line by line, so `timeout_message_carries_output` passed on output the real predictor does
not produce at that moment. The check was right about the stand-in and wrong about the
world. There is now a second stand-in, `STANDIN_SILENT_LOGFILE`, which writes only to its
log file and says nothing on stdout, and a case
`timeout-tail-falls-back-to-predictor-log` that requires the fallback, requires the source
to be named, and requires node 03's stdout capture to be genuinely zero bytes — otherwise
the case would be testing the ordinary path under a different name.

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
