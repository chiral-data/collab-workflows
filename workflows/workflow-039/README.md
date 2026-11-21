---
doc_id: workflow-039
domain: antibody-engineering
doc_type: workflow
version: "0.1.0"
deprecated: false
description: >
  Antibody developability triage. Takes a panel of VH/VL pairs or nanobody
  sequences, numbers and annotates them with ANARCI, models each Fv with
  ABodyBuilder2, computes TAP-style developability metrics, sequence
  liabilities and humanness, proposes AbLang2 CDR mutations, and ranks the
  panel into one HTML report. CPU only.
tags: [antibody-engineering, developability, antibody-numbering, structure-prediction, cpu]
---

# Workflow 039: Antibody Developability Triage

> **Status: implemented, never executed.** All seven nodes are written and their
> dependency-free logic is unit-tested (`python3 tests/test_logic.py`, 95 checks), but the image
> has not been built and no node has been run against the real tools — so no number in this
> workflow has been produced by running it. See
> [antibody_triage_pipeline.md](antibody_triage_pipeline.md) for the design rationale, the open
> items and the verification plan that is still outstanding.

This workflow answers the question a discovery team with a panel in hand actually asks: **which
of these clones do we express?** It takes a few hundred antibody sequences and returns a ranked
shortlist, with the risk flags that put each clone where it is.

Nothing else on `main` treats antibodies as antibodies. Workflows 016 and 017 dock them, 012,
018 and 026 predict their structures as generic proteins, and 013, 033 and 037 design binders de
novo — but none of them knows what a CDR is. This one does.

**Prior work worth knowing about:** the unmerged `workflow-abodybuilder3` branch (June 2026)
rebuilds workflow-010 as an ABodyBuilder3 antibody Fv structure predictor. It overlaps with node
02 of this workflow, which models the same thing with ABodyBuilder2 via ImmuneBuilder. If that
branch lands, node 02 is the piece to reconcile — this workflow would be a consumer of its
structures rather than a second modeller.

> **Developability flags are risk indicators, not predictions of clinical success.** The
> thresholds describe where clinical-stage antibodies sit. Several approved therapeutics violate
> them. A flag is a reason to look, not a reason to kill a clone. See
> [Outputs and interpretation](#outputs-and-interpretation).

## Overview

Six stages, in the order a discovery group would run them by hand:

1. **Annotate** — number every sequence, assign chain type and closest germline, delimit the
   CDRs, scan for sequence liability motifs.
2. **Model** — predict an Fv structure per clone with ABodyBuilder2 (or NanoBodyBuilder2 for
   VHH), keeping the model's own per-residue error estimate.
3. **Developability** — the five TAP-style metrics computed off that structure: total CDR
   length, patches of surface hydrophobicity, patches of positive and negative CDR charge, and
   asymmetry in net heavy/light chain surface charge.
4. **Humanness** — human 9-mer peptide content and germline identity, to flag clones likely to
   need humanization.
5. **Mutation scan** — AbLang2 proposes point mutants over the CDRs, and every proposal is
   re-screened through the liability motifs and developability rules so none trades one risk for
   another.
6. **Rank and report** — one table, one ranking, one self-contained HTML report.

## When to use this workflow

Use it on a **panel**: the output of a phage display campaign, an immunization, or a library
round, where you have tens to hundreds of sequences and limited expression slots. It is a
first-pass filter, and it is cheap enough to run on everything rather than on a hand-picked few.

**Do not use it when:**

- **You want to know what the antibody binds.** No antigen is considered anywhere in this
  workflow — no epitope mapping, no paratope analysis from a complex, no affinity prediction.
  Use workflow-012 or -018 for complex structure, -016 or -017 for docking.
- **You want immunogenicity prediction.** This workflow reports humanness, which is a *proxy*.
  It performs no T-cell epitope prediction — see [the MHC-II
  hole](antibody_triage_pipeline.md#the-mhc-ii-hole) for why not.
- **You have Fabs, full IgGs or scFvs.** Variable domains only; node 01 rejects anything else
  rather than guessing. Trim first, and split scFvs at the linker yourself.
- **You want TAP's own numbers.** Our metrics are reimplemented from the published guidelines and
  are labelled TAP-style throughout. They will not match the OPIG web server exactly.
- **You want aggregation or titre prediction.** Neither is included; see open items 5 and 6.

## Architecture and data flow

```text
panel CSV / FASTA (VH + VL per clone)
   │
   ▼
[01-annotate] ─────────► annotations.json, annotations.csv, liabilities.csv, clones.fasta
   ├───────────────────────────────┬──────────────────────────┐
   ▼                               ▼                          │
[02-model-fv] ──► models/*.pdb   [04-humanness] ──►            │
   │                               humanness.csv               │
   ▼                               │                           │
[03-developability] ──► developability.csv, flags.csv          │
   │                               │                           │
   ▼                               │                           │
[05-mutate-scan] ──► proposals.csv ◄───────────────────────────┘
   │
   ▼
[06-rank] ─────────────► ranked_panel.csv, panel_summary.json
   │
   ▼
[07-report] ───────────► report.html
```

Nodes 02 and 04 are independent and may run in parallel. **No node needs a GPU**, and no node
needs the network: all model weights and the humanness database are baked into the image at
build time.

The arrows above show where each measurement comes from; the `[dependencies]` table in
`.chiral/workflow.toml` is the authority. Nodes 06 and 07 also read the detail files directly —
06 joins node 02's model error and node 05's proposal counts, and 07 reads
`developability.json`, `humanness.json` and `scan_manifest.json` for the per-residue exposure,
the calibration block and the re-screen rules — so both depend on every upstream node that
produces something they display.

## Input requirements

A panel of variable-domain sequences, as **CSV** (`clone_id,vh,vl`) or **FASTA** (record ids
suffixed `_VH`/`_VL`). Nanobodies give a heavy chain only. Extra CSV columns — target, assay
EC50, library round — are carried through to the ranked output untouched.

See [input_files/README.md](input_files/README.md) for the exact format, the sequence
requirements, and the bundled two-clone smoke-test panel.

## Workflow nodes

### Node 01: Annotate Sequences

**Goal:** Every clone numbered, delimited and liability-scanned, or explicitly rejected.

**Process:** Numbers each chain with ANARCI under `numbering_scheme`, assigns chain type and
closest germline within `allowed_species`, delimits the six CDRs under `cdr_definition`, and
scans the sequence for liability motifs. Rejects anything that is not a variable domain.

**Scientific notes:** The CDR definition chosen here propagates into every downstream metric —
CDR total length, the charge and hydrophobicity patch metrics, and the mutation scan's editable
positions all depend on it. It is recorded in the manifest for that reason. The default liability
set is unpaired cysteines, N-glycosylation sequons, NG/NS deamidation, DG/DS isomerization,
Met/Trp oxidation and DP/NP acid cleavage motifs.

**Per-clone failures are recorded, not fatal:** one unnumberable sequence in a 500-clone panel
yields a row with a reason and the run continues.

**Outputs:** `annotations.json`, `annotations.csv`, `liabilities.csv`, `panel_manifest.json`,
`clones.fasta`

### Node 02: Model Fv Structures

**Goal:** One Fv structure per clone, with its error estimate attached.

**Process:** ABodyBuilder2 for paired VH/VL, NanoBodyBuilder2 for VHH, selected per clone when
`receptor_type = auto`. Weights come from the image, never the network. `refine_structures`
controls the OpenMM refinement pass. `--dry-run` resolves every setting without loading a model.

**Scientific notes:** ImmuneBuilder predicts an **ensemble** and writes the per-residue spread
into the B-factor column as an error estimate — so the uncertainty travels with the structure.
Upstream reports CDR-H3 RMSD of 2.81 Å on a 34-antibody benchmark, better than AlphaFold-Multimer
and over a hundred times faster. CDR-H3 is both the loop that matters most and the one modelled
worst; the report shows the error beside every structure-derived metric for that reason.

**Outputs:** `models/*.pdb`, `model_summary.csv`, `modelling_manifest.json`

### Node 03: Developability Metrics

**Goal:** The five TAP-style metrics, per clone.

**Process:** Computes solvent accessibility on each model (Shrake-Rupley), then total CDR
length, surface hydrophobicity patches, positive and negative CDR charge patches, and structural
Fv charge symmetry. Adds solvent-exposed liabilities (a buried NG is a smaller problem than an
exposed one) and an isoelectric point estimate. `tap_threshold_set` controls whether raw values
are flagged.

**Scientific notes:** This is a **reimplementation from the paper's Methods, not TAP itself** —
TAP has no redistributable source. The formulas are implemented as published: surface exposure
at ≥7.5% relative ASA by Shrake-Rupley, the CDR vicinity as exposed CDR and anchor residues plus
exposed residues within 4 Å, pair sums over residues closer than 7.5 Å, Kyte & Doolittle
normalized to [1,2], salt-bridged residues taking glycine's hydrophobicity and zero charge. Two
deviations remain — the exposure reference and the model generation — and both are recorded in
`developability.json` under `calibration`, along with the panel's PSH mean against the published
clinical-stage mean.

**Outputs:** `developability.csv`, `developability.json`, `flags.csv`

### Node 04: Humanness

**Goal:** How human each variable domain looks.

**Process:** promb's `human-oas` database — antibody 9-mers observed in at least 10% of human
subjects — gives the fraction of each domain's peptides that are human (an OASis-like score),
plus the nearest human peptide for each non-human one. Germline identity is computed against the
closest germline from node 01.

**Scientific notes:** **Humanness is not immunogenicity.** No MHC-II binding or T-cell epitope
prediction is performed, because no open, commercially usable class-II predictor exists. The
score is also germline-biased: an unusual-but-fine clone scores worse than a common-but-mediocre
one. This node needs no structure, so it still produces output for clones whose modelling failed.

**Outputs:** `humanness.csv`, `humanness.json`, `peptide_detail.csv`

### Node 05: Mutation Scan

**Goal:** A shortlist of point mutants worth testing.

**Process:** AbLang2 (paired model, VH and VL scored together) gives per-position amino acid
likelihoods over the region named by `scan_region`. Proposals are ranked by likelihood gain,
then **re-screened** through the liability motifs and the sequence-derived developability rules;
anything that introduces a new liability is dropped with a recorded reason. Capped at
`max_mutants_per_clone`. With `mutate_scan = false` the node writes an empty table and records
`enabled: false`, and nodes 06/07 omit the section.

**Scientific notes:** These are **likelihood-based, not affinity-based**. A high-likelihood
mutant is one that looks like a natural antibody, which is not the same as one that binds
better — and with no antigen in the pipeline, nothing here can speak to binding at all. Treat
every proposal as a bench hypothesis.

**Outputs:** `proposals.csv`, `scan_manifest.json`

### Node 06: Rank Panel

**Goal:** A decision, not a pile of numbers.

**Process:** Joins every per-clone metric into one table, applies the pass/warn/fail logic, and
sorts. Marks the top `top_n` as the shortlist. Carries through any extra input columns.

**Scientific notes:** **Nothing is filtered.** Every clone appears in the output, flags and all;
the ranking is a sort, and a user who disagrees with the weighting can re-sort the CSV.

**Outputs:** `ranked_panel.csv`, `panel_summary.json`

### Node 07: Report

**Goal:** One self-contained HTML report.

**Process:** Embeds all data inline and renders the ranked panel with per-metric flags, the
liability map drawn on each sequence, developability distributions across the panel with each
clone's position marked, the proposed mutants, and a 3Dmol viewer of the shortlist with CDRs and
flagged residues highlighted. Bootstrap 5, Plotly and 3Dmol from CDN; no local assets.

**Scientific notes:** The report leads with the interpretation caveats, shows model error beside
every structure-derived metric, and states which thresholds are published and which are
heuristics — reviewers in pharma discount a report that does not.

**Outputs:** `report.html`

## Parameters

### receptor_type

| Value / Range | Description |
|---------------|-------------|
| `auto` (default) | Decide per clone: paired VH/VL → ABodyBuilder2, heavy chain alone → NanoBodyBuilder2. |
| `antibody` | Force ABodyBuilder2; a clone without both chains is an error. |
| `nanobody` | Force NanoBodyBuilder2. |

**Trade-off:** None, except that `auto` can silently treat a clone whose VL failed to parse as a
nanobody. Set it explicitly for a homogeneous panel.

### numbering_scheme / cdr_definition

| Value / Range | Description |
|---------------|-------------|
| `numbering_scheme = imgt` (default) | Also `chothia`, `kabat`, `martin`, `aho`. |
| `cdr_definition = imgt` (default) | What TAP itself uses. Also `chothia`, `kabat`, `north`. |

**Trade-off:** The CDR definition changes CDR total length and every patch metric, so it changes
the flags. It is a convention, not a measurement — compare clones under one definition, never
across two. TAP numbers with IMGT and splits by IMGT CDR definitions, so the thresholds in node
03 were calibrated against IMGT boundaries; changing this makes the CDR-length flag in
particular incomparable. Each definition is evaluated in the numbering scheme it was defined on
(North on Chothia numbering, for instance), so node 01 numbers a second time when this differs
from `numbering_scheme`.

### refine_structures

| Value / Range | Description |
|---------------|-------------|
| `true` (default) | ImmuneBuilder's OpenMM refinement on every model. |
| `false` | Several times faster; leaves minor geometry violations in place. |

**Trade-off:** Refinement fixes local geometry that SASA-based metrics are sensitive to. Turn it
off for a fast look at a large panel, not for numbers you will act on. Upstream has no flag for
this — `save()` always refines — so the unrefined path writes the top-ranked model directly and
node 02 adds the error B-factors itself, keeping the two paths' files comparable.

**Test vs production:** `false` for a smoke test, `true` for a real run.

### tap_threshold_set

| Value / Range | Description |
|---------------|-------------|
| `clinical` (default) | The amber/red bounds of Table 2 of Raybould et al. (2019), derived from 242 clinical-stage therapeutics. |
| `none` | Report raw metric values with no flags. |

**Trade-off:** The bounds are published and are implemented with their citation, but they were
calibrated on TAP's own pipeline — side-chain exposure against a reference we do not have, and
ABodyBuilder v1 models. Our values therefore carry an offset of unknown size, which is why node
03 also reports the panel's PSH mean against the published clinical-stage mean of 123.30 ± 16.60
so a gross mismatch is visible. Set `none` if you would rather read raw numbers than flags with
a caveat. See ["the calibration
problem"](antibody_triage_pipeline.md#the-calibration-problem).

### liability_motifs

| Value / Range | Description |
|---------------|-------------|
| `default` | Unpaired Cys, N-glycosylation sequons, NG/NS deamidation, DG/DS isomerization, Met/Trp oxidation, DP/NP acid cleavage. |
| path to JSON | Your own motif set, which is what most groups will want. |

**Trade-off:** Motif scanning is pattern matching: it over-calls (many sequons are never
glycosylated) and the structural context in node 03 is what separates a real risk from a hit.

### humanness_db / peptide_length

| Value / Range | Description |
|---------------|-------------|
| `human-oas` (default) | Antibody 9-mers seen in ≥10% of human subjects; the OASis-like score. |
| `human-reference` | The human reference proteome (UP000005640). |
| `human-swissprot` | Human SwissProt. |
| `peptide_length = 9` (default) | Ignored for `human-oas`, whose peptide length is fixed at 9. |

**Trade-off:** `human-oas` is the right reference for antibodies — it is built from observed
repertoires rather than the proteome at large. The proteome databases answer a different
question and will score antibodies lower.

### mutate_scan / scan_region / max_mutants_per_clone

| Value / Range | Description |
|---------------|-------------|
| `mutate_scan = true` (default) | Run the scan. |
| `scan_region = cdr` (default) | All six CDRs. Also `cdrh3` (cheapest, highest-value) or `all`. |
| `scan_stepwise_masking = false` (default) | Score positions from a single forward pass per clone. `true` masks each scanned position in turn — more standard for asking a masked language model what belongs at a position, and one forward pass per position instead of per clone. |
| `max_mutants_per_clone = 20` (default) | Cap after re-screening. |

**Trade-off:** Cost scales with clones × positions × 19 candidate residues, and with
`scan_stepwise_masking` also with a forward pass per position. `all` includes framework
positions, where mutations are more likely to hurt expression and humanness than help binding.

**Test vs production:** `cdrh3` to exercise the node, `cdr` for a real run.

### top_n / seed

`top_n = 10` clones in the report's detail section and 3D viewer. `seed = 42` is recorded for
provenance; only the mutation scan's tie-breaking among equally scored proposals is stochastic.

## Outputs and interpretation

### ranked_panel.csv

The deliverable: one row per clone, every metric, every flag, and the rank. **Every clone is
present** — nothing is filtered at any stage. The top `top_n` are marked as the shortlist.

Read the flags as risk indicators. The TAP thresholds describe where clinical-stage antibodies
sit, and several approved therapeutics fall outside them; a red flag means "understand this
before you commit", not "discard".

### developability.csv / developability.json

Per-clone metric values, plus the implementation choices they depend on (SASA algorithm, probe
radius, hydrophobicity scale, CDR definition) recorded in the JSON. **These are TAP-style
values** and will not match the OPIG web server exactly — the dependence on those unpinned
choices is exactly why.

### liabilities.csv

One row per motif hit: clone, chain, position, motif, and whether node 03 found it
solvent-exposed. An exposed sequon in CDR-H3 and a buried one in framework 2 are not the same
risk, and this is where that distinction lives.

### humanness.csv / peptide_detail.csv

The human peptide content per domain, germline identity, and per-peptide detail with the nearest
human peptide for each non-human 9-mer — which is what tells a humanization campaign where to
start. **Not an immunogenicity prediction.**

### proposals.csv

Proposed point mutants with their likelihood gain and their re-screening verdict, including
dropped proposals and why. Likelihood-based, not affinity-based; no antigen was considered.

### models/*.pdb

The Fv models. **The B-factor column holds ImmuneBuilder's per-residue error estimate**, not a
crystallographic B-factor — high values mean the model is unsure there, and CDR-H3 is usually the
highest.

### report.html

Self-contained; opens in any browser with no local assets.

## Quick start

### Running with Docker

Built from [Dockerfile](Dockerfile) as `antibody-triage:2026_11_11`.

```bash
docker build -t antibody-triage:2026_11_11 workflows/workflow-039/
```

Smoke-test the image before running any node:

```bash
docker run --rm antibody-triage:2026_11_11 python -c "import ImmuneBuilder, ablang2, promb, anarci; print('imports ok')"
```

```bash
docker run --rm antibody-triage:2026_11_11 ABodyBuilder2 --help
```

Then run a node, node 02 shown:

```bash
docker run --rm -v /tmp/wf039/inputs:/workspace/inputs -v /tmp/wf039/outputs:/workspace/outputs -v $(pwd)/workflows/workflow-039/02-model-fv:/workspace -w /workspace -e PARAM_RECEPTOR_TYPE=auto -e PARAM_REFINE_STRUCTURES=true antibody-triage:2026_11_11 bash run.sh --dry-run
```

Drop `--dry-run` to model for real. Nodes 02 and 05 both support it.

### Unit tests

The dependency-free logic — CDR and anchor derivation, liability motif scanning, the
hydrophobicity normalization, the TAP pair sums on hand-built geometry, salt-bridge detection,
every Table 2 flag boundary, the mutation re-screen rules and the ranking order — is covered by
tests that need nothing installed:

```bash
python3 workflows/workflow-039/tests/test_logic.py
```

95 checks, no third-party imports. They catch the arithmetic and the boundary conditions, which
is where a reimplementation of someone else's published metric goes wrong quietly. They are not
a substitute for running the pipeline.

### Running on Silva

1. Select this workflow from the workflow list
2. Upload your panel CSV or FASTA
3. Adjust parameters if needed (see Parameters)
4. Click Run

### Test vs production settings

| Setting | Test | Production |
|---------|------|------------|
| input | bundled two-clone panel | your panel |
| `refine_structures` | `false` | `true` |
| `scan_region` | `cdrh3` | `cdr` |
| `tap_threshold_set` | `clinical` | `clinical` |
| `scan_stepwise_masking` | `false` | `true` if the proposals matter |

## Requirements

- **CPU only.** No GPU node. This is the point of the workflow: a few hundred clones is a CPU
  job, which is what makes triaging a whole panel affordable rather than triaging a hand-picked
  ten. Wall-clock figures are deliberately absent until they are measured (verification step 9).
- **Network:** none at run time. ImmuneBuilder, AbLang2 and promb all fetch on first use
  upstream, so all three are baked into the image at build time.
- **Licensing:** every component is BSD-3, MIT or Apache-2.0. No academic-use-only tool, and no
  22 GB database download — see the pipeline doc on
  [the MHC-II hole](antibody_triage_pipeline.md#the-mhc-ii-hole) and on BioPhi's OASis.
- **Image:** conda base (OpenMM, pdbfixer and ANARCI are conda-first upstream), Python 3.11,
  Torch CPU wheels.

## Troubleshooting

**Node 01 rejects my sequence as not a variable domain**
ANARCI did not recognise it as the chain type claimed, or ImmuneBuilder's `len > 70` assertion
failed. Fabs, full IgGs, scFvs and CDR fragments all fail here by design. Trim to the variable
domain; split scFvs at the linker yourself.

**A clone was modelled as a nanobody when it is an antibody**
`receptor_type = auto` falls back to NanoBodyBuilder2 when it sees only a heavy chain, so a VL
that failed to parse looks like a VHH. Check `panel_manifest.json` for the per-clone chain call
and set `receptor_type = antibody` explicitly.

**Node 02 is slow**
OpenMM refinement dominates. Set `refine_structures = false` for a first pass, accepting that
SASA-derived metrics get noisier.

**Node 03's numbers disagree with the TAP web server**
Expected. Ours is a reimplementation from the published guidelines, and the paper does not pin
every implementation choice. Compare flag direction, not values, and check the recorded choices
in `developability.json`.

**Node 03 warns that the panel's PSH mean is far from the published CST mean**
The published clinical-stage mean is 123.30 ± 16.60. More than 3 SD away is either a genuinely
unusual panel or — more likely — the calibration gap: the exposure reference or the
hydrophobicity normalization. Treat the PSH flags as unreliable until it is explained, and set
`tap_threshold_set = none` to read raw values in the meantime.

**Node 03 fails a clone with "numbering mismatch with node 01"**
More than 20% of the model's residues could not be matched to node 01's numbering, so the CDR
assignment would be wrong. It means node 01 and node 02 numbered with different schemes — check
that `numbering_scheme` was not changed between them.

**Node 04 scores a known human antibody as poorly human**
Check `humanness_db` — the proteome databases (`human-reference`, `human-swissprot`) answer a
different question than `human-oas` and will score antibodies lower.

**Every proposal in node 05 was dropped**
The re-screening is working: the high-probability mutations at those positions introduce
liabilities or deepen a risk the clone is already flagged for. Widen `scan_region` and read the
recorded reasons in `proposals.csv` — each dropped row carries the rule that rejected it.

**The image fails to build at the conda step**
Known risk, open item 1 in the pipeline doc: bioconda's `anarci` may lag Python 3.11. Pin
`python=3.10` and confirm promb still installs.

## References

- [Raybould et al., *Five computational developability guidelines for therapeutic antibody
  profiling*, PNAS 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6410772/) — the TAP
  guidelines this workflow's node 03 reimplements. CC BY.
- [brennanaba/ImmuneBuilder](https://github.com/brennanaba/ImmuneBuilder) — ABodyBuilder2 and
  NanoBodyBuilder2, BSD-3, PyPI 1.2.
- [oxpig/ANARCI](https://github.com/oxpig/ANARCI) — antibody numbering, BSD-3.
- [oxpig/AbLang2](https://github.com/oxpig/AbLang2) — paired antibody language model, BSD-3,
  PyPI 0.2.1.
- [MSDLLCPapers/promb](https://github.com/MSDLLCPapers/promb) — human peptide prevalence, MIT,
  PyPI 1.0.2.
- [antibody_triage_pipeline.md](antibody_triage_pipeline.md) — this workflow's design rationale,
  decisions, open items and verification plan.
