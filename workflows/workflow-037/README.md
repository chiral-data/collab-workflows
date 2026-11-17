---
doc_id: workflow-037
domain: protein-design
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  Generative de novo protein design with RFD4-Proteina. Takes a target structure
  (uploaded or fetched from RCSB) plus guided parameters, renders an RFD4 condition
  spec, and generates binders, motif scaffolds, monomers or redesigned binders.
tags: [protein-design, rfdiffusion, generative-model, binder-design, gpu]
---

# Workflow 037: RFD4-Proteina De Novo Design

RFD4-Proteina (RosettaCommons + NVIDIA GenAIR, released 2026-10-23, BSD-3) is a generative protein
design model that supersedes RFdiffusion 1/2/3 with a single model covering monomers, binders, motif
scaffolding, antibodies and enzymes. This workflow takes a target structure and a small set of
guided parameters, renders an RFD4 condition spec, validates it on CPU, and generates designs on one
GPU. It answers "give me candidate backbones that bind or scaffold this thing" — it does not score
them.

## Overview

RFD4-Proteina generates protein structure and sequence together, conditioned on a *condition spec*:
a small YAML file naming an input structure, a `contig` describing the chain layout to build, and a
set of `conditions` that pin parts of the input. The conditioning vocabulary is what makes one model
cover many tasks — `C_CRD` holds coordinates, `C_CRDGRP` holds a set of residues as one rigid group,
`C_SEQ` holds residue identity, and `C_HOT` marks hotspots the design should engage. Each condition
takes selections written as `sel('chain/res_name/res_id/atom_name')`, with `*` as the wildcard and
`[A,B]` for lists, so conditioning can be as coarse as a whole chain (`sel('A')`) or as fine as
named atoms (`sel('A/*/*/[N,CA,C,O]')`).

This workflow wraps that interface rather than exposing it raw: parameters describe the design task
in ordinary terms (target chain, binder length, hotspots) and node 02 compiles them into a spec. The
output is a set of mmCIF structures, each containing a generated chain alongside the target, plus a
geometric metrics table. **There is no oracle here.** Generation runs with `~metrics`, so nothing in
the output estimates whether a design folds as predicted or binds at all — that requires refolding,
which is deliberately out of scope (see [When to use this workflow](#when-to-use-this-workflow)).

## When to use this workflow

Use it when you have a target structure — uploaded, or just an RCSB accession — and you want
candidate backbones against it: a new binder, a scaffold around a functional motif, an unconditioned
monomer, or a regenerated version of an existing binder. You need a GPU from the supported set
(A100, L40S, RTX 6000 Ada, H100, H200, Blackwell). A few designs take minutes; a real campaign means
raising `num_replicates` into the hundreds and filtering afterwards.

Do not use it as a one-stop binder campaign: this workflow generates and stops. It reports no ipTM,
no pLDDT, no self-consistency RMSD, and the metrics it does report are geometry and composition, not
confidence. For an end-to-end scored campaign use **workflow-033** (RFdiffusion → ProteinMPNN →
Boltz2) or **workflow-013** (BoltzGen); to fold a sequence you already have, use **workflow-026**
(ColabFold). Feeding designs from here into workflow-026 returns pLDDT and PAE for each design. Be
clear about what that is: a confidence estimate, **not** self-consistency. Nothing compares the
prediction against the backbone generated here, because the FASTA carries the sequence and the
designed backbone stays behind. Self-consistency RMSD is the standard test for a generative design
and neither workflow performs it today, so a sequence that refolds confidently into a different fold
passes. A step that did perform it would read node 03's `*_generated.cif` and node 04's
`design_metrics.json`, not the exported FASTA, because those are what carry the backbone and the
identity of the designed chain. Workflow-026 also filters nothing — it scores, and the
judgement is yours. The most common mistake is
residue numbering — see [Input requirements](#input-requirements).

## Architecture and data flow

```text
target (upload or PDB ID)
   │
   ▼
[01-acquire-target] ──► target.cif, target_inventory.json
   │
   ▼
[02-build-spec] ─────► design_spec.yaml, spec_preview.cif, target.cif
   │
   ▼
[03-generate]  (GPU) ─► *_generated.cif, generation_manifest.json
   │
   ├──────────────────► [04-analyze] ──► design_metrics.csv/.json
   │                         │
   ├─────────────────────────┼────────► [05-report] ──────► report.html
   │                         │
   └─────────────────────────┴────────► [06-export-fasta] ─► sequences.fasta
```

Nodes run 01 → 02 → 03 → 04, then 05 and 06 both read node 03's structures and node 04's metrics.
Only node 03 requires a GPU. Node 06's `sequences.fasta` is named for where it goes next: drop it
into workflow-026 as its input to refold the designs. Refolding is not validation, though — scoring
a design against its own backbone needs the backbone, so that step reads `*_generated.cif` and
`design_metrics.json` instead. See [Node 06](#node-06-export-fasta).

## Input requirements

One target structure, supplied either way:

- **Upload** a `.cif`, `.pdb` or `.cif.zst` file into the workflow's input directory, or
- **Set `pdb_id`** to an RCSB accession and let node 01 fetch it.

An upload always wins over `pdb_id`. See [input_files/README.md](input_files/README.md) for the
PD-L1 (5O45) sample target used by the correctness gate.

**Residue numbering matters.** `contig` and `conditions` in an RFD4 spec refer to the original
mmCIF `label_seq_id`, not author numbering. RCSB files are author-numbered, so a hotspot copied
from a paper is almost always author numbering — leave `residue_numbering` at `auth` and node 02
translates it. On 5O45 the two differ by 16: author residue 45 is `label_seq_id` 29. Getting this
wrong does not error, it conditions the wrong residues.

## Workflow nodes

### Node 01: Acquire Target

**Goal:** Produce one normalised mmCIF target plus an inventory of its chains and both residue
numberings.

**Process:** Picks the first `.cif.zst`/`.cif`/`.pdb` in the input directory, else fetches
`https://files.rcsb.org/download/<pdb_id>.cif`. A `.cif.zst` is decompressed and a `.pdb` is
converted with Biopython; an mmCIF is copied byte-for-byte, because re-serialising it would renumber
`label_seq_id` and silently invalidate every contig downstream. A converted PDB then has its label
numbering canonicalised — one `label_asym_id` per author chain, `label_seq_id` sequential over
standard residues — because Biopython's writer starts a new `label_asym_id` at every modified
residue, which would otherwise shred a peptide chain into one- and two-residue fragments. The inventory is then read straight
out of `_atom_site`: per polymer chain, the `label_asym_id`/`auth_asym_id` pair, both sequence
ranges, the one-letter sequence, and explicit `auth_to_label` / `label_to_auth` maps. The node fails
if `target_chain` matches no chain under either numbering, listing what is present.

**Scientific notes:** Author numbering is what papers, figures and viewers show; `label_seq_id` is
what the mmCIF polymer sequence uses and what RFD4 reads. They diverge whenever a construct is
numbered to match a parent protein or a crystallised fragment starts partway in. Emitting both,
once, up front is what lets every later node be unambiguous. Non-polymer atoms (ligands, waters)
have no `label_seq_id` and are excluded from the inventory, since they cannot appear in a contig.
PDB input is the weaker path: it carries no polymer definition, so modified residues are left out of
the numbering and a chain containing them will number differently than the same entry's RCSB mmCIF.
Prefer mmCIF when the target has modified residues. Repeated author numbers (insertion codes) make
author->label translation ambiguous; the node warns and maps the first occurrence, and node 02
repeats that warning against the parameter that names the residue.

**Outputs:**
- `target.cif` — the target, normalised to mmCIF
- `target_inventory.json` — per-chain residue inventory with both numberings and the mapping

### Node 02: Build Condition Spec

**Goal:** Render and validate the RFD4 condition spec before any GPU time is spent.

**Process:** Parses `hotspot_residues` / `motif_residues` (forms `A45`, `A45-50`, `45`), translates
them to `label_seq_id` via node 01's maps unless `residue_numbering=label`, and builds the contig and
conditions for the chosen `task`, then dumps the spec as YAML. The spec is keyed by an example id
and each condition is keyed by its value, matching upstream: `C_CRDGRP: {1: [{select: ...}]}`,
`C_SEQ: {true: [{select: ...}]}`. Every task with a target holds it using `C_CRDGRP` over backbone
atoms plus `C_SEQ`, which is upstream's own binder recipe — without it the target is not actually
held. (`monomer` has no target and sets no conditions.) A `contig` parameter, if set, replaces the derived one verbatim. Every referenced residue and chain
is checked against the inventory first, so a typo fails here with a message naming the valid range,
and a chain whose contig span covers unresolved residues produces a warning naming the holes.
The node then writes `spec_preview.cif`: the target plus a comment block recording each condition's
selection, which is replaced by the richer output of `resolve_condition_spec_inputs()` — the same
code path generation uses to turn a spec into an annotated CIF. A spec that parser rejects fails the
node here, on CPU, in seconds.

**Scientific notes:** What each condition constrains is the whole design: `C_HOT` nominates target
*atoms* the binder should engage without dictating how — it is `condition_hotspot_1_atom`, and
naming a residue without naming atoms marks every atom in it, backbone included; `C_CRD` freezes coordinates; `C_CRDGRP`
freezes a set of residues relative to each other while letting the group move; `C_SEQ` freezes
identity. A motif graft therefore needs `C_CRDGRP` *and* `C_SEQ` — geometry without identity
scaffolds the right shape with the wrong chemistry. When a design comes out wrong, `spec_preview.cif`
is the file to open: it shows exactly which residues the model was told to hold. Author
numbering cannot name a residue carrying an insertion code -- `100A` is not a number -- so a
hotspot or motif given in author numbering on such a chain conditions the first occurrence only.
The node says so where the residue is named, listing the residues it could not reach and the label
numbers that would reach them; on an antibody those are the CDR residues themselves.

**Outputs:**
- `design_spec.yaml` — the condition spec handed to `rfd4-generate`
- `spec_preview.cif` — the target annotated with the conditioned residues
- `target.cif` — passed through for node 03

### Node 03: Generate Designs

**Goal:** Generate `num_replicates` designs against the spec.

**Process:** Copies the spec and `target.cif` into a scratch directory (the spec's `input:` is a
bare relative filename) and runs `rfd4-generate` there as an argv list via `subprocess.Popen` without
`shell=True`, with `HF_HUB_OFFLINE=1` so the baked-in weights are never re-fetched. Generated
structures are then flattened out of `inference_samples/` into the node's outputs, with a warning if
that directory is missing or the count does not match `num_replicates`. `--dry-run` prints the
resolved argv and exits 0 without touching the GPU.

**Scientific notes:** `seed` plus the pinned checkpoint revision are what let a run repeat; the flow
model and the autoencoder are a *matched pair*, so pinning one revision for both is not optional.
Note that upstream leaves `torch_deterministic` false and this node does not change it, so a repeat
is empirical rather than guaranteed — same seed and same hardware reproduced designs byte for byte
in testing, but a different GPU or driver may not. `disable_compile` trades a slow `torch.compile` warm-up for slower per-design inference,
which is the right trade for a handful of designs and the wrong one for hundreds.

**Outputs:**
- `*_generated.cif` — one structure per design
- `generation_manifest.json` — resolved argv, seed, checkpoint revision, package version, GPU,
  elapsed time, and the unscreened-designs note

### Node 04: Analyze Designs

**Goal:** Turn the generated structures into a comparable per-design metrics table.

**Process:** Reads `_atom_site` from each design and identifies the designed chain as the one that
was not in the target: node 03 records the target's chain ids in `generation_manifest.json`, and the
chain of a generated file whose id is not among them is the design. The rest are target. It then
computes length, sequence, hydrophobic fraction, net charge, CA radius of gyration, and a CA–CA
interface at an 8 Å cutoff: total contacts, distinct target residues contacted, and minimum
interchain distance. Every row records `designed_chain_source`, so a design that fell back to the
first chain — no manifest, or chain ids that do not separate the chains cleanly — is distinguishable
from one where the question was actually answered. All chains and their lengths are reported
alongside, so a changed output layout is visible rather than silently mis-attributed.

**Scientific notes:** Generation runs with `~metrics`, so every number here is geometric — computed
from coordinates, not predicted. They are useful as filters and sanity checks (a binder with zero
contacts missed the target; an extreme hydrophobic fraction suggests an aggregation-prone surface;
Rg far off the expected range for the length suggests an unfolded or collapsed backbone), and
useless as confidence estimates. Unrecognised residues are scored `X` and excluded from composition
so they cannot dilute the ratios.

**Outputs:**
- `design_metrics.csv` — one row per design
- `design_metrics.json` — the same table plus run provenance, for the report

### Node 06: Export FASTA

**Goal:** The designed sequences in the form workflow-026 reads, one file per design.

**Process:** Reads the sequence of each design out of its mmCIF (not out of the metrics table, so
the structures stay authoritative) and writes one FASTA record per design to `sequences.fasta`,
wrapped at 60 columns, ordered as node 04 ordered them, plus one single-design file
`refold_<design>.fasta` each. Each record holds the designed chain alone, and which chain that is
comes from node 04's `designed_chain` rather than from its position in the file. Record ids are made
filesystem-safe because ColabFold names its outputs after them. Designs with an empty or
all-unknown sequence are skipped with a warning rather than written as junk.

**Scientific notes:** A generated backbone plus sequence is a hypothesis, and the standard test is
self-consistency: fold the sequence independently and measure the prediction against the designed
backbone. **This node is not that handoff and cannot be**, because only the sequence travels — what
comes back from a predictor is confidence, which a wrongly-folded sequence can also earn.
Self-consistency needs the backbone, so the artifacts a validation step has to read are node 03's
`*_generated.cif`, which carries the designed chain and the target together in one file, and node
04's `design_metrics.json`, which names the designed chain so the validator reads it instead of
re-deriving it. The FASTA export stays, because folding a sequence is a useful thing to do on its
own — it is simply not the input a scoring step takes. Unknown residues (`X`) are
counted and warned about, since ColabFold treats them as sequence gaps and they depress confidence
for reasons that have nothing to do with the design.

**Outputs:**
- `sequences.fasta` — every design, one record each; the batch record, not the handoff
- `refold_<design>.fasta` — one design per file; this is what workflow-026 takes, one run per file

### Node 05: Report

**Goal:** One HTML report carrying its own data.

**Process:** Inlines the metrics table and the top `top_n` structures into a single `report.html`
(Bootstrap 5, Plotly and Mol* from CDN, all data as inline JS), with a click-a-row-to-view table, a
length-vs-contacts scatter, a Mol* viewer switchable across designs, and a run-provenance list.
Designs are listed in name order with clashing ones moved to the bottom, and the report says in so
many words that the order is not a ranking. Nothing in this workflow has measured design quality, and
every column in the table is a geometric description an unfolded chain can score well on — ordering
on one would read as a verdict the numbers cannot support. The biosafety notice is rendered as a
top-level alert, not a footnote.

**Scientific notes:** The report carries its own *data*: structures and metrics are embedded
rather than referenced, so it survives being emailed or archived away from the run directory.
That makes file size scale with `top_n` — roughly half a megabyte per structure of this size —
and the node warns past 50 MB.

It is not, however, self-contained in the stricter sense: Bootstrap, Plotly and Mol* are fetched
from a CDN at view time, so without a network the page renders unstyled, the scatter is missing
and the viewer is empty. The table, the run provenance and the biosafety alert all still render —
each block is drawn independently and provenance is drawn before the chart, so a chart library
that fails to load can no longer take the screening disclaimer down with it.

**Outputs:**
- `report.html` — metrics table, plots and a Mol* viewer, no external data files

## Parameters

### task

| Value / Range | Description |
|---------------|-------------|
| `binder` (default) | Design a new chain against the fixed target chain. Contig: `<binder_length>, /0, <target range>`. Optional `hotspot_residues` become `C_HOT`. The common case. |
| `motif_scaffold` | Hold `motif_residues` as a rigid group (`C_CRDGRP`) with fixed identity (`C_SEQ`) and build scaffold around it. `binder_length` becomes the length range of each flank. Requires `motif_residues`. |
| `monomer` | Designs a standalone chain of `binder_length` residues. The target structure is ignored: node 02 writes an annotated scaffold in place of `target.cif`, because a contig made only of expandable segments is rejected upstream. Mechanically verified, but see the quality note below. |
| `binder_redesign` | Discard `redesign_chain` and generate a replacement against the same epitope: the atoms the existing binder contacts (any heavy atom within `epitope_contact_cutoff_a`) become `C_HOT`, thinned to the number upstream would have marked, so the new chain is aimed where the old one bound. Set `binder_length` to the existing binder's length unless you mean to change it, and set `hotspot_residues` to override the derived epitope. **Needs a complex whose binder is an ordinary polypeptide** — see the worked example below. |

**Trade-off:** More conditioning means more control and less diversity. `binder` with no hotspots
explores the whole target surface; adding hotspots focuses it at the cost of discarding
alternative epitopes.

**Quality note on `monomer`.** It is the least constrained task and, at the defaults, the one whose
output most needs filtering. A single 60-90 residue test design came out 70 residues as asked but
extended (Rg 21.7 Å, where a compact 70-mer is nearer 12 Å) and highly charged (net -23, with a
poly-glutamate tail). Nothing in this workflow scores or filters designs, so treat monomer runs as
raw samples: generate many (`num_replicates` in the hundreds), then refold and filter downstream
(workflow-026) rather than reading any single design as a result. The conditioned tasks
(`binder`, `motif_scaffold`, `binder_redesign`) are anchored to a real structure and behave better.

**Test vs production:** The task does not change with scale — `num_replicates` does.

### residue_numbering

| Value / Range | Description |
|---------------|-------------|
| `auth` (default) | Hotspot/motif residues are given in PDB author numbering, as printed in papers and shown by viewers. Node 02 translates to `label_seq_id`. |
| `label` | Residues are already given as mmCIF `label_seq_id`; passed through untranslated. |

**Trade-off:** None — this is a statement about your input, not a quality knob. Getting it wrong
silently conditions the wrong residues.

### epitope_contact_cutoff_a

- **Type:** float (ångström)
- **Default:** 4.5
- **Description:** Heavy-atom distance at which a target atom counts as contacted by the existing
  binder, used only by `binder_redesign` to read off the epitope to aim at.
- **Guidance:** 4.5 Å between any two heavy atoms is upstream's own definition of a hotspot — the
  criterion RFD4 was trained with — so there is a reason to prefer it over a general convention.
  Lower it towards 4 Å for only the tightest contacts; raise it towards 5–6 Å to admit the rim of
  the interface. Node 02 warns outside 3–8 Å: below 3 Å is shorter than a hydrogen bond, and by
  8 Å it is the whole face of the target. Set `hotspot_residues` instead if you know the epitope
  and would rather name it than derive it.

  This replaced an earlier `epitope_cutoff_a`, which measured **CA-CA** at 8.0 Å. The two are not
  interchangeable: 8 Å is the conventional **Cβ-Cβ** threshold rather than the CA-CA one (~10 Å),
  and on the 5O45 test target the CA-CA rule found 10 contact residues where the heavy-atom rule
  finds 17, missing 8 real contacts. Because the same number means different things in the two,
  node 02 refuses the old parameter rather than reinterpreting a value carried over.

### max_hotspot_atoms

- **Type:** integer
- **Default:** 6
- **Description:** Ceiling on how many atoms are marked as `C_HOT`.
- **Guidance:** `C_HOT` is an atom-level annotation. Upstream samples the number of hotspot atoms
  from a geometric distribution averaging about 5, caps it at 10% of the interface atoms, and
  gives none at all 5% of the time. A derived epitope is capped by both this ceiling and that 10%
  rule, whichever is smaller — on 5O45 the interface is 69 atoms, so the cap is 6. Raising this
  conditions the model harder than anything it saw in training; lowering it leaves the design more
  freedom. Node 02 prints how many atoms it marked and what the 10% cap was.

### num_replicates

- **Type:** integer
- **Default:** 4
- **Description:** How many designs to generate in one run.
- **Guidance:** 2–4 to prove the spec is right; 100+ for a real campaign, since generation is cheap
  relative to the downstream filtering and most designs are discarded. Raise `disable_compile` to
  `false` once the batch is large enough to amortise the warm-up.

  **Raise `generate_timeout_s` with it.** All replicates are produced by one `rfd4-generate` call
  under one wall-clock budget, and the default budget is 3600 s. Measured on an L40S at
  `binder_length=50-120`, a replicate costs roughly 35–50 s after load, so 100 replicates needs
  something like 3500–5000 s and will hit the default. See `generate_timeout_s` below.

### generate_timeout_s

- **Type:** integer (seconds)
- **Default:** 3600
- **Description:** Wall-clock budget for the whole `rfd4-generate` call. On expiry the process
  group is ended and node 03 exits 1.
- **Guidance:** The budget covers every replicate in one call, not one design, so it has to scale
  with `num_replicates`. Roughly 35–50 s per replicate at `binder_length=50-120` on an L40S, plus
  60–90 s of model load, plus several minutes of `torch.compile` warm-up if `disable_compile` is
  `false`. Allow headroom: longer binders and larger targets both cost more.

  The budget exists because a CUDA out-of-memory condition does not stop generation. The caching
  allocator logs it as a warning and the process carries on holding the GPU without making
  progress, so without a budget the node would hang indefinitely. When the output contains an OOM
  marker, node 03 says so in the failure message instead of reporting a bare timeout.

### binder_length

- **Type:** string (`N` or `N-M`)
- **Default:** `50-120`
- **Description:** Length, or length range, of the generated segment. Expandable in the contig, so
  the model picks within the range.
- **Guidance:** 50–120 is a sensible mini-binder window. Narrow it when you have a length target
  (expression, format constraints); widen it to let the model choose. For `motif_scaffold` this is
  the per-flank range, and for `binder_redesign` it should match the existing chain.

### hotspot_residues / motif_residues

- **Type:** string (comma-separated; `A45`, `A45-50` and bare `45` all parse)
- **Default:** `""`
- **Description:** Target residues to engage (`C_HOT`) or motif residues to scaffold
  (`C_CRDGRP` + `C_SEQ`). A bare number uses `target_chain`.
- **Guidance:** Give hotspots in the numbering your source uses and set `residue_numbering` to
  match. Two or three surface residues on the epitope of interest is typical; conditioning a large
  patch over-constrains the design.

  Hotspots must name residues on a chain the contig preserves. With the derived contigs that is
  `target_chain` and nothing else — in particular not the binder, which is what is being designed
  and is not in the input. Naming an absent chain is refused, because such a selection parses,
  matches nothing, and would otherwise produce an unaimed design in silence. If you supply your
  own `contig`, node 02 reads it and refuses only chains that contig does not carry.

  You name residues; node 02 marks **atoms**, picking the side-chain atoms that reach furthest
  from CA — up to two per residue, fewer where the residue has fewer (alanine gives one), and CA
  itself for glycine. That matches upstream, whose own example is `sel('A/*/56/[CG,OH]')` — "the
  PD-L1 side chain atoms to engage". It prints how many atoms it marked, and warns past
  `max_hotspot_atoms`.

  The **derived** epitope that `binder_redesign` computes does not apply that side-chain
  restriction, and the difference is deliberate. It draws from every heavy atom within
  `epitope_contact_cutoff_a`, which is upstream's eligible set, so it can mark a main-chain
  carbonyl oxygen — on 4ZQK two of its six marks are backbone O. The named path has no binder to
  measure against, so "which atoms of this residue would a binder engage" has to be answered from
  the residue alone, and the end of the side chain is the only defensible answer. Whether RFD4
  responds to a marked backbone carbonyl as it does to a marked side-chain atom is untested here.

### contig

- **Type:** string
- **Default:** `""`
- **Description:** Advanced escape hatch. When set it replaces the derived contig verbatim. Not
  accepted with `task=monomer`, which uses a scaffold rather than a contig.
- **Guidance:** Use it for layouts the four tasks do not express — several designed segments, or a
  different span of the target than `target_chain` gives you. Residue ids are `label_seq_id`, and
  no translation is applied to this string.

  **A contig that preserves more than one input chain does not currently work.** The task builders
  hold the target chain alone (`C_CRDGRP` + `C_SEQ`), so any other preserved chain is left without
  a coordinate group and upstream groups it with the designed chain, which fails validation:

  ```
  [index_group_molecule_must_be_connected] ... token pair(s) of one index group have no path
  between them in the bond graph though they are one molecule: C/1-B/1, C/1-B/2, ...
  ```

  This is caught on CPU in node 02, before any GPU time, but there is no parameter that fixes it —
  holding a second chain would need a change to the builders. Keep `contig` to one preserved input
  chain.

  Note also that node 04 attributes the **first** chain to the design and the rest to the target,
  so a contig with several designed segments makes its per-design metrics and interface counts
  wrong; node 04 warns when it sees more than two chains.

### target_chain / redesign_chain

- **Type:** string
- **Default:** `A` / `B`
- **Description:** Chain to design against or preserve, and the chain to regenerate for
  `binder_redesign`. Either numbering's chain id is accepted and resolved against the inventory.

### seed / hf_ckpt_revision

- **Type:** integer / string
- **Default:** `42` / `a06e2509c7033b010fc43614e9fd840303588502`
- **Description:** Provenance controls, both recorded in `generation_manifest.json`.
- **Guidance:** Change `seed` to sample different designs from the same spec. Leave the revision
  alone unless you rebuild the image: the flow model and autoencoder are a matched pair, and the
  weights are baked in at this revision.

  `seed` does not make generation bit-reproducible on its own. Upstream's `inference.yaml` leaves
  `torch_deterministic` false — that is the switch that locks cuDNN and cuBLAS to deterministic
  algorithms — and this node does not override it. In practice the same seed on the same GPU and
  driver reproduced designs byte for byte; across different hardware, treat the seed as provenance
  rather than a guarantee.

### batch_size

- **Type:** integer
- **Default:** `1`
- **Description:** How many designs are generated per forward pass.
- **Guidance:** This is the knob to reach for when a larger target exhausts GPU memory: lowering it
  (or leaving it at 1) trades throughput for a smaller peak allocation. Raising it speeds up large
  batches on a card with memory to spare. If generation hits an allocator out-of-memory condition,
  node 03 now says so explicitly rather than stalling, and this is the first parameter to change.

### disable_compile

- **Type:** boolean
- **Default:** `true`
- **Description:** Sets `compile_lengths=null`, skipping the `torch.compile` warm-up.
- **Guidance:** Leave `true` for small runs; set `false` for large batches where the one-off warm-up
  pays for itself.

### top_n

- **Type:** integer
- **Default:** 10
- **Description:** Structures embedded in the report's viewer. The metrics table always covers every
  design.
- **Guidance:** Lower it if the report gets unwieldy — each structure adds roughly half a megabyte.

## Outputs and interpretation

### *_generated.cif

One mmCIF per design, containing the generated chain plus the target chains it was conditioned on.
Coordinates and sequence come from the model jointly. Nothing validates them: the file is a
hypothesis, not a prediction of a complex. Treat the generated chain's sequence as a starting point —
most campaigns re-design the sequence onto the backbone and refold before ordering anything.

### design_metrics.csv

One row per design. `designed_length` and `designed_sequence` describe the generated chain;
`hydrophobic_fraction` (fraction of recognised residues in AVILMFWYC) and `net_charge`
(K+R minus D+E) summarise its composition; `radius_of_gyration_A` is computed over CA atoms.
`n_interface_contacts`, `n_target_residues_contacted` and `min_interchain_ca_distance_A` describe the
interface at an 8 Å CA–CA cutoff — a coarse definition that counts proximity, not interaction
energy. They are **blank for a single-chain design** (`monomer`), where there is no interface to
measure, rather than reported as zero. For a two-chain design, zero contacts means it did not engage
the target and can be dropped immediately; beyond that, these numbers rank nothing. `chain_lengths`
records every chain and `designed_chain_source` records how the designed chain was identified, so an
unexpected layout shows up as a named fallback rather than as quietly wrong metrics.

`chain_breaks` counts consecutive alpha carbons further apart than a peptide bond reaches, and
`compactness_ratio` is `radius_of_gyration_A` over `expected_radius_of_gyration_A`, the latter being
2.2·N^0.38 Å for a compact fold. These two exist to catch output that is not a folded polypeptide at
all — a wrong checkpoint, a failed generation. Nothing else here would: contacts, charge, composition
and even Rg itself are all computable on an unfolded chain and all look ordinary. A denatured chain
sits near ratio 2.1; designs from this workflow have measured 0.94 to 1.49, and the warning fires
above 1.8. They say nothing about whether a folded design is any *good*.

### sequences.fasta

One FASTA record per design, holding that design's generated chain and named after the design file
so results can be traced back.

**Hand workflow-026 one `refold_<design>.fasta` at a time, renamed to `sequences.fasta` — not the
batch file.** Workflow-026 reports one job per run. Its prediction step is fine with several
records, since ColabFold makes each one its own job, but every step after it assumes a single job:
`03_analyze_confidence` globs the score files of all designs together, takes the job id from
whichever sorts first, then sorts by model rank alone and truncates. With three designs and five
models the reported "top 5" are ranks 1,1,1,2,2 drawn from three different designs, filed under one
design's name, with the three designs' lengths read as one job's three chains. It does not fail —
it produces a confident report about a design that does not exist.

Every record's description carries `UNSCREENED - computational design, requires biosecurity
review`. This file is the only artifact built to leave the workflow, so it is the one that has to
carry the warning: the report and the manifest both record it, and neither travels with the FASTA.
The tag survives workflow-026's validator into its `validated_sequences.fasta`.

It is a description rather than a leading `;` comment line on purpose. Biopython still accepts
comment lines but warns that it will raise `ValueError` in a future release, and workflow-026
calls `SeqIO.parse` without a `try`/`except` — a banner would break workflow-026 the next time
its image was rebuilt, and would look like a bug in workflow-026. The tag also contains no colon,
because workflow-026 reads a colon in a description as a chain-break marker.

There is no complex-export mode. Colon-joining the chains is the only encoding that would make
ColabFold predict the complex, because workflow-026 passes the FASTA directory straight to
`colabfold_batch` and treats each record as an independent job, so a chain break signalled in the
header would be ignored. But workflow-026 validates sequences against `ACDEFGHIKLMNPQRSTVWYX` and
rejects anything else, so a colon-joined record fails its first node outright. The two
requirements are incompatible from this side; re-enabling complex export means adding `:` to that
validator's `VALID_AA` first.

### report.html

The metrics table, a length-vs-contacts scatter, a Mol* viewer over the top `top_n` designs, and
the run's provenance (seed, checkpoint revision, package version, GPU). All the run's own data is
embedded, so it opens from anywhere; Bootstrap, Plotly and Mol* are fetched from a CDN at view
time, so styling, the scatter and the viewer need a network. The table, the provenance and the
biosafety alert do not.

## Biosafety

RFD4's HuggingFace Spaces apps screen every input and design with SecureDNA. **The CLI and Python
API — the path this workflow uses — do not.** Designs produced here are therefore unscreened. They
require biosecurity review before any sequence is ordered or synthesised. This workflow deliberately
contains no screening node; screening is the user's responsibility. `report.html` and
`generation_manifest.json` both record that no screening was performed, and every record in
`sequences.fasta` is tagged `UNSCREENED` so the warning travels with the designs when they are
handed to another workflow.

## Quick start

### Running with Docker

Built from [Dockerfile](Dockerfile) as `rfd4-proteina:2026_11_04`:

```bash
docker build -t rfd4-proteina:2026_11_04 workflows/workflow-037/
```

Node 01, fetching PD-L1 from RCSB:

```bash
mkdir -p /tmp/wf037/01/inputs /tmp/wf037/01/outputs && docker run --rm -v /tmp/wf037/01/inputs:/workspace/inputs -v /tmp/wf037/01/outputs:/workspace/outputs -v $(pwd)/workflows/workflow-037/01-acquire-target:/workspace -w /workspace -e PARAM_PDB_ID=5o45 -e PARAM_TARGET_CHAIN=A rfd4-proteina:2026_11_04 bash run.sh
```

Node 03, resolving the argv without using the GPU:

```bash
docker run --rm -v /tmp/wf037/03/inputs:/workspace/inputs -v /tmp/wf037/03/outputs:/workspace/outputs -v $(pwd)/workflows/workflow-037/03-generate:/workspace -w /workspace rfd4-proteina:2026_11_04 python script.py --dry-run
```

Node 03, for real:

```bash
docker run --rm --gpus all -v /tmp/wf037/03/inputs:/workspace/inputs -v /tmp/wf037/03/outputs:/workspace/outputs -v $(pwd)/workflows/workflow-037/03-generate:/workspace -w /workspace -e PARAM_NUM_REPLICATES=2 -e PARAM_SEED=42 rfd4-proteina:2026_11_04 bash run.sh
```

Chain the nodes by copying each node's `outputs/` into the next node's `inputs/`.

### Running on Silva

1. Select this workflow from the workflow list
2. Upload your target structure (or set `pdb_id`)
3. Adjust parameters if needed (see Parameters section)
4. Click Run

### Test vs production settings

| Setting | Test (default) | Production |
|---------|---------------|------------|
| `num_replicates` | 4 | 100+ |
| `generate_timeout_s` | 3600 | scale with `num_replicates` — ~35–50 s each, so 100+ needs 5400+ |
| `disable_compile` | `true` | `false` (amortises the warm-up over a large batch) |
| `top_n` | 10 | 10–25 (report size scales with it) |

A successful test run against `pdb_id=5o45` with `task=binder` produces four structures whose first
chain is 50–120 residues and has a non-zero contact count against chain A, plus a report whose
viewer shows the designed chain against the target.

### Worked example: binder_redesign on 4ZQK

Use **4ZQK** — the human PD-1 / PD-L1 complex — rather than 5O45 for this task:

```
pdb_id=4ZQK  task=binder_redesign  target_chain=A  redesign_chain=B  binder_length=106
```

Chain A is PD-L1 (115 residues, label 1–115 / author 18–132) and chain B is PD-1 (106 residues
resolved over label 1–114, author 33–146). It keeps the same target protein as the `binder`
example, so the two tasks are directly comparable, and the interface is a real protein–protein
one: 77 heavy-atom contacts within 4.5 Å, giving a cap of 7 hotspot atoms.

**Why not 5O45 here.** Its chain B is not a protein binder at all — it is a 15-residue macrocyclic
peptide, 7 of whose residues are non-standard (`MEA`, `9KK`, `SAR`, `CCS`, `NH2`). Asking RFD4 to
replace it with an ordinary polypeptide is not a redesign of that binder so much as a different
molecule in the same place, and node 01 now warns when a chain is largely non-standard for exactly
this reason. It also makes the CA-CA geometry misleading: the macrocycle's alpha carbons sit on a
ring away from the contact surface, which is why a CA-CA epitope rule misses so much of its
interface. 5O45 remains the right demo for `task=binder`, where chain B is ignored.

## Regression suite

`tests/` runs nodes 01 and 02 over ten real RCSB entries, checks that node 02 refuses
the inputs it should, and drives node 04 with structures built for the purpose. Output
is checked against a measurement of the input, taken by a parser that is not the one the
nodes use.

```sh
./tests/run_fixtures.sh
```

Every entry exits 0 on both nodes, which is why the suite does not look at exit codes.
It asks whether a residue the file contains can actually be addressed, whether an NMR
ensemble was noticed, whether a modified residue survived the PDB conversion, whether a
cofactor was named before being dropped. The first run found twelve such answers wrong.

Fixtures are chosen by measuring a downloaded structure, never from reputation, and
`tests/fixtures.json` records the measurement next to each entry. A finding that is
reproduced but not yet fixed is marked `"status": "open"` there and reports without
failing the run; nothing is marked open today. See `tests/README.md`.

**What this suite does not establish.** It tests the plumbing, not the science. Passing means the
pipeline does not silently corrupt, drop, mislabel or misnumber what it is given, and that it
notices output which is not a folded polypeptide. It does not mean any design is good. Nothing here
refolds a sequence, computes self-consistency, or compares a result against a published one — the
structures are chosen for the input shapes they carry, not because the workflow is known to
reproduce a literature result on them. Design quality stays unvalidated until something refolds a
sequence and measures the prediction against the backbone that was designed — which neither this
workflow nor workflow-026 does. The artifacts such a step would read are `*_generated.cif` and
`design_metrics.json`, not the exported FASTA.

## Troubleshooting

**Node 02: "Condition spec rejected by RFD4"**
Upstream's own parser refused the spec, so the contig or a selection is wrong. Check
`target_inventory.json` for the chains and both numberings actually present, and confirm
`residue_numbering` matches what you typed. No GPU time has been spent at this point.

**Node 02: "Contig has no preserved segment"**
Every contig must reference the input at least once; one built only from expandable lengths is
rejected upstream. This is also why `task=monomer` needs an annotated CIF rather than a length.

**Designs come out contacting the wrong part of the target**
Almost always a numbering mix-up. Node 02 prints every condition it emitted, and `spec_preview.cif`
is RFD4's own annotated rendering of the spec — open it and confirm the conditioned residues are
the epitope you meant.

**Node 03: "No GPU available"**
The node requires a GPU and refuses to run on CPU. Confirm the job is scheduled on a GPU node and
that the container sees it (`nvidia-smi` inside the image).

**Node 02 warns that the target chain is missing residues inside its span**
The chain has unresolved loops (or modified residues, on the PDB path), so the contig segment covers
positions with no coordinates. Usually harmless for a binder run against a well-ordered epitope;
worth a look if the holes sit at the interface you are targeting.

**Node 02 warns that rfproteina is not importable**
Only this node's own comment annotation is written and the spec is not validated. That is the
development path outside the image; inside the image it means the image is broken, and node 02
stops rather than continuing.

**Node 02: "hotspot_residues names chain(s) B, but only the target chain (A) is preserved"**
Hotspots name residues on the target that the designed chain should engage. They cannot name the
binder: with `task=binder` the binder does not exist yet, and with `binder_redesign` the existing
one is replaced. Give the residues on `target_chain`. If you supplied your own `contig`, node 02
reads it and this message instead names the chains that contig preserves.

**Node 02: "The spec asks for N hotspot atom(s) ... but the rendered system flags only M"**
At least one selection parsed but matched nothing, so the design would be aimed at part of the
patch you asked for. Usually a residue id outside the contig, or an atom name that residue does
not have. Check every id against `target_inventory.json` from node 01, and confirm each residue
falls inside the contig node 02 printed.

**Node 02: "epitope_cutoff_a has been replaced by epitope_contact_cutoff_a"**
The old parameter was a CA–CA distance defaulting to 8.0 Å; the new one is a heavy-atom contact
distance defaulting to 4.5 Å, which is the criterion RFD4 was trained with. The same number means
different things under the two rules, so the value is not carried over. Remove the old parameter,
and set the new one only if you need something other than 4.5.

**Node 03: "rfd4-generate exceeded its Ns budget"**
All replicates run under one wall-clock budget. Either the run genuinely needs longer — scale
`generate_timeout_s` with `num_replicates`, see that parameter — or generation stalled. If the
message also reports a CUDA out-of-memory condition, it stalled: the allocator logs OOM as a
warning and keeps the GPU without making progress, so reduce `binder_length` or `num_replicates`,
or use a GPU with more memory.

**Node 03: "hf_ckpt_revision ... is not baked into this image"**
The weights are baked in at one revision and this node runs offline, so no other revision can be
fetched. The message lists what the image has. Rebuild the image to change it.

## References

- RFD4-Proteina — RosettaCommons and NVIDIA GenAIR, 2026. BSD-3-Clause.
- [nvidia/RFD4-Proteina model card](https://huggingface.co/nvidia/RFD4-Proteina)
- Design rationale, verification plan and open build-time questions: [rfd4_pipeline.md](rfd4_pipeline.md)
