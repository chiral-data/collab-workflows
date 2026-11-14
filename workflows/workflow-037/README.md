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
(ColabFold). The intended follow-up is to feed designs from here into workflow-026 and check
self-consistency: node 06 writes `sequences.fasta` for exactly that, but the refolding and the
filtering happen there, not here. The most common mistake is
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
into workflow-026 as its input to refold the designs.

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
author->label translation ambiguous; the node warns and maps the first occurrence.

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
residues the binder should engage without dictating how; `C_CRD` freezes coordinates; `C_CRDGRP`
freezes a set of residues relative to each other while letting the group move; `C_SEQ` freezes
identity. A motif graft therefore needs `C_CRDGRP` *and* `C_SEQ` — geometry without identity
scaffolds the right shape with the wrong chemistry. When a design comes out wrong, `spec_preview.cif`
is the file to open: it shows exactly which residues the model was told to hold.

**Outputs:**
- `design_spec.yaml` — the condition spec handed to `rfd4-generate`
- `spec_preview.cif` — the target annotated with the conditioned residues
- `target.cif` — passed through for node 03

### Node 03: Generate Designs

**Goal:** Generate `num_replicates` designs against the spec.

**Process:** Copies the spec and `target.cif` into a scratch directory (the spec's `input:` is a
bare relative filename) and runs `rfd4-generate` there as an argv list via `subprocess.run` without
`shell=True`, with `HF_HUB_OFFLINE=1` so the baked-in weights are never re-fetched. Generated
structures are then flattened out of `inference_samples/` into the node's outputs, with a warning if
that directory is missing or the count does not match `num_replicates`. `--dry-run` prints the
resolved argv and exits 0 without touching the GPU.

**Scientific notes:** `seed` plus the pinned checkpoint revision are what make a run reproducible;
the flow model and the autoencoder are a *matched pair*, so pinning one revision for both is not
optional. `disable_compile` trades a slow `torch.compile` warm-up for slower per-design inference,
which is the right trade for a handful of designs and the wrong one for hundreds.

**Outputs:**
- `*_generated.cif` — one structure per design
- `generation_manifest.json` — resolved argv, seed, checkpoint revision, package version, GPU,
  elapsed time, and the unscreened-designs note

### Node 04: Analyze Designs

**Goal:** Turn the generated structures into a comparable per-design metrics table.

**Process:** Reads `_atom_site` from each design, treats the first chain in the file as the designed
chain (node 02's contig puts the designed segment ahead of the chain break) and the rest as target,
then computes length, sequence, hydrophobic fraction, net charge, CA radius of gyration, and a CA–CA
interface at an 8 Å cutoff: total contacts, distinct target residues contacted, and minimum
interchain distance. All chains and their lengths are reported alongside, so a changed output layout
is visible rather than silently mis-attributed.

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

**Goal:** The designed sequences in the form workflow-026 reads, so refolding is a file copy.

**Process:** Reads the sequence of each design out of its mmCIF (not out of the metrics table, so
the structures stay authoritative) and writes one FASTA record per design to `sequences.fasta`,
wrapped at 60 columns, ordered as node 04 ordered them. `fasta_mode=designed` writes the generated
chain alone; `fasta_mode=complex` writes every chain colon-joined, which ColabFold reads as a chain
break. Record ids are made filesystem-safe because ColabFold names its outputs after them. Designs
with an empty or all-unknown sequence are skipped with a warning rather than written as junk.

**Scientific notes:** This node is the bridge to the scoring this workflow deliberately omits. A
generated backbone plus sequence is a hypothesis; the standard test is self-consistency — fold the
sequence independently and ask whether the prediction matches the designed backbone. Refolding the
designed chain alone (`designed`) answers that. Refolding the complex (`complex`) additionally gives
interface confidence, at more compute and only when there is a target. Unknown residues (`X`) are
counted and warned about, since ColabFold treats them as sequence gaps and they depress confidence
for reasons that have nothing to do with the design.

**Outputs:**
- `sequences.fasta` — one record per design, ready for workflow-026

### Node 05: Report

**Goal:** One self-contained HTML report.

**Process:** Inlines the metrics table and the top `top_n` structures into a single `report.html`
(Bootstrap 5, Plotly and Mol* from CDN, all data as inline JS), with a click-a-row-to-view table, a
length-vs-contacts scatter, a Mol* viewer switchable across designs, and a run-provenance list.
Designs are ordered by interface contacts descending — a convenience, not a ranking. The biosafety
notice is rendered as a top-level alert, not a footnote.

**Scientific notes:** The report is deliberately standalone: structures are embedded rather than
referenced, so it survives being emailed or archived away from the run directory. That makes file
size scale with `top_n` — roughly half a megabyte per structure of this size — and the node warns
past 50 MB.

**Outputs:**
- `report.html` — metrics table, plots and a Mol* viewer, no external data files

## Parameters

### task

| Value / Range | Description |
|---------------|-------------|
| `binder` (default) | Design a new chain against the fixed target chain. Contig: `<binder_length>, /0, <target range>`. Optional `hotspot_residues` become `C_HOT`. The common case. |
| `motif_scaffold` | Hold `motif_residues` as a rigid group (`C_CRDGRP`) with fixed identity (`C_SEQ`) and build scaffold around it. `binder_length` becomes the length range of each flank. Requires `motif_residues`. |
| `monomer` | Designs a standalone chain of `binder_length` residues. The target structure is ignored: node 02 writes an annotated scaffold in place of `target.cif`, because a contig made only of expandable segments is rejected upstream. Mechanically verified, but see the quality note below. |
| `binder_redesign` | Discard `redesign_chain` and generate a replacement against the same epitope: the residues the existing binder contacts (CA-CA within `epitope_cutoff_a`) become `C_HOT`, so the new chain is aimed where the old one bound. Set `binder_length` to the existing binder's length unless you mean to change it, and set `hotspot_residues` to override the derived epitope. |

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

### epitope_cutoff_a

- **Type:** float (ångström)
- **Default:** 8.0
- **Description:** CA-CA distance at which a target residue counts as contacted by the existing
  binder, used only by `binder_redesign` to read off the epitope to aim at.
- **Guidance:** 8 Å over alpha-carbons is a coarse, conventional interface definition, and it is
  this workflow's choice rather than an upstream recommendation (upstream publishes no redesign
  recipe). Lower it towards 6 Å for a tighter, smaller patch — fewer hotspots, more freedom for the
  design; raise it towards 10-12 Å to include the rim of the interface. Below about 4 Å almost
  nothing qualifies and above about 15 Å it stops meaning contact; node 02 warns outside that band.
  Set `hotspot_residues` instead if you know the epitope and would rather name it than derive it.

### num_replicates

- **Type:** integer
- **Default:** 4
- **Description:** How many designs to generate in one run.
- **Guidance:** 2–4 to prove the spec is right; 100+ for a real campaign, since generation is cheap
  relative to the downstream filtering and most designs are discarded. Raise `disable_compile` to
  `false` once the batch is large enough to amortise the warm-up.

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
  match. Three to five surface residues on the epitope of interest is typical; conditioning a large
  patch over-constrains the design.

### contig

- **Type:** string
- **Default:** `""`
- **Description:** Advanced escape hatch. When set it replaces the derived contig verbatim. Not
  accepted with `task=monomer`, which uses a scaffold rather than a contig.
- **Guidance:** Use it for layouts the four tasks do not express (multiple chain breaks, several
  designed segments). Residue ids are `label_seq_id`, and no translation is applied to this string.
  Note that node 04 attributes the **first** chain to the design and the rest to the target, so a
  contig with several designed segments makes its per-design metrics and interface counts wrong;
  node 04 warns when it sees more than two chains.

### target_chain / redesign_chain

- **Type:** string
- **Default:** `A` / `B`
- **Description:** Chain to design against or preserve, and the chain to regenerate for
  `binder_redesign`. Either numbering's chain id is accepted and resolved against the inventory.

### seed / hf_ckpt_revision

- **Type:** integer / string
- **Default:** `42` / `a06e2509c7033b010fc43614e9fd840303588502`
- **Description:** Reproducibility controls, both recorded in `generation_manifest.json`.
- **Guidance:** Change `seed` to sample different designs from the same spec. Leave the revision
  alone unless you rebuild the image: the flow model and autoencoder are a matched pair, and the
  weights are baked in at this revision.

### disable_compile

- **Type:** boolean
- **Default:** `true`
- **Description:** Sets `compile_lengths=null`, skipping the `torch.compile` warm-up.
- **Guidance:** Leave `true` for small runs; set `false` for large batches where the one-off warm-up
  pays for itself.

### fasta_mode

| Value / Range | Description |
|---------------|-------------|
| `designed` (default) | One record per design holding the generated chain only. Refolding it tests self-consistency: does this sequence fold to the shape it was designed as? |
| `complex` | One record per design holding every chain, colon-joined, designed chain first. ColabFold folds the complex, which also reports interface confidence. Slower, and identical to `designed` for a single-chain `monomer` design. |

**Trade-off:** `designed` is the cheaper and more direct test of the design itself; `complex` tells
you about the interaction but costs a multimer prediction per design.

**Test vs production:** Either is fine for a smoke test. For a real campaign, `designed` first to
discard designs that do not fold, then `complex` on the survivors.

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
`n_interface_contacts`, `n_target_residues_contacted` and `min_interchain_distance_A` describe the
interface at an 8 Å CA–CA cutoff — a coarse definition that counts proximity, not interaction
energy. They are **blank for a single-chain design** (`monomer`), where there is no interface to
measure, rather than reported as zero. For a two-chain design, zero contacts means it did not engage
the target and can be dropped immediately; beyond that, these numbers rank nothing. `chain_lengths` records every chain, so an unexpected layout shows
up as a changed first chain rather than as quietly wrong metrics.

### sequences.fasta

One FASTA record per design, named after the design file so results can be traced back. Hand it to
workflow-026 unchanged — the filename is the one that workflow expects. Note that workflow-026
labels any multi-record file a "multimer" and so resolves its `pair_mode=auto` to
`unpaired_paired`; that is harmless for single-chain records and correct for colon-joined ones, so
no adjustment is needed on the way over.

### report.html

Self-contained: the metrics table, a length-vs-contacts scatter, a Mol* viewer over the top `top_n`
designs, and the run's provenance (seed, checkpoint revision, package version, GPU). Opens from
anywhere, with only the CDN assets fetched at view time.

## Biosafety

RFD4's HuggingFace Spaces apps screen every input and design with SecureDNA. **The CLI and Python
API — the path this workflow uses — do not.** Designs produced here are therefore unscreened. They
require biosecurity review before any sequence is ordered or synthesised. This workflow deliberately
contains no screening node; screening is the user's responsibility, and `report.html` and
`generation_manifest.json` both record that no screening was performed.

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
| `disable_compile` | `true` | `false` (amortises the warm-up over a large batch) |
| `top_n` | 10 | 10–25 (report size scales with it) |

A successful test run against `pdb_id=5o45` with `task=binder` produces four structures whose first
chain is 50–120 residues and has a non-zero contact count against chain A, plus a report whose
viewer shows the designed chain against the target.

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

## References

- RFD4-Proteina — RosettaCommons and NVIDIA GenAIR, 2026. BSD-3-Clause.
- [nvidia/RFD4-Proteina model card](https://huggingface.co/nvidia/RFD4-Proteina)
- Design rationale, verification plan and open build-time questions: [rfd4_pipeline.md](rfd4_pipeline.md)
