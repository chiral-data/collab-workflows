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
scaffolding, antibodies and enzymes. This workflow takes a target structure and a small set of guided
parameters, renders an RFD4 condition spec, validates it on CPU, and generates designs on one GPU.
It answers "give me candidate backbones that bind or scaffold this thing" — it does not score them.

## Overview

<!-- TODO: expand — the flow-matching + autoencoder architecture, the condition-spec
     conditioning vocabulary (C_CRDGRP / C_SEQ / C_HOT / C_CRD), what a generated
     backbone is and is not, and the RFD4-Proteina citation. 1-2 paragraphs. -->

## When to use this workflow

<!-- TODO: paragraph 1 — target structure available, want candidate backbones for a
     binder/motif/monomer/redesign task, GPU available (A100/L40S/RTX6000Ada/H100/
     H200/Blackwell). -->

<!-- TODO: paragraph 2 — when NOT to: this workflow is generation only. It emits no
     confidence or self-consistency scores. For folding a sequence use workflow-026
     (ColabFold); for an end-to-end scored binder campaign use workflow-033
     (RFdiffusion → ProteinMPNN → Boltz2) or workflow-013 (BoltzGen). Note the
     label_seq_id vs auth numbering mistake. -->

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
   └─────────────────────────┴────────► [05-report] ──► report.html
```

Nodes run sequentially 01 → 02 → 03 → 04 → 05; node 05 consumes the generated structures from
node 03 directly as well as the metrics from node 04. Only node 03 requires a GPU.

## Input requirements

One target structure, supplied either way:

- **Upload** a `.cif`, `.pdb` or `.cif.zst` file into the workflow's input directory, or
- **Set `pdb_id`** to an RCSB accession and let node 01 fetch it.

An upload always wins over `pdb_id`. See [input_files/README.md](input_files/README.md) for the
PD-L1 (5O45) sample target used by the correctness gate.

**Residue numbering matters.** `contig` and `conditions` in an RFD4 spec refer to the original
mmCIF `label_seq_id`, not author numbering. RCSB files are author-numbered, so a hotspot copied
from a paper is almost always author numbering — leave `residue_numbering` at `auth` and node 02
translates it.

## Workflow nodes

### Node 01: Acquire Target

**Goal:** Produce one normalised mmCIF target plus an inventory of its chains and both residue
numberings.

**Process:** <!-- TODO --> 

**Scientific notes:** <!-- TODO: the label_seq_id ↔ auth_seq_id mapping and why it is emitted. -->

**Outputs:**
- `target.cif` — the target, normalised to mmCIF
- `target_inventory.json` — per-chain residue inventory with both numberings

### Node 02: Build Condition Spec

**Goal:** Render and validate the RFD4 condition spec before any GPU time is spent.

**Process:** <!-- TODO --> 

**Scientific notes:** <!-- TODO: what each condition type constrains; why spec_preview.cif is
     the artifact to inspect when a design comes out wrong. -->

**Outputs:**
- `design_spec.yaml` — the condition spec handed to `rfd4-generate`
- `spec_preview.cif` — annotated CIF showing exactly which atoms are conditioned
- `target.cif` — passed through for node 03

### Node 03: Generate Designs

**Goal:** Generate `num_replicates` designs against the spec.

**Process:** <!-- TODO --> 

**Scientific notes:** <!-- TODO --> 

**Outputs:**
- `*_generated.cif` — one structure per design
- `generation_manifest.json` — resolved argv, seed, checkpoint revision, timings

### Node 04: Analyze Designs

**Goal:** Turn the generated structures into a comparable per-design metrics table.

**Process:** <!-- TODO --> 

**Scientific notes:** <!-- TODO: generation runs with `~metrics`, so every number here is
     geometric; none of it is a confidence estimate. -->

**Outputs:**
- `design_metrics.csv` — one row per design
- `design_metrics.json` — the same table for the report

### Node 05: Report

**Goal:** One self-contained HTML report.

**Process:** <!-- TODO --> 

**Scientific notes:** <!-- TODO --> 

**Outputs:**
- `report.html` — metrics table, plots and a Mol* viewer, no external data files

## Parameters

### task

| Value / Range | Description |
|---------------|-------------|
| `binder` (default) | <!-- TODO --> |
| `motif_scaffold` | <!-- TODO --> |
| `monomer` | <!-- TODO --> |
| `binder_redesign` | <!-- TODO --> |

**Trade-off:** <!-- TODO -->

**Test vs production:** <!-- TODO -->

### residue_numbering

| Value / Range | Description |
|---------------|-------------|
| `auth` (default) | Hotspot/motif residues are given in PDB author numbering, as printed in papers and shown by viewers. Node 02 translates to `label_seq_id`. |
| `label` | Residues are already given as mmCIF `label_seq_id`; passed through untranslated. |

**Trade-off:** None — this is a statement about your input, not a quality knob. Getting it wrong
silently conditions the wrong residues.

<!-- TODO: add sections for pdb_id, target_chain, binder_length, contig, hotspot_residues,
     motif_residues, redesign_chain, num_replicates, seed, disable_compile,
     hf_ckpt_revision, top_n. -->

## Outputs and interpretation

### *_generated.cif

<!-- TODO: what a generated backbone is — coordinates and a designed sequence from the
     model, with no independent validation that it folds or binds. -->

### design_metrics.csv

<!-- TODO: define each column, its typical range, and its caveats. -->

## Biosafety

RFD4's HuggingFace Spaces apps screen every input and design with SecureDNA. **The CLI and Python
API — the path this workflow uses — do not.** Designs produced here are therefore unscreened. They
require biosecurity review before any sequence is ordered or synthesised. This workflow deliberately
contains no screening node; screening is the user's responsibility and the report repeats this
notice.

## Quick start

### Running with Docker

Built from [Dockerfile](Dockerfile) as `rfd4-proteina:2026_11_04`.

```bash
docker build -t rfd4-proteina:2026_11_04 workflows/workflow-037/
```

<!-- TODO: per-node docker run commands once the nodes are implemented. -->

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

<!-- TODO: one sentence on the included test data and what a successful run looks like. -->

## Troubleshooting

**`sel()` selector or contig error from node 02**
The spec references residues or chains that do not exist under the numbering in use. Check
`target_inventory.json` for the chains and both numberings actually present, and confirm
`residue_numbering` matches what you typed.

<!-- TODO: add GPU/driver and weight-cache entries after the build checks run. -->

## References

- RFD4-Proteina — RosettaCommons and NVIDIA GenAIR, 2026. BSD-3-Clause. <!-- TODO: paper/DOI -->
- [nvidia/RFD4-Proteina model card](https://huggingface.co/nvidia/RFD4-Proteina)
