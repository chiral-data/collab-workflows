---
doc_id: workflow-001-wip
domain: protein-protein-docking
doc_type: workflow
version: "0.1.0"
deprecated: false
description: >
  (WIP) Protein-protein docking combining Boltz structure prediction with
  DiffDock-PP docking. Not yet functional.
tags: [boltz, diffdock-pp, protein-protein-docking, wip]
---

# Workflow 001: Protein-Protein Docking (WIP)

This workflow is **incomplete and not yet functional**. It is intended to combine Boltz structure prediction with DiffDock-PP protein-protein docking in a two-node pipeline.

## Current status

- **Node 01 (Boltz):** Has a container image (`ghcr.io/chiral-data/boltz:2025_09_05`) and a placeholder `run.sh` that only calls `boltz --help`. No actual prediction logic.
- **Node 02 (DiffDock-PP):** Has a `job.toml` with no container image and an empty `run.sh`. Not implemented.

## Alternatives

For protein-protein docking, use one of these functional workflows instead:

- **workflow-016** — DiffDock-PP with ESM-2 embeddings (deep learning scoring)
- **workflow-017** — LightDock with glowworm swarm optimization (classical docking)
- **workflow-012** — Boltz-2 structure prediction (predicts complex structure from sequence, does not dock pre-existing structures)
