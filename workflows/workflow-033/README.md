---
doc_id: workflow-033
domain: protein-binder-design
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  De novo protein binder design pipeline that composes RFdiffusion, ProteinMPNN,
  and Boltz2 NVIDIA NIM services to generate, sequence-design, co-fold, and rank
  novel protein binders against a target structure and hotspot epitope.
tags: [protein-binder-design, rfdiffusion, proteinmpnn, boltz2, nvidia-nim, structure-prediction]
---

# Workflow 033: Protein Binder Design — De Novo NIM Pipeline

End-to-end de novo protein binder design campaign using three NVIDIA BioNeMo NIM services. Given a target protein PDB and a set of hotspot residues, the workflow generates novel binder backbones with RFdiffusion, designs sequences with ProteinMPNN, co-folds each binder–target complex with Boltz2, and ranks candidates by interface confidence, pLDDT, and self-consistency RMSD.

---

## When to use this workflow

Use this workflow when you want to design novel protein binders against a known target structure and have a set of interface residues (hotspots) to condition the design on. Requires an `NVIDIA_API_KEY` for hosted NIM access.

Do not use this workflow for small-molecule docking (use workflow-025), single-chain structure prediction (use workflow-026), or when you have a known binder sequence to refine.

---

## Architecture and data flow

```text
global_params.json
       |
[01: target-prep] ──────────────────────────────────────┐
       |                                                 |
  target.pdb                                        target.pdb
  chain_seq.txt                                     target_a3m.txt
  hotspots.json                                          |
  target_a3m.txt                                         |
       |                                                 |
[02: generate-backbones]                                 |
       |                                                 |
  backbones/bb*.pdb                                      |
       |                                                 |
[03: design-sequences]                                   |
       |                                                 |
  sequences/*.json                                       |
       |                                                 ▼
[04: cofold-score] ◄─────────────────────────────────────┘
       |
  complexes/*.cif
  scores.json
  manifest.json
       |
[05: report]
       |
  report.html  summary.json  ranked_binders.csv
```

Nodes 01–03–04 run sequentially with node 01 feeding both 02 and 04.

---

## Workflow nodes

### Node 01: target-prep

**Goal:** Prepare the target structure and compute a multiple sequence alignment for the target chain.

**Process:** Fetches the target PDB from RCSB using `target_pdb_id`, extracts the specified chain, remaps hotspot author residue IDs to 1-based sequence indices (required for Boltz2 pocket constraints), and calls the MSA-search NIM to build a target A3M alignment.

**Scientific notes:** The target A3M is used in node 04 to condition Boltz2 with evolutionary context on the target chain, improving co-fold confidence. Hotspot residue remapping is necessary because RFdiffusion uses chain+author strings (`C56`) while Boltz2 uses 1-based sequence indices.

**Outputs:**
- `target.pdb` — target chain ATOM records only
- `chain_seq.txt` — one-letter sequence (Cα order)
- `hotspots.json` — remapped hotspot list with author and sequence indices
- `target_a3m.txt` — MSA alignment from MSA-search NIM
- `prep_report.json` — PDB ID, chain, sequence length, hotspot count, MSA depth

### Node 02: generate-backbones

**Goal:** Generate N de novo binder backbone structures conditioned on the target hotspots.

**Process:** Calls the RFdiffusion NIM N times (one call per backbone, distinct seeds) with the target PDB, a contig string combining the fixed target segment and a generated binder segment of the specified length range, and the hotspot residue list. Each call returns one backbone PDB.

**Scientific notes:** RFdiffusion diffuses protein backbones by reversing a noising process conditioned on the target structure and hotspot geometry. The contig format `E1-200/0 60-90` keeps target chain residues 1–200 fixed and generates a 60–90 residue binder chain de novo.

**Outputs:**
- `backbones/bb000.pdb … bb049.pdb` — one RFdiffusion output PDB per backbone
- `backbone_list.json` — index of all backbones with seed and path

### Node 03: design-sequences

**Goal:** Design amino acid sequences for each backbone using ProteinMPNN.

**Process:** Calls the ProteinMPNN NIM once per backbone, selecting only the binder chain for redesign (target chain is fixed). Returns k sequences per backbone with per-sequence NLL scores. The native/WT row is dropped from the FASTA output; only designed sequences are retained.

**Scientific notes:** ProteinMPNN uses a graph neural network to design sequences that are likely to fold into the given backbone. Lower NLL indicates better sequence–backbone compatibility. The binder chain ID is re-read from each backbone PDB — RFdiffusion may rename chains.

**Outputs:**
- `sequences/bb000.json … bb049.json` — designed sequences and NLL scores per backbone
- `sequence_manifest.json` — all candidates with backbone ID, sequence, NLL

### Node 04: cofold-score

**Goal:** Co-fold the top shortlisted binder–target complexes with Boltz2 and score each design.

**Process:** Shortlists the top `cofold_shortlist_n` candidates by ProteinMPNN NLL. For each, calls Boltz2 with a two-chain protein complex payload (binder as single-sequence with no MSA; target with A3M from node 01). Extracts interface confidence (ipTM from `confidence_scores`), binder pLDDT (mean B-factor over binder chain Cα from the returned mmCIF), and self-consistency RMSD (Kabsch Cα RMSD between Boltz2-predicted binder and the original RFdiffusion backbone).

**Scientific notes:** Binder chain is intentionally given no MSA or template — it is de novo with no known homologs. Target chain receives the A3M to improve its fold confidence without biasing the binder prediction. Self-consistency RMSD < 2.0 Å indicates the designed sequence is predicted to fold back into the backbone it was designed for (Bennett et al. 2023 filter).

**Outputs:**
- `complexes/bb000_seq00.cif … ` — Boltz2 mmCIF co-fold structures
- `scores.json` — per-candidate ipTM, binder_plddt, self_consistency_rmsd, proteinmpnn_nll
- `manifest.json` — full campaign manifest (lineage, scores, artifact paths, filter status)

### Node 05: report

**Goal:** Filter, rank, and report all designs.

**Process:** Loads the manifest, applies the three-metric filter (ipTM, pLDDT, RMSD), generates scrambled negative controls at report time (same composition, shuffled sequence), computes success rate, and writes a self-contained HTML report with a ranked table and score distributions.

**Scientific notes:** Success rate (n_passed / n_candidates) is the primary campaign metric — not individual top scores, which are sensitive to pipeline configuration. Scrambled controls validate score distribution separation: a working pipeline should show designed binders scoring clearly above scrambled sequences.

**Outputs:**
- `report.html` — self-contained HTML: ranked table, score distributions, success rate, controls comparison
- `summary.json` — `{n_candidates, n_passed, success_rate, top_iptm, top_sequence}`
- `ranked_binders.csv` — every candidate (pass and fail) with all scores and artifact paths

---

## Parameters

### `n_backbones`

| Value | Description |
|-------|-------------|
| `10–20` (fast test) | Quick sanity check; low diversity |
| `50` (default) | Balanced screening run |
| `100–200` | Production campaign; higher diversity |

**Trade-off:** More backbones = more sequence diversity and higher chance of finding a passer, but linear cost in RFdiffusion and ProteinMPNN NIM calls.

### `cofold_shortlist_n`

| Value | Description |
|-------|-------------|
| `20` | Fast; validates top designs only |
| `50` (default) | Covers ~1 design per backbone |
| `n_backbones × seqs_per_backbone` | Full pool; expensive, risks rate limits |

**Trade-off:** Co-folding is the most expensive stage (~20–30 s per Boltz2 call on hosted). The shortlist caps cost while focusing Boltz2 on the most sequence-compatible designs.

### `nim_mode`

| Value | Description |
|-------|-------------|
| `hosted` (default) | NVIDIA-managed API; requires `NVIDIA_API_KEY`; rate-limited |
| `local` | Self-hosted Docker NIMs; no auth header; no rate limit |

---

## Outputs and interpretation

### `ranked_binders.csv`

Each row is one candidate. Key columns:

| Column | Meaning | Pass threshold |
|--------|---------|---------------|
| `iptm` | Interface confidence (Boltz2 `confidence_scores[0]`) | ≥ 0.8 |
| `binder_plddt` | Mean pLDDT over binder chain Cα | ≥ 80 |
| `self_consistency_rmsd` | Kabsch Cα RMSD vs RFdiffusion backbone (Å) | ≤ 2.0 |
| `proteinmpnn_nll` | Sequence–backbone NLL (lower = better) | — |
| `passed_filter` | True if all three thresholds met | — |

### `report.html`

Self-contained HTML with ranked design table, score distribution histograms, and scrambled control comparison. Open in any browser — no server needed.

---

## Quick start

### Running on Silva

1. Select workflow-033 from the workflow list
2. Set `target_pdb_id`, `target_chain`, and `hotspot_residues` for your target
3. Set `NVIDIA_API_KEY` in platform secrets
4. Click Run

### Test vs production settings

| Setting | Test | Production |
|---------|------|------------|
| `n_backbones` | 10 | 100–200 |
| `seqs_per_backbone` | 4 | 8 |
| `cofold_shortlist_n` | 10 | 50–100 |
| `diffusion_steps` | 20 | 50 |

Default parameters run in ~30–60 min on hosted NIMs. Test settings run in ~5 min.

---

## References

- Watson J.L. et al. "De novo design of protein structure and function with RFdiffusion." *Nature* 620:1089–1100, 2023.
- Dauparas J. et al. "Robust deep learning–based protein sequence design using ProteinMPNN." *Science* 378:49–56, 2022.
- Wohlwend J. et al. "Boltz-2: Towards Accessible and Scalable Joint Structure and Affinity Prediction." *bioRxiv*, 2025.
- Bennett N.R. et al. "Improving de novo protein binder design with deep learning." *Nat Commun* 14:2625, 2023.
- [NVIDIA BioNeMo Agent Toolkit](https://github.com/NVIDIA-BioNeMo/bionemo-agent-toolkit)
