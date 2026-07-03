---
doc_id: workflow-029
domain: protein-binder-design
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  De novo protein binder design pipeline: generates binder backbone structures
  with RFdiffusion, designs sequences with ProteinMPNN (SolubleMPNN), folds and
  filters with ColabFold (AlphaFold2 multimer), and ranks by predicted binding
  free energy with PRODIGY.
tags: [binder-design, rfdiffusion, proteinmpnn, colabfold, prodigy, protein-protein-interaction]
---

# Workflow 029: RFdiffusion → ProteinMPNN → PRODIGY Binder Design

De novo protein binder design pipeline. Starting from a receptor PDB, it generates binder backbone structures with RFdiffusion, designs amino acid sequences with ProteinMPNN (SolubleMPNN mode), folds each complex with ColabFold, filters by four structural self-consistency metrics, and ranks passing designs by predicted binding free energy with PRODIGY. The benchmark target is IL-7Rα (PDB 2B5I) from Watson et al. 2023, which yielded 32 confirmed binders in 95 tested, with the best Kd of 40 nM.

## Overview

The pipeline implements the standard computational binder design protocol described in Watson et al. 2023 (*Nature*). RFdiffusion samples binder backbone conformations conditioned on the receptor surface, optionally guided by hotspot residues. ProteinMPNN (SolubleMPNN model) then designs sequences for the binder chain while holding the receptor chain fixed; the soluble model was selected because standard ProteinMPNN yields 0% wet-lab expression for designed binders while SolubleMPNN yields 93.1% (Adaptyv Bio PO102 benchmarks). ColabFold folds each binder:receptor complex as an AlphaFold2 multimer and computes confidence scores. Designs passing all four structural filters proceed to PRODIGY for interface ΔG and Kd estimation.

PRODIGY accuracy on natural crystal structures is r = 0.73, RMSE = 1.89 kcal/mol (Vangone & Bonvin 2016). Performance on AlphaFold2-predicted structures is lower and unbenchmarked; treat ΔG values as relative rankings within a campaign, not absolute affinity predictions.

## When to use this workflow

Use this workflow when you have a receptor PDB (single chain or multi-chain with a known receptor chain) and want to computationally design de novo protein binders against it. It is appropriate when the target has a defined surface patch to design against (optional hotspot residues improve targeting), when you want an end-to-end pipeline from structure to ranked candidate list, and when you intend to validate top hits experimentally by SPR, ITC, or co-crystallisation.

Do not use this workflow for designing binders to intrinsically disordered targets, membrane proteins where the binding interface is not structurally defined, or if you want to redesign an existing binder sequence (use ProteinMPNN alone for that). For enzyme active-site design or symmetric oligomer design, use RFdiffusion2 instead — the original RFdiffusion used here is specifically optimised for PPI binder design. GPU is required for nodes 02, 03, and 04; CPU-only runs are not supported.

## Architecture and data flow

```text
input_files/target.pdb
        |
[01: Validate Target] ── validated_target.pdb ──────────────────────────┐
        |                                                                │
[02: RFdiffusion]                                                        │
        | backbones/design_0.pdb, design_1.pdb, ...                     │
[03: ProteinMPNN]                                                        │
        | sequences/design_000.fasta, design_001.fasta, ...             │
        | (each FASTA has num_seq_per_target records)                    │
[04: ColabFold Fold & Filter] <──────────────────────────────────────────┘
        | folded/design_*_rank001.pdb
        | folded/filter_report.json
[05: PRODIGY Score]
        | results/prodigy_all_designs.json
        | results/ranked_designs.csv
        | results/top10/*.pdb
[06: Report]
        | report/index.html
```

Nodes 04 and 06 each depend on two upstream nodes. Node 04 takes sequences from 03 and the validated receptor from 01. Node 06 takes PRODIGY results from 05 and filter metrics from 04.

**Chain convention (critical):** Node 01 always writes the receptor as chain A regardless of the input chain ID. RFdiffusion then outputs complexes where **chain A = binder (poly-Gly backbone), chain B = receptor**. This convention is preserved in all downstream nodes: ProteinMPNN designs chain A, ColabFold folds the binder:receptor pair, and PRODIGY `--selection A B` means binder=A, receptor=B.

## Input requirements

Place one PDB file in `input_files/`:

- **`target.pdb`** — Receptor structure in PDB format. Must contain at least one protein chain. The chain specified by `receptor_chain` (default `A`) is extracted and cleaned. All other chains are discarded. Node 01 renumbers residues and strips heteroatoms by default.

The IL-7Rα test target (PDB 2B5I) is included. To download a different target:
```bash
wget https://files.rcsb.org/download/XXXX.pdb -O input_files/target.pdb
```

Recommended structure quality: resolution ≤ 3.0 Å, no large gaps in the binding interface, no non-standard residues unless `strip_heteroatoms` is set to `true` (default).

## Workflow nodes

### Node 01: Validate Target

**Goal:** Clean the input receptor PDB and normalize it to a canonical format for downstream tools.

**Process:** Reads `target.pdb` with BioPython, selects the chain specified by `receptor_chain`, and optionally strips heteroatoms. Validates any hotspot residues against the original PDB numbering before renumbering (since users provide hotspots based on the source structure), then renumbers residues starting from 1 using a detach-renumber-add pattern to avoid ID collisions. Writes the output always as chain A regardless of the input chain identity.

**Scientific notes:** Chain normalization to A is required because RFdiffusion hardcodes the receptor as chain A in its contig specification. Stripping heteroatoms removes ligands and waters that can interfere with RFdiffusion's backbone sampling. Residue renumbering prevents gaps in residue numbering that can confuse downstream tools.

**Outputs:**
- `validated_target.pdb` — Cleaned receptor PDB, chain A, residues numbered from 1.

---

### Node 02: RFdiffusion Backbone Design

**Goal:** Generate binder backbone structures conditioned on the receptor surface.

**Process:** Counts receptor residues in chain A of the validated PDB, constructs a Hydra-format contig string `[A1-{receptor_len}/0 {binder_length}-{binder_length}]`, and runs `scripts/run_inference.py` in the RFdiffusion container. Optionally appends `ppi.hotspot_res` to guide backbone sampling toward specified interface residues. Outputs one PDB per design.

**Scientific notes:** The contig format specifies receptor chain A residues 1 through N (held fixed) plus a newly generated binder of exactly `binder_length` residues (chain A in the output, receptor becomes chain B). Hotspot residues from the Watson et al. 2023 IL-7Rα paper (A55, A58, A102, A105, A106) are known to improve binder success rates on this target. Without hotspots, RFdiffusion samples the full receptor surface.

**Outputs:**
- `backbones/design_*.pdb` — Complex PDBs: chain A = binder (poly-Gly), chain B = receptor.

---

### Node 03: ProteinMPNN Sequence Design

**Goal:** Design amino acid sequences for each binder backbone.

**Process:** For each backbone PDB, runs ProteinMPNN with `--pdb_path_chains A` to design chain A (binder) while holding chain B (receptor) fixed. Uses the SolubleMPNN model when `use_soluble_model` is `true`. Renames outputs to zero-padded IDs (`design_000.fasta`, `design_001.fasta`, ...) sorted by backbone integer index. Each FASTA file contains `num_seq_per_target` sequence records.

**Scientific notes:** SolubleMPNN adds a solubility-biasing term that dramatically improves wet-lab expression rates compared to the standard model. Sampling temperature of 0.1 (default) produces conservative, high-probability sequences; higher temperatures (0.3–0.5) increase diversity at the cost of average quality. The zero-padded renaming establishes the canonical design ID used across all downstream nodes.

**Outputs:**
- `sequences/design_*.fasta` — One FASTA per backbone, each containing `num_seq_per_target` sequence records.

---

### Node 04: ColabFold Fold and Filter

**Goal:** Fold each designed sequence as a binder:receptor complex and apply four structural self-consistency filters.

**Process:** Iterates over every sequence record in every FASTA file (creating one ColabFold job per sequence, not per file). For each sequence, concatenates the binder sequence and the receptor sequence extracted from `validated_target.pdb` with a `:` separator to form the ColabFold multimer input. Runs `colabfold_batch` with AlphaFold2 multimer v3. Extracts metrics from `*_scores_rank_001_*.json`: iPTM directly; pLDDT for the binder as `mean(plddt[:binder_length])`; interface PAE as the mean of the off-diagonal cross-chain block of the PAE matrix. Computes Cα RMSD by aligning folded chain A against the backbone PDB by sequence position (not residue number), warning and skipping if the length difference exceeds 5 residues.

**Scientific notes:** The four filter thresholds follow published standards: iPTM ≥ 0.6 and PAE < 10 Å (Zhang et al. 2025; BinderFlow 2025), pLDDT > 80 (Bennett et al. 2023), Cα RMSD ≤ 1.5 Å (Bennett et al. 2023 scRMSD filter). Designs that fold to a structure matching the RFdiffusion backbone (low RMSD) and show confident interface prediction (high iPTM, low PAE) are the most promising experimental candidates.

**Outputs:**
- `folded/design_*_rank001.pdb` — Folded complex PDBs for passing designs only.
- `folded/filter_report.json` — Per-design records with all four metrics, pass/fail status, and fail reason.

---

### Node 05: PRODIGY Affinity Scoring

**Goal:** Score all filter-passing designs with PRODIGY and rank by predicted binding free energy.

**Process:** Loads `filter_report.json` and collects passing design IDs. If none passed, writes empty output files and exits 0 with guidance on relaxing filter thresholds. For each passing PDB, runs `prodigy --selection A B -q` (A = binder, B = receptor) and parses tab-separated stdout for ΔG and Kd. Sorts by ΔG ascending, flags designs above `dg_cutoff` as weak binders, and exports the top N PDBs with optional PyMOL `.pml` scripts using relative paths.

**Scientific notes:** PRODIGY is a fast structure-based predictor trained on X-ray crystal structures (Vangone & Bonvin 2016). Its accuracy degrades on predicted structures; use ΔG rankings as a relative guide, not absolute affinity predictions. Designs with ΔG < −8 kcal/mol are considered strong binder candidates by convention; experimentally confirmed binders typically range from −7 to −15 kcal/mol.

**Outputs:**
- `results/prodigy_all_designs.json` — Full scored and ranked list with all metrics.
- `results/ranked_designs.csv` — CSV version for spreadsheet analysis.
- `results/top10/*.pdb` — Top N folded complex PDB files.
- `results/top10/*.pml` — (optional) PyMOL scripts for interface visualisation.

---

### Node 06: Campaign Report

**Goal:** Generate a self-contained interactive HTML dashboard summarising the campaign.

**Process:** Reads `filter_report.json`, `prodigy_all_designs.json`, and any PDB files in `results/top10/` (glob — no assumed count). Renders a Plotly funnel chart, ΔG ranking bar chart, iPTM vs ΔG scatter plot, Molstar 3D viewer (pinned to v5.9.0) for top structures, a full ranking table, and a methodology note with filter thresholds and PRODIGY caveats. If `prodigy_all_designs.json` is empty, renders a prominent banner without crashing.

**Outputs:**
- `report/index.html` — Fully self-contained HTML report (Bootstrap 5.3.3, Plotly 2.35.2, Molstar 5.9.0; all data embedded inline).

## Parameters

### receptor_chain

| Value / Range | Description |
|---|---|
| `"A"` (default) | Extract chain A from the input PDB. For most single-chain deposited structures. |
| Any chain letter | Extract that chain. For multi-chain deposits (e.g., 2B5I has chains A–F). |

**Trade-off:** Only one chain is used as the receptor. If the functional binding surface spans multiple chains, you must provide a pre-assembled structure.

### hotspot_residues

- **Type:** string
- **Default:** `""`
- **Description:** Comma-separated residue IDs on the receptor (after renumbering) to bias backbone sampling, e.g. `"A55,A58,A102"`.
- **Guidance:** Empty means RFdiffusion samples the full receptor surface. Specifying hotspots from a known binding partner dramatically improves success rate when the epitope is defined (Watson et al. 2023 report ~5× improvement on IL-7Rα with the IL-7 contact patch hotspots).

### num_designs

- **Type:** integer
- **Default:** `100`
- **Description:** Number of RFdiffusion backbone structures to generate.
- **Guidance:** 100 is the minimum for a reasonable hit rate. Watson et al. 2023 screened 95–200 designs per target. Use 10 for local testing.

### binder_length

- **Type:** integer
- **Default:** `80`
- **Description:** Length in residues of the designed binder.
- **Guidance:** 60–100 is the validated range for single-domain binders. Shorter designs have higher expression rates but smaller interface area. Watson et al. 2023 used 65–100 for IL-7Rα.

### num_seq_per_target

- **Type:** integer
- **Default:** `8`
- **Description:** Number of ProteinMPNN sequences per backbone.
- **Guidance:** 8 provides good sequence diversity per backbone. Increasing to 16 or 32 finds better sequences at the cost of proportionally more ColabFold compute.

### min_iptm / min_plddt_binder / max_bb_rmsd / max_pae_interaction

| Parameter | Default | Description |
|---|---|---|
| `min_iptm` | `0.6` | Interface pTM threshold (Zhang 2025; BinderFlow 2025 standard) |
| `min_plddt_binder` | `80.0` | Binder chain mean pLDDT (Bennett 2023 standard) |
| `max_bb_rmsd` | `1.5` | Cα RMSD vs backbone Å (Bennett 2023 scRMSD ≤ 2 Å) |
| `max_pae_interaction` | `10.0` | Interface PAE Å (Zhang 2025; BinderFlow 2025 standard) |

**Test vs production:**

| Setting | Default | Test Override | Production |
|---|---|---|---|
| `num_designs` | 100 | 3–10 | 100–200 |
| `num_seq_per_target` | 8 | 8 | 8–16 |
| `num_recycle` | 6 | 6 | 6–12 |
| `msa_mode` | `single_sequence` | `single_sequence` | `mmseqs2_uniref_env` |
| `min_iptm` | 0.6 | 0.05 (relaxed) | 0.6 |
| `min_plddt_binder` | 80.0 | 30.0 (relaxed) | 80.0 |
| `max_bb_rmsd` | 1.5 | 25.0 (relaxed) | 1.5 |
| `max_pae_interaction` | 10.0 | 35.0 (relaxed) | 10.0 |

## Outputs and interpretation

### filter_report.json

Per-design JSON array with fields: `design_id`, `iptm`, `plddt_binder`, `bb_rmsd`, `pae_interaction`, `pass`, `fail_reasons`. Designs with `pass: true` proceed to PRODIGY scoring. The `fail_reasons` field is a list of all metrics that failed (e.g. `["iptm", "plddt_binder"]`); use it to tune filter thresholds.

### prodigy_all_designs.json / ranked_designs.csv

Ranked list sorted by ΔG ascending. The `weak_binder` flag is `true` when ΔG > `dg_cutoff` (default −8 kcal/mol). Typical strong binders have ΔG < −8 kcal/mol; Kd is derived from ΔG by the Gibbs equation. Published PRODIGY accuracy: r = 0.73, RMSE = 1.89 kcal/mol on crystal structures; lower on predicted structures.

### report/index.html

Interactive single-file HTML dashboard. Open in any modern browser. Contains: design funnel chart showing attrition at each step; ΔG ranking bar chart; iPTM vs ΔG scatter to visualise the correlation between structural confidence and predicted affinity; Molstar 3D viewer for the top N complex structures; full sortable ranking table; methodology caveats.

## Quick start

### Running with Docker

```bash
# Node 01: validate receptor
docker run --rm \
  -v $(pwd)/workflows/workflow-029/input_files:/workspace/inputs \
  -v /tmp/node01_out:/workspace/outputs \
  -v $(pwd)/workflows/workflow-029/01_validate_target:/workspace \
  -w /workspace \
  -e PARAM_TARGET_CHAIN=A \
  biopython:2026_06_30 bash run.sh
```

### Running on Silva

1. Place `target.pdb` in `input_files/`
2. Set `SILVA_WORKFLOW_HOME` to the parent directory of `workflows/workflow-029/`
3. Launch Silva and select "RFdiffusion → ProteinMPNN → PRODIGY Binder Design"
4. Adjust parameters if needed (see Parameters section)
5. Click Run — GPU is required for nodes 02, 03, and 04

### Test vs production settings

| Setting | Default | Quick Test | Production |
|---|---|---|---|
| `num_designs` | `100` | `3–10` | `100–200` |
| `num_seq_per_target` | `8` | `8` | `8–16` |
| `min_iptm` | `0.6` | `0.05` (relaxed) | `0.6` |
| `min_plddt_binder` | `80.0` | `30.0` (relaxed) | `80.0` |
| `max_bb_rmsd` | `1.5` | `25.0` (relaxed) | `1.5` |
| `max_pae_interaction` | `10.0` | `35.0` (relaxed) | `10.0` |
| `msa_mode` | `single_sequence` | `single_sequence` | `mmseqs2_uniref_env` |
| `hotspot_residues` | `""` | `""` | `"A55,A58,A102,A105,A106"` (for 2B5I) |

The included test input (2B5I, IL-7Rα receptor chain A) with `num_designs=3` and relaxed filter thresholds completes in approximately 20–30 minutes on a single GPU. With few backbones, designs are unlikely to pass production thresholds — relaxed thresholds allow all designs through for end-to-end pipeline validation. A production run with `num_designs=100` takes approximately 12–18 hours (bottleneck: ColabFold, ~700 predictions).

## Troubleshooting

**No designs pass the fold filter**
Relax thresholds: increase `max_bb_rmsd` to 2.0, decrease `min_iptm` to 0.5. Also increase `num_designs` — a larger backbone pool improves the odds of finding self-consistent designs.

**PRODIGY exits with non-zero code**
Ensure the input PDB contains both chain A and chain B. ColabFold multimer may occasionally output single-chain models for very short binders — check the folded PDB before running PRODIGY.

**ColabFold OOM on GPU**
Reduce `num_seq_per_target` to 2–4 to lower peak GPU memory per batch.

## References

- Watson J.L. et al. "De novo design of protein structure and function with RFdiffusion." *Nature* 620:1089–1100, 2023. DOI: 10.1038/s41586-023-06415-8
- Bennett N.R. et al. "Atomically accurate de novo design of single-domain antibodies." *Nat Commun* 2023. DOI: 10.1038/s41467-023-40000-1
- Zhang J. et al. "De novo design of protein binders targeting β-stranded peptides." *Nat Commun* 2025. DOI: 10.1038/s41467-025-58095-5
- Vangone A. & Bonvin A.M.J.J. "Contacts-based prediction of binding affinity in protein–protein complexes." *eLife* 5:e07454, 2016. DOI: 10.7554/eLife.07454
- Dauparas J. et al. "Robust deep learning–based protein sequence design using ProteinMPNN." *Science* 378:49–56, 2022. DOI: 10.1126/science.add2187
- [RFdiffusion GitHub](https://github.com/RosettaCommons/RFdiffusion)
- [ProteinMPNN GitHub](https://github.com/dauparasj/ProteinMPNN)
- [ColabFold documentation](https://github.com/sokrypton/ColabFold)
- [PRODIGY documentation](https://github.com/haddocking/prodigy)
