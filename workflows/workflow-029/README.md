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

De novo protein binder design pipeline. Starting from a receptor PDB, it generates binder backbone structures with RFdiffusion, designs amino acid sequences with ProteinMPNN (SolubleMPNN mode), folds each complex with ColabFold, filters by four structural self-consistency metrics, and ranks passing designs by predicted binding free energy with PRODIGY. The bundled example target is interleukin-2 (PDB 2B5I, chain A), with hotspots on its IL-2Rα-binding face.

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

**Chain convention (critical):** Node 01 always writes the receptor as chain A regardless of the input chain ID. RFdiffusion, whose contig lists the receptor first, emits that receptor as chain A and the generated binder as chain B; node 02 relabels them before publishing, so every backbone it writes is **chain A = binder (poly-Gly backbone), chain B = receptor**. This convention holds in all downstream nodes: ProteinMPNN designs chain A, ColabFold folds the binder:receptor pair, and PRODIGY `--selection A B` means binder=A, receptor=B.

## Input requirements

Place one PDB file in `input_files/`:

- **`target.pdb`** — Receptor structure in PDB format. Must contain at least one protein chain. The chain specified by `receptor_chain` (default `A`) is extracted and cleaned. All other chains are discarded. Node 01 strips heteroatoms by default and preserves the deposited residue numbering, which is what lets node 02 describe the receptor's real chain breaks.

The IL-2 test target (PDB 2B5I, chain A) is included. To download a different target:
```bash
wget https://files.rcsb.org/download/XXXX.pdb -O input_files/target.pdb
```

Recommended structure quality: resolution ≤ 3.0 Å, no non-standard residues unless `strip_heteroatoms` is set to `true` (default). Gaps in the deposited coordinates are detected and expressed to RFdiffusion as real index jumps, so a target with disordered loops is usable — but a gap *inside* the surface you mean to target still removes the epitope you wanted.

## What this workflow sends outside your machine

**By default, node 04 sends your binder and receptor sequences to a third-party server.**

`msa_mode` defaults to `mmseqs2_uniref_env`, which builds each multiple sequence alignment by querying the public ColabFold MSA server at `api.colabfold.com`. One query goes out per predicted complex, carrying the designed binder sequence and the receptor sequence extracted from your target. Nothing else leaves the machine: coordinates, scores and reports stay local.

For the bundled IL-2 demo this is harmless. **For a proprietary target it may not be** — a run at the production profile (`num_designs=100`, `num_seq_per_target=8`) is roughly 800 predictions and about 1,600 sequence queries to a service run by someone else.

**To avoid it,** set `msa_mode = "single_sequence"`. Nothing then leaves the machine. The cost is real and is why that is not the default: without co-evolutionary signal AlphaFold2-multimer cannot place two chains relative to each other, so iPTM and interface PAE stop being meaningful and the two interface filters become unreachable. You get a pipeline that runs offline and a funnel that cannot report an interface.

**Why the default is what it is.** `mmseqs2_uniref_env` is what makes the interface metrics interpretable at all, and it is also what gives AlphaFold2 a usable model of the receptor. In `single_sequence` mode the receptor scores around 40 pLDDT, and supplying its own coordinates as a template does not rescue it — measured at 45 for the single segment node 04 folds, because a fragment carved out of a larger fold is not a folding unit on its own. So the MSA is doing work no local substitute was found to replace, and the default stays, disclosed here rather than left as a quiet dependency.

**Throughput, which is a second consequence of the same dependency.** The queries are fast once they start — 16 to 60 seconds each — but submissions queue behind other users of a shared service. In one measured 16-prediction run, 45 of 89 minutes were queue wait, and a single submission stalled for 27 minutes before starting. The delay is intermittent rather than progressive, so it is contention, not a throttle tightening against you. There is no resume: a long run interrupted partway through is lost. Both are reasons to precompute alignments with `colabfold_search` against local databases, or to self-host the API and point at it, before any large campaign.

## Workflow nodes

### Node 01: Validate Target

**Goal:** Clean the input receptor PDB and normalize it to a canonical format for downstream tools.

**Process:** Reads `target.pdb` with BioPython, selects the chain specified by `receptor_chain`, and optionally strips heteroatoms. Resolves any hotspot residues against the input PDB numbering, which with the default `renumber_residues = false` is an identity mapping. Then walks the chain and records its contiguous segments and the breaks between them, measuring the Cα–Cα distance across each one. Writes the output always as chain A regardless of the input chain identity.

A residue is split off into a new segment at every discontinuity in the numbering, because a contig range may only name residues RFdiffusion has indexed. A Cα–Cα distance above **4.5 Å** is additionally recorded as a break: an ideal trans peptide bond spans 3.80 Å, so the threshold leaves headroom for strained geometry without admitting a real break. A break with no accompanying gap in the numbering is reported but cannot be acted on — the index jump a contig encodes is `next - prev - 1`, which is zero in that case. Residues with no Cα atom, and residues carrying an insertion code, are left out of the segments with a warning, matching the residues RFdiffusion is able to address.

All of this travels to node 02 as `REMARK 999` records prepended to the output, alongside the hotspot records:

```
REMARK 999 HOTSPOT A38 A38                  original → resolved, one per hotspot
REMARK 999 HOTSPOTS A38,A42,A45,A68         the summary line node 02 reads
REMARK 999 NUMBERING ORIGINAL               or RENUMBERED
REMARK 999 SEGMENT 1 A 6 75 70              index, chain, first, last, residue count
REMARK 999 SEGMENT 2 A 80 98 19
REMARK 999 SEGMENT 3 A 103 133 31
REMARK 999 SEGMENTS A6-75/A80-98/A103-133   the summary line node 02 reads
REMARK 999 BREAK 75 80 4 7.24 BOTH          prev, next, missing, Cα–Cα, kind
REMARK 999 BREAK 98 103 4 9.28 BOTH
```

`kind` is `NUMBERING`, `DISTANCE` or `BOTH`, naming which criterion fired. The `SEGMENTS` line is written even for a gapless target, so node 02 can tell "node 01 handled this" from "older input". A hotspot absent from the structure is named in a warning and dropped from the set.

**Scientific notes:** Chain normalization to A is required because RFdiffusion hardcodes the receptor as chain A in its contig specification. Stripping heteroatoms removes ligands and waters that can interfere with RFdiffusion's backbone sampling.

The deposited numbering is preserved because the gaps in it are load-bearing. RFdiffusion is given real coordinates either way — the atoms never move — but what tells the model two residues are bonded is the difference between their positional indices. Renumbering sequentially from 1 closes every gap, so the residues flanking a disordered loop arrive as consecutive, which is to say as bonded; on the bundled target that presented pairs 7.2 Å and 9.3 Å apart as peptide bonds, and on chain D a pair 33.9 Å apart. Preserving the numbering lets node 02 name each real segment, and RFdiffusion then advances its index by `next_start - prev_end - 1` between `/`-separated segments of one chain. `renumber_residues = true` restores the old behaviour and warns, naming each break it fuses; it is there for a target whose numbering a contig cannot express at all — non-positive or duplicated residue numbers.

**Outputs:**
- `validated_target.pdb` — Cleaned receptor PDB, chain A, residues numbered as in the input, carrying the `REMARK 999` hotspot, segment and break records above.

---

### Node 02: RFdiffusion Backbone Design

**Goal:** Generate binder backbone structures conditioned on the receptor surface.

**Process:** Reads the receptor's segments from node 01's `REMARK 999 SEGMENTS` line, constructs a Hydra-format contig string `[{segments}/0 {binder_length}-{binder_length}]` — for the bundled target, `[A6-75/A80-98/A103-133/0 80-80]` — and runs `scripts/run_inference.py` in the RFdiffusion container. A target carrying no such line, such as one written by an older node 01, has its segments derived from its own Cα records on gaps in the numbering alone; for a renumbered target that yields the single span `[A1-{receptor_len}/0 ...]` this node used to build unconditionally. Optionally appends `ppi.hotspot_res` to guide backbone sampling toward specified interface residues, taking the IDs from node 01's `REMARK 999 HOTSPOTS` line so they match the coordinates it wrote; a target carrying no such line falls back to the raw parameter. Outputs one PDB per design.

**Scientific notes:** The contig lists each contiguous run of the receptor (held fixed) plus a newly generated binder of exactly `binder_length` residues. The separator between two runs carries the meaning: `/` keeps them in one chain and advances RFdiffusion's positional index by the real size of the gap, whereas `/0` declares a 200-index break — a separate chain, free to drift away. The receptor's segments are therefore joined with `/`, and only the binder is preceded by `/0`.

RFdiffusion assigns output chains in contig order, so it writes the receptor as A and the binder as B; the node relabels them afterwards so the binder is chain A. The receptor keeps its input numbering in the output, gaps included, while the designed chain is numbered from 1.

For the bundled target, `A38,A42,A45,A68` are the IL-2 residues that contact IL-2Rα in this very structure — Arg38, Phe42, Tyr45 and Glu68, the CD25-binding face. They are passed through as given, since node 01 no longer renumbers. Without hotspots, RFdiffusion samples the full receptor surface.

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

**Process:** Iterates over every sequence record in every FASTA file (creating one ColabFold job per sequence, not per file). **This is the node that makes outbound network calls** — one MSA query per prediction under the default `msa_mode`; see the disclosure section above. For each sequence, concatenates the binder sequence and the receptor sequence extracted from `validated_target.pdb` with a `:` separator to form the ColabFold multimer input. When node 01 reports the receptor in more than one segment, the sequence is trimmed to the single segment carrying the hotspots — the largest one if no hotspots were given. Runs `colabfold_batch` with AlphaFold2 multimer v3. Extracts metrics from `*_scores_rank_001_*.json`: iPTM directly; pLDDT for the binder as `mean(plddt[:binder_length])`; interface PAE as the mean of the off-diagonal cross-chain block of the PAE matrix. Computes Cα RMSD by aligning folded chain A against the backbone PDB by sequence position (not residue number), warning and skipping if the length difference exceeds 5 residues.

**Scientific notes:** The receptor is trimmed rather than folded whole because ColabFold has no way to express a gap inside a chain: concatenating the fragments would fold two lobes as though a peptide bond joined them, which is the same fault node 02 avoids with its contig. The cost is real — surface outside the chosen segment is absent, so a binder that contacts two lobes is scored against one. Splitting the receptor across `:` instead would make iPTM a global measure over binder plus *n* receptor chains, which silently changes what the two interface filters mean.

The four filter thresholds follow published standards: iPTM ≥ 0.6 and PAE < 10 Å (Zhang et al. 2025; BinderFlow 2025), pLDDT > 80 (Bennett et al. 2023), Cα RMSD ≤ 1.5 Å (Bennett et al. 2023 scRMSD filter). Designs that fold to a structure matching the RFdiffusion backbone (low RMSD) and show confident interface prediction (high iPTM, low PAE) are the most promising experimental candidates.

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

**What is in the bundled 2B5I.** The deposit is the quaternary IL-2 receptor complex: chain A is interleukin-2, B is IL-2Rβ, C is the common γ chain, and D is IL-2Rα (CD25). The example targets chain A, so the designed binder competes with the receptors for IL-2. Any of the other chains is a legitimate target instead. Each carries disordered loops, which node 01 measures and node 02 expresses to RFdiffusion as real gaps:

| Chain | Residues | Segments | Breaks |
|---|---|---|---|
| A (IL-2) | 120 | `A6-75/A80-98/A103-133` | 4 missing at 7.2 Å, 4 at 9.3 Å |
| B (IL-2Rβ) | 196 | `A6-24/A31-207` | 6 missing at 12.5 Å |
| C (γc) | 191 | `A34-224` | none |
| D (IL-2Rα) | 122 | `A1-47/A53-64/A103-165` | 5 missing at 4.9 Å, 38 at 33.9 Å |

Chain D is usable, but its 38-residue gap leaves the target presenting as two lobes that are 33.9 Å apart; hotspots should sit within one of them, since node 04 folds the receptor one segment at a time.

### hotspot_residues

- **Type:** string
- **Default:** `""`
- **Description:** Comma-separated residue IDs on the receptor to bias backbone sampling, e.g. `"A38,A42,A45"`. Given in the numbering of **your input PDB**, and used as given — node 01 preserves that numbering, so no translation happens unless `renumber_residues` is on. They also decide which segment node 04 folds when the receptor has chain breaks.
- **Guidance:** Empty means RFdiffusion samples the full receptor surface. Specifying hotspots from a known binding partner improves success rate when the epitope is defined; Watson et al. 2023 report roughly a 5× improvement from targeting a known contact patch. Where the deposit contains the natural partner, as this one does, the epitope can be read straight off the structure.

### msa_mode

| Value | Description |
|---|---|
| `mmseqs2_uniref_env` (default) | Build each MSA by querying the public ColabFold server. Required for the interface metrics to mean anything. **Sends your sequences to a third party** — see the disclosure section above. |
| `single_sequence` | No network, no MSA. iPTM and interface PAE are not meaningful and `min_iptm` / `max_pae_interaction` become unreachable. |

**Trade-off:** this is the one parameter that decides whether the workflow is self-contained. The default is not self-contained, deliberately, because the interface half of the funnel is worthless without an MSA. Choose `single_sequence` when the target is confidential and you intend to judge designs on the binder metrics alone — binder pLDDT and Cα RMSD remain valid in that mode.

### num_designs

- **Type:** integer
- **Default:** `100`
- **Description:** Number of RFdiffusion backbone structures to generate.
- **Guidance:** 100 is the minimum for a reasonable hit rate. Watson et al. 2023 screened 95–200 designs per target. Use 10 for local testing.

### binder_length

- **Type:** integer
- **Default:** `80`
- **Description:** Length in residues of the designed binder.
- **Guidance:** 60–100 is the validated range for single-domain binders. Shorter designs have higher expression rates but smaller interface area. Watson et al. 2023 used 65–100.

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
| `msa_mode` | `mmseqs2_uniref_env` | `mmseqs2_uniref_env` | `mmseqs2_uniref_env` |
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
  -v $(pwd)/workflows/workflow-029/input_files:/work/inputs:ro \
  -v /tmp/node01_out:/work/outputs \
  -v $(pwd)/workflows/workflow-029/01_validate_target:/src:ro \
  -w /work \
  -e PARAM_RECEPTOR_CHAIN=A \
  -e PARAM_HOTSPOT_RESIDUES=A38,A42,A45,A68 \
  biopython:2026_06_30 python3 /src/script.py
```

### Running locally with Silva

1. Place `target.pdb` in `input_files/`
2. Set `SILVA_WORKFLOW_HOME` to the parent directory of `workflows/workflow-029/`
3. Launch Silva and select "RFdiffusion → ProteinMPNN → PRODIGY Binder Design"
4. Adjust parameters if needed (see Parameters section)
5. Click Run — GPU is required for nodes 02, 03, and 04

### Running from the browser

Workflow 029 needs a GPU for nodes 02, 03 and 04, and a full production run takes most of a day. Start with a quick test run to confirm the pipeline works end to end against your target, then re-run with production settings.

**1. Pick the workflow.** Select **RFdiffusion → ProteinMPNN → PRODIGY Binder Design** from the workflow list.

**2. Upload your target.** Upload your receptor structure to `input_files/` as `target.pdb`. If you upload nothing, the run uses the bundled IL-2 example target (PDB 2B5I, chain A) — useful for a first end-to-end test, but not for your own campaign. One PDB only; the workflow reads a single receptor.

**3. Set the parameters that depend on your target.**

| Parameter | Set it to |
|---|---|
| `receptor_chain` | The chain ID of the receptor in *your* uploaded PDB (node 01 rewrites it to chain A afterwards) |
| `hotspot_residues` | The surface patch you want the binder to hit, e.g. `A38,A42,A45`. Leave empty and RFdiffusion picks for you — usable, but noticeably less targeted |
| `binder_length` | `80` suits most PPI targets; shorten only if you have a reason |

**4. Choose quick test or production.** Use the Quick Test column of the table below for your first run: `num_designs=3` with relaxed filters lets everything through, so you can see each node produce output. Budget 15 minutes at `num_seq_per_target=2` and around 45 at the column's default of 8 — node 04 dominates either way. Production (`num_designs=100`, default filters) takes about 12–18 hours, nearly all of it in node 04.

**5. Start the run** and follow it node by node in the log view. Nodes 02, 03 and 04 request a GPU automatically and have no CPU fallback, so the run needs a GPU-capable worker.

#### What a healthy run looks like

| Node | What you should see | Share of runtime |
|---|---|---|
| 01 Validate Target | `validated_target.pdb`, carrying the hotspot mapping, the receptor's segments and any chain breaks as `REMARK 999` records; a hotspot that isn't found is warned about and dropped, not fatal | Seconds |
| 02 RFdiffusion | `backbones/` holding a `design_N.pdb` and a `design_N.trb` per `num_designs`, plus a `traj/` folder of diffusion trajectories | ~4 min for 3 designs |
| 03 ProteinMPNN | `sequences/design_NNN.fasta`, one per backbone. Each holds `num_seq_per_target` designed sequences *plus the native input sequence* as its first record | Seconds |
| 04 ColabFold Fold & Filter | `folded/design_NNN_seq_N_rank001.pdb`, one per designed sequence, and `filter_report.json` scoring every design whether it passed or not | The bulk of the run |
| 05 PRODIGY | `results/ranked_designs.csv`, `results/prodigy_all_designs.json`, and `results/top{top_n}/` with a `.pdb` and a PyMOL `.pml` per exported design | Seconds |
| 06 Report | `report/index.html`, self-contained, around 1 MB | Seconds |

`top{top_n}/` is named from `top_n` at runtime, and holds every design that passed when fewer pass than `top_n` — so a quick test that passes 6 designs gives a `top10/` with 6 in it, not 10.

The figures above come from a measured run: the bundled 2B5I target, `num_designs=3`, `num_seq_per_target=2`, relaxed quick-test filters, single GPU (RTX 3060 12 GB). End to end in 15 minutes, of which node 04 took about 11 — six ColabFold predictions. Node 04 scales with `num_designs × num_seq_per_target`, so the Quick Test column's `num_seq_per_target=8` is roughly four times that node again.

#### Collecting results

Download `report/index.html` first — it is self-contained and opens in any browser, with the design funnel, the ranking table, an iPTM vs ΔG scatter and a 3D viewer for the top N. Then take `results/ranked_designs.csv` for the full ranking and `results/top10/*.pdb` for the structures you'd hand to a wet lab.

Read ΔG as a *relative* ranking within one campaign. PRODIGY was validated on crystal structures (r = 0.73, RMSE = 1.89 kcal/mol); its accuracy on predicted structures is unbenchmarked, so the ordering is meaningful but the absolute kcal/mol values are not.

#### If nothing passes the filter

With a small `num_designs` this is expected — that's why the quick-test column relaxes all four thresholds. On a production run, raise `num_designs` before loosening thresholds: a larger backbone pool is the reliable way to find self-consistent designs. See Troubleshooting below.

### Test vs production settings

| Setting | Default | Quick Test | Production |
|---|---|---|---|
| `num_designs` | `100` | `3–10` | `100–200` |
| `num_seq_per_target` | `8` | `8` | `8–16` |
| `min_iptm` | `0.6` | `0.05` (relaxed) | `0.6` |
| `min_plddt_binder` | `80.0` | `30.0` (relaxed) | `80.0` |
| `max_bb_rmsd` | `1.5` | `25.0` (relaxed) | `1.5` |
| `max_pae_interaction` | `10.0` | `35.0` (relaxed) | `10.0` |
| `msa_mode` | `mmseqs2_uniref_env` | `mmseqs2_uniref_env` | `mmseqs2_uniref_env` |
| `hotspot_residues` | `""` | `""` | `"A38,A42,A45,A68"` (for the bundled 2B5I chain A) |

The included test input (2B5I chain A, interleukin-2) with `num_designs=3` and relaxed filter thresholds completes in 15 minutes on a single GPU at `num_seq_per_target=2`, or roughly 45 at the default of 8. With few backbones, designs are unlikely to pass production thresholds — relaxed thresholds allow all designs through for end-to-end pipeline validation. A production run with `num_designs=100` takes approximately 12–18 hours (bottleneck: ColabFold, ~700 predictions).

## Troubleshooting

**No designs pass the fold filter**
First check `msa_mode`. With `single_sequence`, AlphaFold2-multimer has no co-evolutionary signal and cannot place two chains relative to each other: iPTM collapses to ~0.1 and interface PAE sits near the matrix ceiling for every design, including known native complexes. The interface filters are unreachable in that mode. Otherwise, relax thresholds: increase `max_bb_rmsd` to 2.0, decrease `min_iptm` to 0.5. Also increase `num_designs` — a larger backbone pool improves the odds of finding self-consistent designs.

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
