---
doc_id: workflow-013
domain: protein-design
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  AI-based protein and peptide binder design using BoltzGen. Takes a target
  protein sequence and design specification, generates binder candidates via
  continuous diffusion, and produces a ranked interactive dashboard.
tags: [boltzgen, binder-design, protein-design, diffusion-model, peptide-design, antibody-design]
---

# Workflow 013: BoltzGen Binder Design

AI-based protein and peptide binder design using [BoltzGen](https://github.com/HannesStark/boltzgen). This workflow takes a target protein sequence with a design specification and uses a continuous diffusion model built on Boltz-2 to generate binder candidates — peptides, antibodies (Fab), or general protein binders — ranked by predicted interface quality with an interactive HTML dashboard.

## Overview

BoltzGen is a generative model for designing proteins and peptides that bind to biomolecular targets. It uses a continuous diffusion process (shared with Boltz-2) with a transformer architecture operating at both atom and token levels to progressively denoise 3D coordinates. Rather than encoding residue identity with discrete labels, BoltzGen uses a geometric encoding based on the 14-atom representation, enabling joint training for both structure prediction and design. The model generates an initial set of candidate designs (`num_designs`), then filters and ranks them to produce a final set (`budget`) of diverse, high-quality binders.

BoltzGen solves the *inverse design* problem — generating new binders from scratch — in contrast to Boltz-2 which solves the *forward* problem of predicting structure and affinity for existing complexes (Stark et al., 2025).

## When to use this workflow

Use this workflow when you want to computationally design a new protein or peptide binder against a target of known sequence. It supports de novo peptide design (specifying length range and binding site residues) and scaffold-based antibody Fab design (providing a scaffold library). The workflow requires GPU acceleration.

Do not use this workflow for structure prediction of existing protein complexes — use workflow-012 (Boltz-2) instead. For protein-protein docking of two known structures, use workflow-016 (DiffDock-PP) or workflow-017 (LightDock). For small-molecule virtual screening against a protein target, use workflow-004 (AutoDock Vina) or workflow-002/003 (Smina).

## Architecture and data flow

```text
[01: Target Upload] ──> [02: Binder Design] ──> [03: Report]
        |                       |                     |
  design_spec.yaml         *.cif              boltzgen_dashboard_*.html
  target_summary.json   aggregate_metrics_
                         analyze.csv
```

Nodes run sequentially: 01 → 02 → 03.

## Input requirements

- **Design specification:** A BoltzGen YAML file defining the target sequence and design parameters. Two design types are supported:
  - **De novo** — specify the design chain's length range and target binding site residues
  - **Scaffold library** — provide a directory of scaffold YAML files (e.g., Fab heavy/light chain frameworks)
- **Target structure:** Optional PDB/CIF file if referenced by the design spec. Not required when the target is defined by sequence in the YAML.
- **Sample data:** A test input (`beetletert.yaml`) is included for a 12–20 residue peptide binder against BeetleTERT (PDB 5CQG). Additional examples are in `examples/`.

### Design spec format

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: <target_sequence>
  - protein:
      id: B
      msa: empty
      design:
        type: de_novo          # or "scaffold_library"
        min_length: 12         # de novo only
        max_length: 20         # de novo only
        path: fab_scaffolds    # scaffold_library only
        binding_site:          # optional
          chain: A
          residues: [145, 146, ...]
```

## Workflow nodes

### Node 01: Target Upload

**Goal:** Validate the design specification and prepare inputs for BoltzGen.

**Process:** Parses the design spec YAML (supporting both `sequences` and `entities` formats), validates that required fields are present, identifies design chains and their design type (de_novo or scaffold_library). Copies the design spec to `design_spec.yaml`, copies the target structure file if provided, and copies any referenced scaffold directories. Writes a `target_summary.json` with entity information and design chain IDs.

**Scientific notes:** The design spec defines which chains are fixed (target) and which are to be designed (binder). For de novo design, BoltzGen samples sequences within the specified length range. For scaffold-based design, the model uses provided framework structures and designs the variable regions (e.g., CDR loops for antibodies).

**Outputs:**
- `design_spec.yaml` — validated design specification
- `target_summary.json` — entity summary with design chain IDs
- Target structure and scaffold files (if provided)

### Node 02: Binder Design

**Goal:** Run BoltzGen to generate and rank binder candidates.

**Process:** Invokes `boltzgen run` with the design spec, protocol, num_designs, budget, and `--devices 1` (single GPU). BoltzGen generates `num_designs` initial candidates via continuous diffusion, then filters and ranks them to produce `budget` final designs. Collects output CIF files (designed structures), the aggregate metrics CSV, and any overview PDF from the BoltzGen output directory.

**Scientific notes:** BoltzGen's pipeline includes diffusion-based structure generation, inverse folding, refolding validation, and affinity prediction. Each design is ranked on every metric independently, then assigned its worst rank across all metrics — this "worst-rank" approach selects designs that perform consistently well across all dimensions rather than excelling in just one. The default `num_designs=50` and `budget=10` are suitable for quick testing; production use recommends `num_designs` of 5,000–60,000 and `budget` of 10–100.

**Outputs:**
- `*.cif` — designed binder structures in mmCIF format
- `aggregate_metrics_analyze.csv` — per-design metrics (iPTM, pLDDT, PTM, bb_rmsd, delta_sasa, sequence)
- `results_overview.pdf` — optional summary PDF

### Node 03: Report

**Goal:** Generate an interactive HTML dashboard for analyzing design results.

**Process:** Reads the aggregate metrics CSV and CIF files. Parses per-design metrics (iPTM, pLDDT, PTM, bb_rmsd, delta_sasa_refolded, designed sequence), assigns quality categories based on configurable thresholds, ranks designs by iPTM score, and generates an interactive HTML dashboard with Plotly charts (ranking bar chart, per-design iPTM and RMSD bars, quality pie chart, iPTM vs RMSD scatter) and a 3Dmol.js 3D structure viewer for the top 5 designs.

**Scientific notes:** The dashboard uses iPTM (interface predicted TM-score) as the primary ranking metric, reflecting predicted quality of the binder-target interface. Quality thresholds are: excellent ≥ 0.9, good ≥ 0.7, moderate ≥ 0.5, poor < 0.5.

**Outputs:**
- `boltzgen_dashboard_{target}_{timestamp}.html` — self-contained interactive dashboard with charts, 3D viewer, and detailed results table

## Parameters

### protocol

- **Type:** enum
- **Default:** `"protein-anything"`
- **Node:** 02
- **Description:** BoltzGen design protocol determining the type of binder to generate.

| Value | Description |
|-------|-------------|
| `protein-anything` (default) | Design proteins to bind proteins or peptides. Includes a design folding validation step. |
| `protein-protein` | Protein-protein specific design. |
| `protein-peptide` | Protein-peptide specific design. |

### num_designs

- **Type:** integer
- **Default:** `50`
- **Node:** 02
- **Description:** Number of initial binder candidates to generate during the diffusion sampling phase.

**Trade-off:** More designs improve the chance of finding high-quality binders but increase GPU runtime proportionally.

**Test vs production:** The default of 50 is for quick testing only. For production use, set to 5,000–60,000.

### budget

- **Type:** integer
- **Default:** `10`
- **Node:** 02
- **Description:** Number of final designs to keep after filtering and ranking. Must be ≤ `num_designs`.
- **Guidance:** 10 is suitable for testing. For production campaigns, use 10–100 depending on how many candidates you want to validate experimentally.

### design_spec

- **Type:** file
- **Default:** `.chiral/test_inputs/beetletert.yaml`
- **Node:** 01
- **Description:** Path to the BoltzGen design specification YAML file.

### target_structure

- **Type:** file
- **Default:** (empty)
- **Node:** 01
- **Description:** Optional target structure file (PDB/CIF). Only needed if the design spec references a structure file rather than defining the target by sequence.

## Outputs and interpretation

### iPTM (interface predicted TM-score)

The primary ranking metric. Measures predicted quality of the binder-target interface.

| Range | Interpretation |
|-------|---------------|
| ≥ 0.9 | Excellent — high-confidence predicted binding |
| 0.7–0.9 | Good — likely binds, consider experimental validation |
| 0.5–0.7 | Moderate — may benefit from additional design rounds |
| < 0.5 | Poor — unlikely to bind as designed |

### pLDDT (predicted local distance difference test)

Per-residue structure confidence, averaged over the binder chain. Values > 0.7 indicate confident structure prediction; values > 0.9 indicate very high confidence.

### bb_rmsd (backbone RMSD)

RMSD between the designed structure and its refolded validation structure. Lower values indicate that the design is self-consistent — it refolds to the intended structure. Values < 2.5 Å are recommended.

### delta_sasa (change in solvent-accessible surface area)

Measures the buried surface area upon binding. Larger values indicate more extensive binding interfaces, which generally correlate with stronger binding affinity.

### PTM (predicted TM-score)

Overall structural quality of the designed complex. Higher values indicate better predicted global fold quality.

## Quick start

### Running with Docker

Nodes 01 and 02 require GPU acceleration:

```bash
docker pull ghcr.io/chiral-data/boltzgen:2026_02_13    # Nodes 01, 02
docker pull ghcr.io/chiral-data/boltz_report:2026_02_13 # Node 03
```

### Running on Silva

1. Select "BoltzGen Binder Design" from the workflow list
2. Upload your design specification YAML (or use the default BeetleTERT test)
3. Set `protocol` appropriate for your design type
4. Adjust `num_designs` and `budget` for your desired quality level
5. Click Run

### Test vs production settings

| Setting | Test (default) | Production |
|---------|---------------|------------|
| `design_spec` | `beetletert.yaml` | Your design spec |
| `protocol` | `protein-anything` | As appropriate |
| `num_designs` | `50` | `5,000`–`60,000` |
| `budget` | `10` | `10`–`100` |

A successful test run with the BeetleTERT defaults generates 10 peptide binder candidates (12–20 residues) and completes in ~15–30 minutes on a single GPU.

## Examples

### Peptide binder against BeetleTERT

`examples/peptide/beetletert.yaml` — de novo design of a 12–20 residue peptide binder targeting residues 145–150, 180–184 of BeetleTERT (PDB 5CQG).

### Fab antibody binder against PD-L1

`examples/antibody/pdl1.yaml` — scaffold-based Fab design against PD-L1 (PDB 7UXQ). Uses the scaffold library in `examples/antibody/fab_scaffolds/` for heavy and light chain frameworks.

## References

- Stark, H. et al. "BoltzGen: Toward Universal Binder Design." *bioRxiv*, 2025. DOI: https://doi.org/10.1101/2025.11.20.689494
- [BoltzGen GitHub](https://github.com/HannesStark/boltzgen)
- [Boltz-2](https://github.com/jwohlwend/boltz)
