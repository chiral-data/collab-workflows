# Workflow 013: BoltzGen Binder Design

AI-based protein/peptide binder design using [BoltzGen](https://github.com/HannesStark/boltzgen).

## Overview

This workflow takes a target protein structure and design specification, then uses BoltzGen to generate binder candidates (antibodies or peptides) with GPU acceleration. It produces ranked designs, an interactive HTML dashboard, and 3D visualizations via Mol*.

## Nodes

| Node | Name | Description |
|------|------|-------------|
| 01 | Target Upload | Validate target structure and design spec |
| 02 | Binder Design | Run BoltzGen design (GPU) |
| 03 | Report | Generate interactive HTML dashboard |
| 04 | Mol* Visualization | Select top designs for 3D visualization |

Nodes 03 and 04 run in parallel after Node 02 completes.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `protocol` | enum | protein-anything | Design protocol |
| `num_designs` | integer | 50 | Number of binder candidates to generate |
| `budget` | integer | 10 | Final designs after filtering |

## Examples

### Antibody (Fab) binder against PD-L1

From `examples/antibody/`:
- `pdl1.yaml` — Design spec for Fab binder vs PD-L1 target (7uxq)
- `fab_scaffolds/` — Scaffold library for heavy/light chain frameworks

### Peptide binder against BeetleTERT

From `examples/peptide/`:
- `beetletert.yaml` — Design spec for 12-20 residue peptide vs BeetleTERT (5cqg)
- De novo design with specified binding site residues

## Input Format

BoltzGen design specification YAML:

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
        type: de_novo
        min_length: 12
        max_length: 20
        binding_site:
          chain: A
          residues: [145, 146, ...]
```

## Outputs

- **CIF files**: Designed binder structures
- **Metrics CSV**: `aggregate_metrics_analyze.csv` with confidence, pLDDT, iPTM, binding energy
- **Dashboard**: Interactive HTML report with quality assessment and plots
- **Mol* files**: Top-ranked designs for 3D viewer

## Containers

| Node | Image |
|------|-------|
| 01, 02, 04 | `ghcr.io/chiral-data/boltzgen:2026_02_13` |
| 03 | `ghcr.io/chiral-data/boltz_report:2026_02_13` |
