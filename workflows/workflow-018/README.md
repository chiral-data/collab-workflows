# Workflow 018: Boltz-2 vs Chai-1 Comparison

Side-by-side comparison of [Boltz-2](https://github.com/jwohlwend/boltz) and [Chai-1](https://github.com/chaidiscovery/chai-lab) protein structure predictions from a single input sequence.

## Overview

This workflow takes a protein sequence in Boltz-2 YAML format, runs both Boltz-2 and Chai-1 predictions in parallel, and produces an interactive HTML report comparing confidence metrics (pLDDT, pTM, iPTM) and predicted structures side by side.

## Nodes

| Node | Name | Description |
|------|------|-------------|
| 01 | Sequence Upload | Validate input YAML and copy to workspace |
| 02 | Preprocessing | Convert YAML to Boltz-2 and Chai-1 input formats |
| 03 | Boltz-2 Prediction | Run Boltz-2 structure prediction (GPU) |
| 04 | Chai-1 Prediction | Run Chai-1 structure prediction (GPU) |
| 05 | Visualization | Generate interactive HTML comparison report |

Nodes 03 and 04 run in parallel after Node 02 completes.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `diffusion_samples` | integer | 2 | Number of diffusion samples for Boltz-2 |
| `recycling_steps` | integer | 3 | Recycling steps for Boltz-2 structure refinement |
| `use_msa_server` | boolean | true | Use ColabFold MSA server for Boltz-2 |
| `num_trunk_recycles` | integer | 3 | Trunk recycling steps for Chai-1 |
| `num_diffusion_timesteps` | integer | 200 | Diffusion timesteps for Chai-1 (50 for testing) |

## Input Format

Boltz-2 YAML sequence file:

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: MVLSPADKTNVKAAWGK...
```

## Outputs

- **PDB files**: Predicted structures from Boltz-2
- **CIF files**: Predicted structures from Chai-1
- **Report**: `comparison_report.html` — interactive dashboard with confidence metrics, per-residue pLDDT plots, and PAE heatmaps

## Example

From `input/`:
- `prot.yaml` — Human hemoglobin alpha chain (142 residues, UniProt P69905)

## Containers

| Node | Image |
|------|-------|
| 01, 02, 03, 05 | `ghcr.io/chiral-data/boltz:2025_09_05` |
| 04 | `ghcr.io/chiral-data/chai:2026_03_19` |
