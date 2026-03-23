# Workflow 012: Boltz-2 Structure Prediction

AI-based 3D protein structure prediction using [Boltz-2](https://github.com/jwohlwend/boltz).

## Overview

This workflow takes a protein sequence in Boltz-2 YAML format and predicts its 3D structure using diffusion-based modeling with GPU acceleration. It generates multiple structural models ranked by confidence and an interactive HTML dashboard.

## Nodes

| Node | Name | Description |
|------|------|-------------|
| 01 | Sequence Upload | Validate input YAML with protein sequence(s) |
| 02 | Structure Prediction | Run Boltz-2 prediction (GPU) |
| 03 | Report | Generate interactive HTML dashboard with quality metrics |

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `diffusion_samples` | integer | 10 | Number of structure models to generate |
| `recycling_steps` | integer | 5 | Number of recycling steps |
| `use_msa_server` | boolean | true | Use MSA server for evolutionary information |

## Input Format

Boltz-2 YAML format (see `input_files/prot.yaml`):

```yaml
version: 1
sequences:
  - protein:
      id: A
      sequence: MVLSPADKTNVKAAWGKVGA...
```

## Outputs

- **PDB files**: Predicted 3D structures (`*_model_N.pdb`)
- **Confidence metrics**: JSON files with confidence, pLDDT, PTM, iPTM scores
- **Quality matrices**: PAE, PDE, pLDDT as NPZ files
- **Dashboard**: Interactive HTML report with quality assessment

## Containers

| Node | Image |
|------|-------|
| 01, 02 | `chiral.sakuracr.jp/boltz:2025_09_05` |
| 03 | `chiral.sakuracr.jp/boltz_report:2026_02_13` |

## Example

The test input (`input_files/prot.yaml`) contains a single protein chain of 141 residues (hemoglobin alpha subunit).
