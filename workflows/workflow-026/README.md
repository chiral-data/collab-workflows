# Workflow 024: ColabFold Structure Prediction

AlphaFold2-based 3D protein structure prediction using [ColabFold](https://github.com/sokrypton/ColabFold) with MSA from MMseqs2.

## Overview

This workflow takes a protein sequence in FASTA format, runs multiple sequence alignment via the MMseqs2 server, performs AlphaFold2 structure prediction with GPU acceleration, analyzes prediction confidence, and generates visualizations.

## Nodes

| Node | Name | Description |
|------|------|-------------|
| 01 | FASTA Validation | Validate input FASTA: character set, length bounds, monomer/multimer detection |
| 02 | Structure Prediction | MSA + AlphaFold2 inference via `colabfold_batch` (GPU) |
| 03 | Confidence Analysis | Parse score JSONs, compute pLDDT/PAE metrics, emit `confidence_summary.json` |
| 04 | Structure Visualization | pLDDT-colored structure PNG and PAE heatmap |

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `min_length` | integer | 10 | Minimum sequence length (aa) |
| `max_length` | integer | 2500 | Maximum sequence length (aa) |
| `num_models` | integer | 5 | Number of AlphaFold2 models to run (1–5) |
| `num_recycle` | integer | 3 | Recycling iterations (use 6–12 for multimers) |
| `msa_mode` | string | `mmseqs2_uniref_env` | MSA mode |
| `pair_mode` | string | `auto` | Pairing mode: `auto`, `unpaired`, `unpaired_paired` |
| `use_templates` | boolean | false | Enable template search |
| `host_url` | string | `https://api.colabfold.com` | MSA server URL (`--host-url`) |
| `server_timeout_seconds` | integer | 300 | MSA server timeout |
| `top_n_models` | integer | 5 | Models to include in confidence summary |
| `color_by` | string | `pLDDT` | Structure color scheme |

## Input Format

Standard FASTA with one or more sequences:

```
>protein_name
MVLSPADKTNVKAAWGKVGA...
```

For multimers, provide multiple FASTA records or use `:` in the header.

## Outputs

| File | Description |
|------|-------------|
| `validated_sequences.fasta` | Cleaned, validated input sequences |
| `*_unrelaxed_rank_*.pdb` | Ranked predicted structures |
| `*_scores_rank_*.json` | Per-residue pLDDT, PAE matrix, PTM/iPTM scores |
| `confidence_summary.json` | Aggregated confidence metrics (see schema below) |
| `structure_plddt.png` | Structure colored by pLDDT |
| `pae_matrix.png` | PAE heatmap (with chain boundaries for multimers) |

### `confidence_summary.json` Schema

```json
{
  "colabfold_version": "1.6.1",
  "job_id": "string",
  "mode": "monomer | multimer",
  "models": [
    {
      "rank": 1,
      "model_name": "string",
      "plddt_mean": 0.0,
      "plddt_per_residue": [0.0],
      "pae_mean": 0.0,
      "max_pae": 0.0,
      "pae": [[0.0]],
      "interface_pae_mean": null,
      "ptm": 0.0,
      "iptm": null
    }
  ]
}
```

`interface_pae_mean` and `iptm` are `null` for monomers.

## Containers

| Node | Image |
|------|-------|
| 01 | `python:3.11-slim` |
| 02 | `ghcr.io/sokrypton/colabfold:1.6.1-cuda12` |
| 03 | `python:3.11-slim` |
| 04 | `condaforge/mambaforge:latest` |

## Example

The test input (`input_files/sequences.fasta`) is T4 lysozyme — 164 residues, canonical AlphaFold2 benchmark monomer. Expected: `plddt_mean > 90`, `pae_mean < 5 Å`, TM-score > 0.95 vs PDB 2LZM.

### Dry-run (no GPU required)

To verify the `colabfold_batch` command that would be run without executing it:

```bash
cd workflows/workflow-024/02_predict
cp ../input_files/sequences.fasta ./inputs/validated_sequences.fasta
python script.py --dry-run
```
