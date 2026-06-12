---
doc_id: structure-prediction
domain: structure-prediction
doc_type: workflow
version: "0.1.0"
deprecated: true
description: >
  Legacy Boltz-2 structure prediction workflow with two case studies
  (Spike RBD and FMC63 antibody). Superseded by workflow-012.
tags: [structure-prediction, boltz-2, legacy, spike-rbd, antibody]
---

# Structure Prediction: Boltz-2 Case Studies (Legacy)

Legacy standalone workflow demonstrating Boltz-2 protein structure prediction on two Nobel Prize-connected case studies: the SARS-CoV-2 Spike RBD and the FMC63 bispecific antibody construct. This workflow has been **superseded by workflow-012**, which provides proper Silva orchestration, input validation, and automated dashboard generation.

## Overview

This workflow runs Boltz-2 structure prediction via manual bash scripts and JSON configuration files, predating the Chiral workflow engine. It contains two independent prediction cases with pre-configured parameters and actual test results. While functional, it lacks the `.chiral/workflow.toml` orchestration, input validation, and automated reporting of workflow-012 (Wohlwend et al., 2024).

## When to use this workflow

**Use workflow-012 instead for all new structure predictions.** This legacy workflow is preserved for reference and reproducibility of the original case study results. Use it only if you need to reproduce the specific Spike RBD or FMC63 predictions with the original parameters and container image.

Do not use this workflow for new predictions — use workflow-012 (Boltz-2 with Silva integration) or workflow-018 (Boltz-2 vs Chai-1 comparison).

## Architecture and data flow

```text
sequences/*.fasta ──> create_boltz_inputs.py ──> inputs/*.fasta ──> job_script.sh ──> results/
                         (extract RBD /                                  |
                          format headers)                          *.pdb, *.json, *.npz
```

Each case study (1_mRNA, 2_antibody) is an independent prediction run with its own configuration. There is no automated orchestration — scripts are run manually.

## Input requirements

- **Format:** FASTA files with Boltz-2 format headers (`>chain_id|protein|`).
- **Preprocessing:** The `create_boltz_inputs.py` script converts raw UniProt FASTA to Boltz-2 format, including RBD region extraction (residues 330–524) for the Spike protein.
- **Placement:** Input FASTA files go in `<case>/inputs/`.

## Case studies

### Case 1: SARS-CoV-2 Spike RBD (1_mRNA)

**Goal:** Predict the 3D structure of the Spike protein Receptor-Binding Domain.

**Process:** Extracts the RBD region (194 residues, positions 331–524) from the full SARS-CoV-2 Spike protein (UniProt P0DTC2, 1273 residues), converts to Boltz-2 format, and runs prediction with 10 diffusion samples and 5 recycling steps.

**Scientific notes:** The Spike RBD is the key domain that binds human ACE2 receptor, making it a primary target for COVID-19 vaccines and therapeutics. The 2023 Nobel Prize in Physiology or Medicine was awarded to Karikó and Weissman for nucleoside base modifications that enabled effective mRNA vaccines — including COVID-19 vaccines targeting this domain.

**Parameters:**
- `recycling_steps`: 5
- `diffusion_samples`: 10
- `sampling_steps`: 200
- `step_scale`: 1.638
- `use_msa_server`: true

**Test results:** 10 PDB models generated in ~4 minutes on GPU (50 output files total: 10 PDB + 10 confidence JSON + 30 NPZ matrices).

### Case 2: FMC63 Bispecific Antibody (2_antibody)

**Goal:** Predict the 3D structure of the FMC63 CAR-T antibody construct.

**Process:** Takes the FMC63-28Z sequence (GenBank ADM64594.1, 489 residues from the full 767 aa construct), converts to Boltz-2 format, and runs prediction with 15 diffusion samples and 7 recycling steps (more than the Spike RBD due to larger size and multi-domain architecture).

**Scientific notes:** FMC63 is a single-chain variable fragment (scFv) that targets CD19, used in FDA-approved CAR-T cell therapies (Kymriah, Yescarta) for B-cell malignancies. The 2018 Nobel Prize in Physiology or Medicine recognized Allison and Honjo for immune checkpoint inhibition, a related breakthrough in cancer immunotherapy.

**Parameters:**
- `recycling_steps`: 7
- `diffusion_samples`: 15
- `sampling_steps`: 200
- `step_scale`: 1.638
- `use_msa_server`: true

**Test results:** 15 PDB models generated (60 output files total: 15 PDB + 15 confidence JSON + 45 NPZ matrices).

## Parameters

### recycling_steps

- **Type:** integer
- **Default:** 5 (Spike RBD) / 7 (FMC63)
- **Description:** Refinement iterations through the model trunk.
- **Guidance:** Higher values for larger, multi-domain proteins. 5 is sufficient for single-domain proteins under 200 residues.

### diffusion_samples

- **Type:** integer
- **Default:** 10 (Spike RBD) / 15 (FMC63)
- **Description:** Number of independent structural models to generate.
- **Guidance:** More samples for larger or flexible targets. 10 is a good starting point.

### sampling_steps

- **Type:** integer
- **Default:** 200
- **Description:** Number of diffusion denoising steps per sample.
- **Guidance:** 200 is the standard setting. Reducing to 100 approximately halves runtime with some quality loss.

### step_scale

- **Type:** float
- **Default:** 1.638
- **Description:** Temperature-like parameter controlling sampling diversity.
- **Guidance:** Default value from Boltz-1. Lower values produce more diverse conformations; higher values produce more conservative predictions. Recommended range: 1.0–2.0.

## Outputs and interpretation

### PDB models

Predicted 3D structures in PDB format, named `<name>_model_N.pdb`. Multiple models from the same prediction explore different conformational samples. The model with the highest confidence score is typically the best prediction.

### Confidence metrics

Each model has a corresponding `confidence_<name>.json` with:
- **confidence_score:** 0.8 × pLDDT + 0.2 × iPTM (range 0–1, higher is better)
- **pTM:** Overall fold quality (> 0.8 = high confidence)
- **iPTM:** Interface quality for complexes (> 0.7 = reliable interface)
- **complex_plddt:** Average per-residue confidence (> 0.7 = confident)
- **complex_pde:** Predicted distance error in Å (lower is better)

### Quality matrices (NPZ)

- **pae_*.npz:** Predicted Aligned Error matrix — low values indicate confident relative positioning between residue pairs
- **pde_*.npz:** Predicted Distance Error matrix
- **plddt_*.npz:** Per-residue pLDDT values

## Quick start

### Running with Docker

```bash
docker run --gpus all -v $(pwd)/1_mRNA:/workspace ghcr.io/chiral-data/boltz_dok_nvidia_2 bash inputs/job_script.sh
```

**Note:** For new predictions, use workflow-012 with Silva instead.

### Container

| Image | Description |
|-------|-------------|
| `ghcr.io/chiral-data/boltz_dok_nvidia_2` | Legacy container with Boltz-2 and GPU support |

## References

- Wohlwend J, Corso G, Passaro S et al. "Boltz-1: Democratizing Biomolecular Interaction Modeling." *bioRxiv*, 2024. DOI: https://doi.org/10.1101/2024.11.19.624167
- [Boltz source code](https://github.com/jwohlwend/boltz)
