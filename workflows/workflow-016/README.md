---
doc_id: workflow-016
domain: protein-protein-docking
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  Antibody-antigen docking using DiffDock-PP diffusion models, with
  interface analysis and RMSD-based comparison to reference structures.
tags: [antibody, antigen, docking, diffdock-pp, diffusion-model, protein-protein]
---

# Workflow 016: DiffDock-PP Antibody-Antigen Docking

A protein-protein docking pipeline that uses DiffDock-PP, a diffusion generative model for rigid-body protein-protein docking, to predict antibody-antigen binding poses. The workflow splits a co-crystal complex, generates ranked docking poses with confidence scores, analyzes the predicted binding interface, and validates predictions against the reference structure via RMSD metrics.

## Overview

DiffDock-PP applies a score-based diffusion model to the SE(3) space of rigid-body transformations (translations and rotations) to generate protein-protein docking poses. It represents proteins as heterogeneous geometric graphs at residue level, with node features derived from ESM-2 language model embeddings. A separate confidence model ranks the generated poses. The pipeline is designed for re-docking benchmarks and prospective antibody-antigen docking, producing ranked poses, interface contact analysis, and quantitative comparison to experimental structures (Ketata et al., 2023).

## When to use this workflow

Use this workflow when you have an antibody-antigen co-crystal structure (PDB file) and want to benchmark DiffDock-PP's docking accuracy, or when you have separate antibody and antigen structures and want to predict their binding pose. The input must be a PDB file containing protein chains only (standard amino acids). The workflow auto-detects antibody chains (H/L convention) but supports manual chain specification.

Do not use this workflow for small-molecule docking — use workflow-004 (AutoDock Vina) or workflow-003 (Smina) instead. Do not use this workflow for flexible docking where backbone conformational changes are expected — DiffDock-PP performs rigid-body docking only. For general protein-protein docking beyond antibody-antigen systems, the workflow will still function but the auto-detection heuristic assumes antibody H/L chain naming. For protein structure prediction (not docking), use workflow-012 (Boltz-2) or workflow-018 (Boltz-2 vs Chai-1).

## Architecture and data flow

```text
input.pdb ──> [01: Complex Splitting] ──> [02: Structure Prep] ──> [03: Feature Extraction]
                     │                          │                         │
               original_complex.pdb       processed PDBs           ESM-2 embeddings
                     │                          │                         │
                     │                          ├─────────────────────────┘
                     │                          ▼
                     │                    [04: DiffDock-PP Inference]
                     │                          │
                     │                    ranked poses
                     │                     ┌────┴────┐
                     │                     ▼         ▼
                     │              [05: Interface] [06: Docking Comparison] ◄──┘
                     └──────────────────────────────────────┘
```

Nodes 01 → 02 → 03 run sequentially. Node 04 depends on nodes 02 and 03. Nodes 05 and 06 both depend on node 04 and can run in parallel.

## Input requirements

- **Format:** PDB file containing an antibody-antigen co-crystal structure with standard amino acid residues.
- **Constraints:** The antibody chains should follow standard naming (H for heavy chain, L for light chain) for auto-detection. Non-standard residues, water molecules, and heteroatoms are removed during structure preparation.
- **Placement:** Place the PDB file in `input_files/` and set the `input_pdb` parameter to its filename.
- **Sample data:** `input_files/5B8C.pdb` — an antibody-antigen co-crystal structure used as the default test case.

## Workflow nodes

### Node 01: Complex Splitting

**Goal:** Separate the antibody-antigen co-crystal complex into individual protein components.

**Process:** Parses the input PDB using BioPython, auto-detects antibody chains by looking for H+L chain identifiers (falls back to VH+VL naming, then first two chains). Extracts antibody and antigen into separate PDB files. Records chain metadata (residue counts, atom counts, detection method) in `chain_info.json`. The user can override auto-detection by specifying `antibody_chains` explicitly.

**Scientific notes:** Standard antibody nomenclature designates the heavy chain as "H" and the light chain as "L". The Fv (fragment variable) region containing the CDR loops is sufficient for docking — the constant domains (Fc) do not participate in antigen recognition and are typically excluded from docking calculations.

**Outputs:**
- `antibody.pdb` — extracted antibody chains
- `antigen.pdb` — extracted antigen chains
- `original_complex.pdb` — copy of the input complex for later comparison
- `chain_info.json` — chain detection metadata and statistics

### Node 02: Structure Preparation

**Goal:** Clean and standardize protein structures for DiffDock-PP inference.

**Process:** Removes water molecules (HOH, WAT), heteroatoms (ligands, ions), and non-standard amino acid residues. Renumbers residues per chain starting from 1 to ensure consistent indexing. Retains only standard amino acid residues classified by BioPython's `is_aa(standard=True)`.

**Scientific notes:** DiffDock-PP operates on protein backbone geometry and requires clean, standard protein structures. Non-standard residues and solvent molecules from crystallography can introduce artifacts in the geometric graph representation. Renumbering ensures compatibility with the model's residue-level featurization.

**Outputs:**
- `processed_antibody.pdb` — cleaned antibody structure
- `processed_antigen.pdb` — cleaned antigen structure

### Node 03: Feature Extraction

**Goal:** Extract amino acid sequences and prepare ESM-2 language model embeddings for DiffDock-PP.

**Process:** Converts PDB residues to FASTA sequences using standard IUPAC three-to-one letter mapping. Generates ESM-2 embeddings using the `esm2_t33_650M_UR50D` model (33 Transformer layers, 1024-dimensional per-residue embeddings). Outputs separate FASTA and feature tensor files for antibody and antigen.

**Scientific notes:** ESM-2 embeddings capture evolutionary and structural information learned from approximately 250 million protein sequences (Lin et al., 2023). These embeddings provide rich sequence-level context that complements the 3D structural features (alpha-carbon coordinates) used by DiffDock-PP's geometric graph network, improving docking accuracy over sequence-naive approaches.

**Outputs:**
- `antibody.fasta`, `antigen.fasta` — extracted sequences
- `antibody_features.pt`, `antigen_features.pt` — ESM-2 embedding tensors
- `sequence_info.json` — chain lengths and amino acid composition

### Node 04: DiffDock-PP Inference

**Goal:** Generate ranked antibody-antigen docking poses using the DiffDock-PP diffusion model.

**Process:** Runs DiffDock-PP inference using two model checkpoints: a score model (`large_model_dips`) for pose generation via iterative denoising over SE(3) transformations, and a confidence model (`confidence_model_dips`) for pose ranking. Generates `num_samples` candidate poses through `inference_steps` denoising steps, ranks them by confidence score (descending), and exports the top-ranked pose as `rank1.pdb`.

**Scientific notes:** DiffDock-PP models protein-protein docking as a generative diffusion process over the space of rigid-body transformations (translations and rotations). Starting from a random initial placement, the score model iteratively denoises the ligand protein's pose relative to the fixed receptor. The confidence model is a separate classification network trained to predict pose quality. Higher confidence scores correlate with lower RMSD to the native structure (Spearman ρ ≈ 0.68), but the scores are not calibrated probabilities of success. The model uses E3-equivariant graph neural networks (e3nn) with KNN-based graph construction.

**Outputs:**
- `rank1.pdb` — top-ranked docking pose
- `confidence_scores.json` — confidence scores for all generated poses

### Node 05: Interface Analysis

**Goal:** Characterize the predicted antibody-antigen binding interface contacts.

**Process:** Identifies interface residues using a 5.0 Å distance cutoff between antibody and antigen atoms. Classifies contacts by interaction type: hydrogen bonds (< 3.5 Å involving polar residues), hydrophobic (both residues hydrophobic), electrostatic (opposite charges), and polar interactions. Reports contact statistics, interface residue lists, and the top 20 closest contacts.

**Scientific notes:** The binding interface between antibody CDR loops and antigen epitope typically involves 15–25 contact residues per side. A balanced mix of hydrogen bonds and hydrophobic contacts is characteristic of high-affinity antibody-antigen interactions. The 5.0 Å distance cutoff is standard for defining protein-protein interface contacts in structural biology.

**Outputs:**
- `interface_analysis.txt` — human-readable interface summary
- `contact_residues.json` — machine-readable interface residue lists
- `final_complex.pdb` — combined antibody-antigen complex

### Node 06: Docking Comparison

**Goal:** Validate the predicted docking pose against the original experimental co-crystal structure.

**Process:** Superimposes the predicted complex onto the reference structure by aligning antibody Cα atoms (the receptor is the fixed anchor). Computes full RMSD (all antigen Cα atoms after superposition) and interface RMSD (antigen residues within 5.0 Å of the antibody in the reference). Classifies prediction quality and generates a PyMOL visualization script for manual inspection.

**Scientific notes:** RMSD after receptor superposition measures how accurately the model predicted the antigen's binding pose. Quality thresholds follow structural biology conventions: < 2.0 Å (excellent), 2.0–5.0 Å (good), 5.0–10.0 Å (moderate), > 10.0 Å (poor). For reference, the CAPRI (Critical Assessment of PRediction of Interactions) criteria classify protein-protein docking predictions as acceptable when LRMSD < 10 Å and fnat ≥ 0.1; note that CAPRI's LRMSD and fnat differ from the receptor-superposition RMSD reported here, so these thresholds are approximate conventions rather than strict CAPRI classifications. Interface RMSD is more relevant than full RMSD for assessing binding specificity, as it focuses on the residues that form the actual contact.

**Outputs:**
- `rmsd_analysis.json` — full RMSD, interface RMSD, atom counts, pairing statistics
- `best_match.pdb` — copy of the top-ranked predicted pose
- `final_comparison_report.txt` — human-readable quality assessment
- `align_structures.pml` — PyMOL script for visual inspection

## Parameters

### input_pdb

- **Type:** string
- **Default:** `complex.pdb`
- **Description:** Filename of the input antibody-antigen complex PDB in `input_files/`.
- **Guidance:** Set to the name of your PDB file. The sample input is `5B8C.pdb`.

### antibody_chains

- **Type:** string
- **Default:** `""` (auto-detect)
- **Description:** Comma-separated antibody chain IDs (e.g., `H,L`).
- **Guidance:** Leave empty for standard H/L chain naming. Set explicitly if your PDB uses non-standard chain identifiers (e.g., `A,B` for heavy and light chains).

### num_samples

| Value | Description |
|-------|-------------|
| `10` (default) | Generates 10 candidate poses. Sufficient for testing and most benchmarks. |
| `20`–`40` | More poses increase the chance of sampling a near-native conformation. |

**Trade-off:** More samples improve coverage of the conformational space but increase runtime linearly. For benchmarking, 10 samples is standard; for prospective docking where accuracy is critical, 20–40 is recommended.

**Test vs production:** Default of 10 is suitable for testing. For publication-quality results, use 20–40.

### inference_steps

| Value | Description |
|-------|-------------|
| `20` (default) | Number of denoising diffusion steps per sample. Standard for DiffDock-PP. |
| `40` | More steps can improve pose quality at the cost of longer runtime. |

**Trade-off:** More denoising steps allow finer refinement of the pose but with diminishing returns beyond 20 steps.

## Outputs and interpretation

### Confidence scores

DiffDock-PP's confidence model produces a score for each generated pose. Higher scores indicate higher predicted quality. The scores are used for ranking poses (best first) but are not calibrated probabilities — a confidence of 0.8 does not mean 80% chance of success. Use the top-ranked pose (rank1.pdb) for downstream analysis.

### RMSD metrics

- **Full RMSD:** Backbone Cα RMSD of the antigen after superposing the antibody chains. Measures overall positioning accuracy. Values < 5.0 Å indicate good global agreement with the reference structure.
- **Interface RMSD:** RMSD computed only for antigen residues at the binding interface (< 5.0 Å from antibody in the reference). More informative than full RMSD for assessing binding mode accuracy, as it focuses on the contact region.

Quality classification:
| RMSD (Å) | Quality |
|-----------|---------|
| < 2.0 | Excellent — very close to experimental structure |
| 2.0–5.0 | Good — reasonable agreement |
| 5.0–10.0 | Moderate — significant differences, binding site may be partially correct |
| > 10.0 | Poor — substantial deviation from native pose |

### Interface contacts

The interface analysis reports the number and type of antibody-antigen contacts. A well-predicted binding interface typically shows 15–25 contact residues per side with a mix of hydrogen bonds (30–40%) and hydrophobic contacts (20–30%). Dominance of a single interaction type may indicate a non-specific or artifactual interface.

## Quick start

### Running with Docker

The workflow uses the `diffdock_abag:v1` Docker image (NVIDIA CUDA 11.8, PyTorch 2.0.1, ESM-2, DiffDock-PP). GPU is required for inference.

### Running on Silva

1. Select workflow-016 from the workflow list
2. Upload your antibody-antigen complex PDB to `input_files/`
3. Set `input_pdb` to your filename; adjust `antibody_chains` if needed
4. Click Run

### Test vs production settings

| Setting | Test (default) | Production |
|---------|---------------|------------|
| `input_pdb` | `complex.pdb` (sample 5B8C) | Your target complex |
| `num_samples` | `10` | `20`–`40` |
| `inference_steps` | `20` | `20`–`40` |

A successful test run with the sample 5B8C complex produces ranked docking poses, an interface analysis report, and RMSD comparison against the original crystal structure.

## References

- Ketata, M.A., Laue, C., Mammadov, R., Stark, H., Wu, M., Corso, G., Marquet, C., Barzilay, R. & Jaakkola, T.S. "DiffDock-PP: Rigid Protein-Protein Docking with Diffusion Models." *ICLR 2023 Workshop on Machine Learning for Drug Discovery*. arXiv:2304.03889.
- [DiffDock-PP source code](https://github.com/ketatam/DiffDock-PP)
- Lin, Z. et al. "Evolutionary-scale prediction of atomic-level protein structure with a language model." *Science* 379(6637):1123–1130, 2023. DOI: https://doi.org/10.1126/science.ade2574
- Basu, S. & Wallner, B. "DockQ: A Quality Measure for Protein-Protein Docking Models." *PLOS ONE* 11(8):e0161879, 2016. DOI: https://doi.org/10.1371/journal.pone.0161879
