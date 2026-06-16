---
doc_id: workflow-026
domain: structure-prediction
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  AlphaFold2-based 3D protein structure prediction from a FASTA sequence
  using ColabFold with MMseqs2 MSA, confidence analysis, and PyMOL visualization.
tags: [colabfold, alphafold2, structure-prediction, msa, multimer]
---

# Workflow 026: ColabFold Structure Prediction

AlphaFold2-based 3D protein structure prediction using [ColabFold](https://github.com/sokrypton/ColabFold) with MSA from MMseqs2. Supports both monomer and multimer inputs and produces ranked PDB structures, per-residue confidence scores, and visualizations.

## Overview

This workflow takes a protein sequence (or complex of sequences) in FASTA format, validates it, builds a multiple sequence alignment via the public MMseqs2 server, and runs AlphaFold2 inference with GPU acceleration. Up to 5 ranked models are produced per run, each accompanied by a score JSON containing per-residue pLDDT and an N×N PAE matrix. A confidence analysis node aggregates those scores into a structured summary, and a visualization node renders a pLDDT-colored structure PNG and a PAE heatmap with chain boundaries for multimers.

ColabFold (Mirdita et al., 2022) combines the speed of MMseqs2-based MSA search with the accuracy of AlphaFold2 (Jumper et al., 2021), enabling structure prediction in minutes rather than hours. All four nodes run in a single unified Docker image (`colabfold:2026_06_09`) that includes ColabFold, AlphaFold2 weights, PyMOL, and all Python dependencies.

## When to use this workflow

Use this workflow when you have a protein sequence (50–2500 amino acids) and need a predicted 3D structure. It handles both monomers (single FASTA record) and homomeric or heteromeric complexes (multiple FASTA records or `:` separator in the header). The workflow is appropriate for generating structural hypotheses, identifying likely fold, or preparing a starting model for molecular docking or design.

Do not use this workflow if you already have an experimental structure — use the experimental model directly instead. For sequences longer than 2500 residues, colabfold_batch will fail validation; contact the dev team to discuss chunking strategies. If you need energy-minimized (relaxed) structures, this workflow produces unrelaxed PDBs only. For binding affinity estimation after structure prediction, chain into workflow-025.

## Architecture and data flow

```text
input_files/sequences.fasta
         |
         v
[01: FASTA Validation] ──> validated_sequences.fasta
         |
         v
[02: Structure Prediction] ──> *_unrelaxed_rank_*.pdb
    (GPU, colabfold_batch)      *_scores_rank_*.json
                                validated_sequences.fasta (pass-through)
         |                               |
         v                               |
[03: Confidence Analysis] <─────────────┘
    (score JSONs + FASTA)
         |
         v
   confidence_summary.json
         |
         v
[04: Structure Visualization]
    (rank 001 PDB + summary)
         |
         v
   structure_plddt.png
   pae_matrix.png
```

Node 02 depends on 01; node 03 depends on 02; node 04 depends on both 02 and 03. Nodes 01–03 run sequentially; node 04 waits for both 02 and 03 before starting.

## Input requirements

A standard FASTA file (`sequences.fasta`) placed in `input_files/`:

```
>protein_name
MVLSPADKTNVKAAWGKVGA...
```

For multimers, provide one record per chain as separate FASTA entries:

```
>ChainA
MNIFEMLRIDE...
>ChainB
ACDEFGHIKLM...
```

Alternatively, join chains with `:` in a single sequence header (ColabFold convention). Sequence lengths must be between `min_length` (default 10) and `max_length` (default 2500) amino acids. Only standard amino acid characters plus `X` (unknown) are accepted.

A sample test input is provided at `input_files/sequences.fasta`: T4 lysozyme (164 residues, PDB reference 2LZM), a canonical AlphaFold2 benchmark monomer.

## Workflow nodes

### Node 01: FASTA Validation

**Goal:** Validate the input FASTA file and detect monomer vs. multimer mode before committing GPU time.

**Process:** Reads `sequences.fasta` from `./inputs/`, checks each sequence against the allowed amino acid character set (`ACDEFGHIKLMNPQRSTVWYX`), enforces `PARAM_MIN_LENGTH` and `PARAM_MAX_LENGTH` bounds, and detects mode: a single sequence is monomer; multiple records or a `:` in the header are treated as multimer. Writes the validated sequences to `validated_sequences.fasta` unchanged (normalization only, no sequence modification).

**Scientific notes:** Early validation prevents colabfold_batch from failing mid-run due to malformed input. The `X` character is permitted because it represents unknown amino acids in some database sequences and ColabFold handles it internally. Multimer detection at this stage allows downstream nodes to select the correct MSA pairing mode.

**Outputs:**
- `validated_sequences.fasta` — cleaned, validated FASTA passed to node 02

### Node 02: Structure Prediction

**Goal:** Build a multiple sequence alignment and run AlphaFold2 inference to produce ranked predicted structures.

**Process:** Copies the FASTA to a `cf_input/` directory (colabfold_batch requires a directory input), then invokes `colabfold_batch cf_input cf_output` with flags derived from parameters: `--num-models`, `--num-recycle`, `--msa-mode`, `--pair-mode`, `--host-url`, and optionally `--templates`. For `pair_mode = auto` the node resolves to `unpaired` for monomers and `unpaired_paired` for multimers. After the run, all output files are copied flat to `./outputs/`. The script exits with a clear error if no GPU is detected via `nvidia-smi`.

**Scientific notes:** MMseqs2-based MSA search (via the public ColabFold server at `api.colabfold.com`) is orders of magnitude faster than jackhmmer/HHblits used by the original AlphaFold2. Paired MSA for multimers captures inter-chain coevolution signals critical for accurate complex geometry. Recycling (controlled by `num_recycle`) refines the structure by feeding the output back as input; 3 iterations are sufficient for monomers but 6–12 are recommended for multimers.

**Outputs:**
- `*_unrelaxed_rank_001_*.pdb` through `rank_005` — ranked predicted structures (unrelaxed)
- `*_scores_rank_001_*.json` through `rank_005` — per-residue pLDDT arrays, PAE matrices, PTM/iPTM scores
- `validated_sequences.fasta` — pass-through for node 03 to read chain lengths

### Node 03: Confidence Analysis

**Goal:** Parse ColabFold score JSONs and emit a structured confidence summary with aggregated metrics.

**Process:** Globs `./inputs/*_scores_rank_*.json`, sorts by rank, and takes the top N models (`PARAM_TOP_N_MODELS`). For each model it extracts `plddt_mean`, `pae_mean`, `max_pae`, `ptm`, and `iptm` (multimer only). Multimer mode is detected when any score JSON contains an `iptm` key. For multimers, chain lengths are read from `validated_sequences.fasta` via BioPython and used to compute `interface_pae_mean` as the mean of the off-diagonal PAE matrix blocks (inter-chain residue pairs only). Results are written to `confidence_summary.json`.

**Scientific notes:** pLDDT (predicted local-distance difference test) is a per-residue confidence score from 0–100; values above 90 indicate high confidence, 70–90 moderate, below 70 low confidence typically associated with disordered regions. PAE (predicted aligned error) measures the expected positional error in Å when the structure is aligned on a reference residue; low off-diagonal PAE between chains indicates a confident prediction of relative chain orientation.

**Outputs:**
- `confidence_summary.json` — aggregated confidence metrics (see schema in [Outputs and interpretation](#outputs-and-interpretation))

### Node 04: Structure Visualization

**Goal:** Render a pLDDT-colored structure image and a PAE heatmap for the top-ranked model.

**Process:** Globs `./inputs/*_unrelaxed_rank_001_*.pdb` for the top-ranked structure and reads `./inputs/confidence_summary.json` for per-residue pLDDT values, the PAE matrix, and chain lengths. For the structure image: attempts to import PyMOL and if available renders a headless ray-traced PNG with residues colored by pLDDT B-factors (`red_yellow_cyan_blue` spectrum, 0–100); falls back to a matplotlib color-scale figure if PyMOL is unavailable at runtime. For the PAE heatmap: uses `matplotlib.pyplot.imshow` with the `RdYlGn_r` colormap and draws white `axvline`/`axhline` at chain boundaries for multimers using cumulative chain lengths from `chain_lengths` in the summary.

**Scientific notes:** pLDDT coloring maps confidence directly onto 3D structure — predominantly blue/teal regions are reliable, red/orange regions are likely disordered or poorly predicted. PAE heatmaps are the primary tool for assessing multimer interface confidence: a well-predicted interface shows low PAE (green) in both off-diagonal quadrants symmetrically.

**Outputs:**
- `structure_plddt.png` — structure colored by pLDDT (PyMOL ray-traced or matplotlib fallback)
- `pae_matrix.png` — PAE heatmap with chain boundary lines for multimers

## Parameters

### `num_models`

| Value / Range | Description |
|---------------|-------------|
| `5` (default) | Run all 5 AlphaFold2 model weights. Best coverage of conformational space; use for any result you intend to report. |
| `1`–`2` | Faster runs for pipeline validation or exploratory work where you only need a rough structure. |

**Trade-off:** More models increases confidence in the top-ranked prediction and surfaces ensemble diversity, at proportionally higher GPU time.

**Test vs production:** Default of 5 is production-quality. Use `2` for quick validation runs.

### `num_recycle`

| Value / Range | Description |
|---------------|-------------|
| `3` (default) | 3 recycling iterations. Sufficient for most monomers. |
| `6`–`12` | Recommended for multimers; helps AlphaFold2 converge on correct inter-chain geometry. |

**Trade-off:** Higher recycle count increases run time roughly linearly; returns diminish above 12.

**Test vs production:** For monomer test runs, `3` is fine. Set to `6` or higher for any multimer intended for reporting.

### `msa_mode`

| Value / Range | Description |
|---------------|-------------|
| `mmseqs2_uniref_env` (default) | Searches UniRef + environmental databases. Best sensitivity; recommended for all production runs. |
| `mmseqs2_uniref` | UniRef only. Slightly faster; may miss hits from metagenomic sequences. |
| `single_sequence` | No MSA. Very fast but significantly lower accuracy. Use only for testing or highly novel sequences with no homologs. |

**Trade-off:** Broader MSA search improves prediction accuracy but increases server query time.

### `pair_mode`

| Value / Range | Description |
|---------------|-------------|
| `auto` (default) | Resolves to `unpaired` for monomers and `unpaired_paired` for multimers. |
| `unpaired` | MSA columns are not paired across chains. Use for monomers or when server pairing fails. |
| `unpaired_paired` | Paired MSA for multimers captures inter-chain coevolution. Required for accurate complex geometry. |

**Trade-off:** `unpaired_paired` is strictly better for multimers but requires pairable sequences in the database; `auto` handles this transparently.

### `use_templates`

- **Type:** boolean
- **Default:** `false`
- **Description:** Pass `--templates` to colabfold_batch to incorporate PDB templates into the MSA.
- **Guidance:** Enable if the target has close PDB homologs (>30% identity) and you want template-guided modeling. Has minimal effect on well-conserved proteins where MSA alone is informative.

### `host_url`

- **Type:** string
- **Default:** `https://api.colabfold.com`
- **Description:** MMseqs2 server URL passed as `--host-url` to colabfold_batch.
- **Guidance:** Override with a local MMseqs2 server URL for air-gapped environments or to avoid rate limits on the public server.

### `server_timeout_seconds`

- **Type:** integer
- **Default:** `300`
- **Description:** Timeout in seconds for MSA server requests. The script exits with a clear error on timeout rather than hanging.
- **Guidance:** Increase if the public ColabFold server is under heavy load and requests are failing.

### `top_n_models`

- **Type:** integer
- **Default:** `5`
- **Description:** Number of ranked models to include in `confidence_summary.json`.
- **Guidance:** Reduce to `1` if you only care about the top-ranked structure and want a smaller output file.

### `color_by`

- **Type:** string (enum: `pLDDT`, `chain`, `rainbow`)
- **Default:** `pLDDT`
- **Description:** Color scheme for the structure visualization PNG.
- **Guidance:** `pLDDT` is the most scientifically informative for assessing prediction quality. `chain` is useful for multimers to visually distinguish subunits. `rainbow` is for publication figures.

### `min_length` / `max_length`

- **Type:** integer
- **Defaults:** `10` / `2500`
- **Description:** Sequence length bounds enforced during FASTA validation. Sequences outside this range cause node 01 to exit with an error.
- **Guidance:** Increase `max_length` with caution — colabfold_batch memory usage scales quadratically with sequence length.

## Outputs and interpretation

### `confidence_summary.json`

Aggregated confidence metrics for the top N ranked models. Schema:

```json
{
  "colabfold_version": "1.6.1",
  "job_id": "string",
  "mode": "monomer | multimer",
  "chain_lengths": [164],
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

`interface_pae_mean` and `iptm` are `null` for monomers. `chain_lengths` lists the length of each chain in residues (e.g., `[164]` for a monomer, `[200, 150]` for a two-chain complex).

### `plddt_mean`

Mean per-residue pLDDT across the full sequence. Values above 90 indicate a high-confidence prediction suitable for structural analysis; 70–90 is moderate confidence; below 70 suggests disordered or poorly predicted regions. For T4 lysozyme (the included test case), expect `plddt_mean > 90`.

### `pae_mean`

Mean predicted aligned error across all residue pairs (Å). Lower is better. Values below 5 Å indicate a confident global fold. For T4 lysozyme, expect `pae_mean < 5 Å`. High PAE between two chains in a multimer suggests uncertain relative orientation.

### `ptm` / `iptm`

PTM (predicted TM-score) estimates overall fold confidence on a 0–1 scale; values above 0.5 indicate a predicted fold likely to be correct. iPTM (interface PTM) is a multimer-only score measuring confidence in inter-chain contacts; values above 0.6 suggest a reliable complex prediction. Both are `null` for monomers.

### `interface_pae_mean`

For multimers: mean PAE over all cross-chain residue pairs — matrix elements where residue i belongs to one chain and residue j belongs to a different chain. This is distinct from the matrix diagonal; it is the set of PAE blocks that lie outside all intra-chain diagonal blocks. A low value indicates confident prediction of relative chain orientation and binding interface geometry. `null` for monomers.

### `structure_plddt.png`

PyMOL ray-traced image (or matplotlib fallback) of the top-ranked structure colored by pLDDT using a `red_yellow_cyan_blue` spectrum (red = low confidence, blue = high confidence). A predominantly blue/teal structure indicates a reliable prediction.

### `pae_matrix.png`

PAE heatmap using the `RdYlGn_r` colormap (red = high error, green = low error). The diagonal represents intra-residue alignment (always ~0). For multimers, white lines mark chain boundaries; low PAE in off-diagonal quadrants indicates a confident interface prediction.

## Quick start

### Requirements

- [Docker](https://docs.docker.com/get-docker/) with the `colabfold:2026_06_09` image available
- NVIDIA GPU with [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html) installed
- [Silva](https://chiral.bio/silva) for orchestrated pipeline execution

### Running with Docker

The workflow uses the unified image `colabfold:2026_06_09`. To test an individual node, mount its directory as the workspace and run `run.sh` directly. Example for node 02:

```bash
mkdir -p /tmp/test-02/inputs /tmp/test-02/outputs
cp workflows/workflow-026/input_files/sequences.fasta /tmp/test-02/inputs/validated_sequences.fasta

docker run --rm --gpus all \
  -v /tmp/test-02/inputs:/workspace/inputs \
  -v /tmp/test-02/outputs:/workspace/outputs \
  -v $(pwd)/workflows/workflow-026/02_predict:/workspace \
  -w /workspace \
  -e PARAM_NUM_MODELS=2 \
  -e PARAM_NUM_RECYCLE=3 \
  colabfold:2026_06_09 bash run.sh
```

### Running on Silva

1. Select **workflow-026: ColabFold Structure Prediction** from the workflow list
2. Upload `sequences.fasta` to `input_files/`
3. Adjust parameters if needed (see Parameters section above)
4. Press **Enter** to run

### Dry-run (no GPU required)

To verify the `colabfold_batch` command that would be executed without running inference:

```bash
cd workflows/workflow-026/02_predict
cp ../input_files/sequences.fasta ./inputs/validated_sequences.fasta
python script.py --dry-run
```

This prints the fully resolved `colabfold_batch` command and exits 0 without touching a GPU.

### Test vs production settings

| Setting | Test | Production |
|---------|------|------------|
| `num_models` | `2` | `5` (default) |
| `num_recycle` | `3` (default) | `3` (monomer) / `6–12` (multimer) |
| `msa_mode` | `mmseqs2_uniref_env` (default) | `mmseqs2_uniref_env` (default) |
| `use_templates` | `false` (default) | `true` (if PDB homologs exist) |

The included test input (`input_files/sequences.fasta`) is T4 lysozyme — 164 residues, canonical AlphaFold2 benchmark. A correct test run produces `plddt_mean > 90`, `pae_mean < 5 Å`, and a TM-score > 0.95 vs PDB [2LZM](https://www.rcsb.org/structure/2LZM).

## Troubleshooting

**Node 02 exits with "No GPU detected"**
The container requires `--gpus all` (Docker) or `use_gpu = true` in the job.toml (Silva). Confirm the NVIDIA Container Toolkit is installed and `nvidia-smi` runs inside the container.

**MSA server timeout**
The public ColabFold server at `api.colabfold.com` can be slow under heavy load. Increase `server_timeout_seconds` or set `host_url` to a local MMseqs2 server.

**Node 01 exits with "Invalid amino acid character"**
The input FASTA contains a non-standard character. Check for DNA/RNA sequences, lowercase letters, or gap characters (`-`, `.`). Only `ACDEFGHIKLMNPQRSTVWYX` are accepted.

## References

- Mirdita M, Schütze K, Moriwaki Y, Heo L, Ovchinnikov S, Steinegger M. "ColabFold: making protein folding accessible to all." *Nature Methods* 19:679–682, 2022. DOI: https://doi.org/10.1038/s41592-022-01488-1
- Jumper J et al. "Highly accurate protein structure prediction with AlphaFold." *Nature* 596:583–589, 2021. DOI: https://doi.org/10.1038/s41586-021-03819-2
- Steinegger M, Söding J. "MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets." *Nature Biotechnology* 35:1026–1028, 2017. DOI: https://doi.org/10.1038/nbt.3988
- [ColabFold GitHub](https://github.com/sokrypton/ColabFold)
