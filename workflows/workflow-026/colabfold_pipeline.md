# ColabFold Structure Prediction Pipeline

---

## 01_validate_fasta

**Purpose:** Validate input FASTA file format, sequence length bounds, allowed amino acid character set, and single vs. multi-chain (multimer) configuration.

**Inputs:** `sequences.fasta` — user-provided FASTA with one or more protein sequences; multimer chains separated by `:` or as multi-FASTA.

**Outputs:** `validated_sequences.fasta` — cleaned, normalized FASTA passed to MSA building.

**Dependencies:** none

**Docker image:** `python:3.11-slim`

**Python/system packages:** `biopython`

**Key parameters:**

| Parameter | Type | Default | Notes |
|-----------|------|---------|-------|
| `input_file` | string | `"sequences.fasta"` | |
| `min_length` | int | `10` | |
| `max_length` | int | `2500` | |
| `mode` | string | `"auto"` | Auto-detects monomer vs. multimer from sequence count |

---

## 02_predict

**Purpose:** Build a multiple sequence alignment via the MMseqs2 server and run AlphaFold2 inference via `colabfold_batch` in a single step. For multimers, uses paired MSA mode to capture inter-chain coevolution signals. The intermediate alignment files are not surfaced as outputs — MSA and inference are kept together as a single atomic step.

**Inputs:** `validated_sequences.fasta` (from `01_validate_fasta`)

**Outputs:**
- `predictions/<job_id>_unrelaxed_rank_001_*.pdb` through `rank_005` — ranked predicted structures
- `predictions/<job_id>_scores_rank_001.json` through `rank_005` — per-residue pLDDT arrays, PAE matrices, predicted TM-scores

**Dependencies:** `01_validate_fasta`

**Docker image:** `ghcr.io/sokrypton/colabfold:1.6.1-cuda12`

**Python/system packages:** `colabfold` (included in image)

**Key parameters:**

| Parameter | Type | Default | Notes |
|-----------|------|---------|-------|
| `num_models` | int | `5` | |
| `--num-recycle` | int | `3` for monomers | Recommended 6–12 for multimers — set explicitly per run type |
| `msa_mode` | string | `"mmseqs2_uniref_env"` | Options: `mmseqs2_uniref_env`, `mmseqs2_uniref`, `single_sequence` |
| `pair_mode` | string | `"unpaired_paired"` | Used for multimers; set to `"unpaired"` for monomers |
| `--templates` | bool | `false` | |
| `use_gpu` | bool | `true` | Pipeline fails explicitly if no GPU is available rather than silently falling back to CPU |
| `msa_server_url` | string | `"https://api.colabfold.com"` | Override for local MMseqs2 server |
| `server_timeout_seconds` | int | `300` | Raises explicit error on timeout rather than hanging |

---

## 03_analyze_confidence

**Purpose:** Parse per-residue pLDDT and PAE from ColabFold score JSONs; compute aggregated confidence metrics per model and emit a structured summary with a defined schema.

**Inputs:** `predictions/<job_id>_scores_rank_001.json` through `rank_005` (from `02_predict`)

**Outputs:** `analysis/confidence_summary.json` (schema below)

**Dependencies:** `02_predict`

**Docker image:** `python:3.11-slim`

**Python/system packages:** `biopython`, `numpy`, `pandas`

**Key parameters:**

| Parameter | Type | Default | Notes |
|-----------|------|---------|-------|
| `top_n_models` | int | `5` | Number of ranked models to include in summary |

**Output schema for `confidence_summary.json`:**

```json
{
  "colabfold_version": "1.6.1",
  "job_id": "string",
  "mode": "monomer" | "multimer",
  "models": [
    {
      "rank": 1,
      "model_name": "string",
      "plddt_mean": 0.0,
      "plddt_per_residue": [0.0],
      "pae_mean": 0.0,
      "max_pae": 0.0,
      "pae": [[0.0]],
      "interface_pae_mean": 0.0 | null,
      "ptm": 0.0,
      "iptm": 0.0 | null
    }
  ]
}
```

> **Notes:**
> - `interface_pae_mean` and `iptm` are `null` for monomers.
> - `pae` is always N×N where N = total residue count across all chains. Field names (`pae`, `max_pae`) match ColabFold's native score JSON output directly.
> - `iptm` is present in the ColabFold score JSON for multimers and is passed through as-is.

---

## 04_pymol_render

**Purpose:** Render the top-ranked predicted structure colored by pLDDT and generate a PAE heatmap. Falls back to py3Dmol if `pymol-open-source` is unavailable.

**Inputs:**
- `predictions/<job_id>_unrelaxed_rank_001_*.pdb` (from `02_predict`)
- `analysis/confidence_summary.json` (from `03_analyze_confidence`)

**Outputs:**
- `visualizations/structure_plddt.png` — structure colored by pLDDT with color scale legend
- `visualizations/pae_matrix.png` — PAE heatmap; annotated with chain boundaries for multimers

**Dependencies:** `03_analyze_confidence`

**Docker image:** `condaforge/mambaforge:latest`

**Python/system packages:**
- Primary: `pymol-open-source` via conda-forge, pinned to v3.1.0
- Fallback (used automatically if PyMOL install fails): `py3Dmol`, `matplotlib`, `numpy`

> **Note:** `pymol-open-source` is the community open-source build available on conda-forge. This is distinct from the commercial Schrödinger `pymol-bundle`.

**Key parameters:**

| Parameter | Type | Default | Notes |
|-----------|------|---------|-------|
| `color_by` | string | `"pLDDT"` | Options: `pLDDT`, `chain`, `rainbow` |
| `renderer` | string | `"auto"` | Auto-selects PyMOL if available, py3Dmol otherwise |
| `view_angles` | string | `"default"` | |

---

## Sample Input & Test Case

**Protein:** T4 lysozyme (T4L)

**Rationale:** T4 lysozyme is the canonical AlphaFold2 benchmark monomer. It is small (164 residues), folds rapidly, has hundreds of crystal structures in the PDB (reference: 2LZM), and ColabFold is known to predict it with mean pLDDT > 90. This makes it ideal for confirming the pipeline runs correctly end-to-end before committing GPU time to a real target.

**`sequences.fasta` contents:**

```
>T4_lysozyme
MNIFEMLRIDEGLRLKIYKDTEGYYTIGIGHLLTKSPSLNAAKSELDKAIGRNTNGVITKDEAEKLFNQDVDAAVRGILRNAKLKPVYDSLDAVRRAALINMVFQMGETGVAGFTNSLRMLQQKRWDEAAVNLAKSRWYNQTPNRAKRVITTFRTGTWDAYKNL
```

### Expected Outputs and Confidence Thresholds

**1. `predictions/T4_lysozyme_unrelaxed_rank_001_*.pdb`**
- File exists and is valid PDB format
- Chain A only, 164 residues

**2. `predictions/T4_lysozyme_scores_rank_001.json`**
- `plddt` array length: 164
- All values between 0 and 100

**3. `analysis/confidence_summary.json`**
- `plddt_mean` for rank 1 model: > 90.0
- `pae_mean` for rank 1 model: < 5.0 Å
- `interface_pae_mean`: `null` (monomer run)
- `iptm`: `null` (monomer run)

**4. `visualizations/structure_plddt.png`**
- File exists and is non-zero bytes
- Structure is predominantly blue/teal (high pLDDT) with little to no red

### Validation

Compare the top-ranked predicted structure against PDB 2LZM using TM-score. A correct prediction should yield TM-score > 0.95:

```bash
tmalign predictions/T4_lysozyme_unrelaxed_rank_001_*.pdb reference/2LZM_A.pdb
```

### Suggested Run Parameters for This Test Case

| Parameter | Value | Notes |
|-----------|-------|-------|
| `num_models` | `2` | Reduced for speed during pipeline validation; use 5 for real runs |
| `--num-recycle` | `3` | |
| `msa_mode` | `mmseqs2_uniref_env` | |
| `--templates` | `false` | |
| `mode` | `monomer` | Auto-detected from single sequence in FASTA |
