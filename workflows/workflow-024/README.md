# Workflow 024: qPCR Primer/Probe Design

End-to-end design and validation of qPCR primer/probe sets for any target organism.

**Original tool:** https://github.com/ajaypavan1004/qpcr_pipeline

## Overview

Given a target organism name, this workflow:

1. Downloads target genome sequences and close-relative sequences from NCBI Entrez (or accepts user-provided FASTAs)
2. Validates the FASTA inputs
3. Scans for Regions of Interest (ROIs) unique to the target using k-mer subtraction against close-relative sequences
4. Designs a forward primer, reverse primer, and hydrolysis TaqMan probe for each ROI using Primer3, enforcing Tm, GC, amplicon size, self-dimer, and hairpin constraints
5. Validates the top primer sets for specificity using local BLAST+, then generates JSON, TSV, and human-readable TXT reports

Supports whole-genome mode and gene-targeted mode (e.g. `cpn60`, `rpoB`, `rnaseh`).

## Pipeline

```
00_download_sequences → 01_validate_inputs → 02_roi_selection → 03_primer_design → 04_blast_report
```

## Nodes

| Node | Name                    | Description                                                          | Key Outputs                                    |
|------|-------------------------|----------------------------------------------------------------------|------------------------------------------------|
| 00   | Download Sequences      | Fetch target + exclusion FASTAs from NCBI, or copy from input_files/ | `target.fasta`, `exclusion.fasta`              |
| 01   | Validate Inputs         | Check FASTA format, sequence count, length, nucleotide content       | `validation_summary.json`                       |
| 02   | ROI Selection           | K-mer uniqueness scan; top-N non-overlapping ROI windows             | `rois.json`                                    |
| 03   | Primer Design           | Primer3 design + thermodynamic constraint check                      | `primer_sets.json`                             |
| 04   | BLAST Validation/Report | Local BLAST+ specificity check + JSON/TSV/TXT reports                | `results.json`, `results.tsv`, `report.txt`    |

## Input Files

Place these in `input_files/` to run in **offline mode** (skips NCBI download):

| File | Description |
|------|-------------|
| `target.fasta` | One or more target organism sequences (FASTA format) |
| `exclusion.fasta` | Close-relative sequences for specificity screening (FASTA format). An empty file is accepted — all ROI windows score 1.0 |

Sample test files for *Rickettsia rickettsii* are included.

If these files are absent, Node 00 downloads from NCBI using the `organism` parameter.

## Parameters

### global_params.json (workflow root)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `organism` | string | `"Rickettsia rickettsii"` | Target organism scientific name |
| `email` | string | `"user@example.com"` | Email for NCBI Entrez (required) |
| `ncbi_api_key` | string | `""` | NCBI API key (optional; raises rate limit to 10 req/sec) |
| `target_gene` | string | `""` | Gene name for gene-targeted mode (e.g. `cpn60`, `rpoB`) |
| `max_target_seqs` | integer | `3` | Max target sequences to fetch from NCBI |
| `max_relative_seqs` | integer | `5` | Max close-relative sequences for exclusion |
| `roi_window` | integer | `500` | ROI window size in bp |
| `roi_step` | integer | `100` | Sliding step in bp |
| `min_uniqueness` | float | `0.80` | Min k-mer uniqueness score (0–1) |
| `top_rois` | integer | `3` | Top ROIs to pass to Primer3 |
| `primer_min_size` | integer | `19` | Min primer length (nt) |
| `primer_opt_size` | integer | `20` | Optimal primer length (nt) |
| `primer_max_size` | integer | `26` | Max primer length (nt) — increase to 28 for difficult organisms |
| `primer_min_tm` | float | `59.0` | Min primer Tm (°C) |
| `primer_opt_tm` | float | `60.0` | Optimal primer Tm (°C) |
| `primer_max_tm` | float | `62.0` | Max primer Tm (°C) |
| `amplicon_min` | integer | `70` | Min amplicon size (bp) |
| `amplicon_max` | integer | `200` | Max amplicon size (bp) |
| `blast_sets` | integer | `3` | Number of top sets to BLAST |
| `skip_blast` | boolean | `false` | Skip BLAST validation (design-only mode) |

## Output Files

| File | Description |
|------|-------------|
| `results.json` | Machine-readable full results: all primer sets with Tm, GC, BLAST hits |
| `results.tsv` | Spreadsheet-friendly summary of all primer sets |
| `report.txt` | Human-readable ranked report with per-oligo stats and BLAST results |

## Running

### With NCBI download (online mode)

1. Set `organism` and `email` in `global_params.json`
2. Set `SILVA_WORKFLOW_HOME` to the parent of `workflow-024/`
3. Launch Silva and select `workflow-024`
4. Press Enter

### With local FASTAs (offline mode)

1. Place `target.fasta` and `exclusion.fasta` in `input_files/`
2. Set `organism` in `global_params.json` (used for BLAST specificity evaluation)
3. Follow steps 2–4 above

### Gene-targeted mode

For organisms where whole-genome ROIs are impractical (e.g. large or fragmentary assemblies), set `target_gene` to a specific gene (e.g. `cpn60`, `rpoB`). The pipeline will fetch gene sequences for both target and exclusion, and auto-adjust the ROI window to fit shorter sequences.

## Requirements

- Docker
- Silva (https://github.com/chiral-data/silva)
- Internet access (online mode only)

Build the Docker image:

```bash
cd workflow-024
./build.sh
```
