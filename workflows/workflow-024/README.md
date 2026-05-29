# Workflow 024: qPCR Primer/Probe Design

End-to-end design and validation of qPCR primer/probe sets for any target organism.

**Original Tool:** https://github.com/ajaypavan1004/qpcr_pipeline

## Overview

Given a target organism and sequence source, this workflow:

1. Resolves target and exclusion FASTA sequences — from a user-supplied file (`target_fasta`), from `input_files/target.fasta` + `exclusion.fasta` (full offline mode), or by downloading from NCBI Entrez using `organism`
2. Validates the FASTA inputs (format, sequence count, length, nucleotide alphabet)
3. Scans for Regions of Interest (ROIs) unique to the target using a sliding-window k-mer subtraction against close-relative exclusion sequences
4. Designs a forward primer, reverse primer, and hydrolysis TaqMan probe for each ROI using Primer3, enforcing Tm, GC%, amplicon size, probe Tm-delta, self-dimer, and hairpin constraints
5. Validates the top primer sets for specificity using a local organism-scoped BLAST+ database, then generates JSON, TSV, and human-readable TXT reports

Supports whole-genome mode, gene-targeted mode (e.g. `cpn60`, `rpoB`), and reference-FASTA mode via `target_fasta`.

## Pipeline

```
00_download_sequences → 01_validate_inputs → 02_roi_selection → 03_primer_design → 04_blast_report
```

## Nodes

| Node | Name                    | Description                                                                               | Key Outputs                                 |
|------|-------------------------|-------------------------------------------------------------------------------------------|---------------------------------------------|
| 00   | Download Sequences      | Resolve target + exclusion FASTAs from inputs/ or NCBI                                   | `target.fasta`, `exclusion.fasta`           |
| 01   | Validate Inputs         | Check FASTA format, sequence count, length, nucleotide content                            | `validation_summary.json`                   |
| 02   | ROI Selection           | Sliding-window k-mer uniqueness scan; select top-N non-overlapping ROI windows            | `rois.json`                                 |
| 03   | Primer Design           | Primer3 design + thermodynamic constraint filtering (exit code 2 = no primers, not error) | `primer_sets.json`                          |
| 04   | BLAST Validation/Report | Local BLAST+ specificity check + JSON/TSV/TXT/HTML reports (exit code 2 = no sets pass)  | `results.json`, `results.tsv`, `report.txt`, `report.html` |

## Input Files

Place files in `input_files/` for offline or hybrid operation.

| File | Description |
|------|-------------|
| `<name>.fasta` | Any target sequence FASTA — reference via `target_fasta` in global_params.json |
| `exclusion.fasta` | Close-relative sequences for specificity screening. An empty file is accepted — all ROI windows score 1.0 |

### Operating modes (Node 00 resolution order)

1. **`target_fasta` mode** (hybrid): `target_fasta` is set in `global_params.json` → Node 00 copies `inputs/<target_fasta>` as the target. For exclusion, it uses `inputs/exclusion.fasta` if present, otherwise fetches close relatives from NCBI using `organism`.
2. **Full offline mode**: `target_fasta` is empty, `organism` is empty, and both `target.fasta` and `exclusion.fasta` are present in `input_files/` → Node 00 copies both files through.
3. **NCBI download mode**: fetches target sequences (whole-genome or gene-targeted via `target_gene`) and close-relative exclusion sequences from NCBI Entrez.

## Parameters

### global_params.json (workflow root)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `organism` | string | `""` | Target organism scientific name (for NCBI queries and BLAST database scoping) |
| `email` | string | `"user@example.com"` | Email for NCBI Entrez API (required for NCBI modes) |
| `ncbi_api_key` | string | `""` | NCBI API key (optional; raises rate limit to 10 req/sec) |
| `target_fasta` | string | `""` | Filename of a FASTA in `input_files/` to use as the target (skips NCBI target fetch) |
| `target_gene` | string | `""` | Gene name for gene-targeted NCBI mode (e.g. `cpn60`, `rpoB`). Ignored when `target_fasta` is set. |
| `max_target_seqs` | integer | `3` | Max target sequences to fetch from NCBI |
| `max_relative_seqs` | integer | `5` | Max close-relative sequences for exclusion |
| `roi_window` | integer | `500` | ROI window size in bp |
| `roi_step` | integer | `100` | Sliding step in bp |
| `min_uniqueness` | float | `0.80` | Min k-mer uniqueness score (0–1) |
| `top_rois` | integer | `3` | Top ROIs to pass to Primer3 |
| `primer_min_size` | integer | `19` | Min primer length (nt) |
| `primer_opt_size` | integer | `20` | Optimal primer length (nt) |
| `primer_max_size` | integer | `26` | Max primer length (nt) |
| `primer_min_tm` | float | `59.0` | Min primer Tm (°C) |
| `primer_opt_tm` | float | `60.0` | Optimal primer Tm (°C) |
| `primer_max_tm` | float | `62.0` | Max primer Tm (°C) |
| `primer_min_gc` | float | `40.0` | Min primer GC% |
| `primer_max_gc` | float | `60.0` | Max primer GC% |
| `amplicon_min` | integer | `70` | Min amplicon size (bp) |
| `amplicon_max` | integer | `200` | Max amplicon size (bp) |
| `blast_sets` | integer | `3` | Number of top primer sets to BLAST |
| `skip_blast` | boolean | `false` | Skip BLAST validation (design-only mode) |

## Output Files

| File | Description |
|------|-------------|
| `results.json` | Machine-readable full results: all primer sets with Tm, GC%, amplicon size, BLAST hits |
| `results.tsv` | Spreadsheet-friendly summary of all primer sets |
| `report.txt` | Human-readable ranked report with per-oligo stats and BLAST results |
| `report.html` | Self-contained visual HTML report: ranked primer sets with sequences, Tm, GC%, amplicon size, BLAST results, and pass/fail status |

## Exit Codes

| Exit Code | Meaning |
|-----------|---------|
| `0` | Pipeline completed and at least one primer set passed all constraints. |
| `2` | Pipeline completed successfully but **no primer sets were designed or passed all constraints**. This is a scientific result, not a technical failure — the pipeline ran to completion and all outputs were written. Nodes 03 and 04 both return exit code 2 for their respective no-result conditions; the harness treats this as success. |

### What to do when exit code 2 is returned

Inspect the output files in `04_blast_report/outputs/` to understand which constraints failed:

- **`report.txt`** — Human-readable ranked report. Each primer set shows which thermodynamic, GC, and BLAST constraints it passed or failed.
- **`results.tsv`** — Spreadsheet-friendly summary. Sort and filter by Tm, GC%, amplicon size.
- **`results.json`** — Machine-readable full results with all constraint values.

Then relax the offending parameters in `global_params.json` and re-run:

| Symptom (from report.txt) | Parameter to relax |
|---------------------------|--------------------|
| All sets fail Tm check | Widen `primer_min_tm` / `primer_max_tm` range |
| All sets fail amplicon size | Widen `amplicon_min` / `amplicon_max` |
| No ROIs found upstream | Lower `min_uniqueness` or increase `roi_window` |
| All sets fail BLAST specificity | Review exclusion sequences, or set `skip_blast: true` to inspect designs first |

## Requirements

- Docker
- Silva (https://github.com/chiral-data/silva)
- Internet access (NCBI modes only)

Build the Docker image:

```bash
cd workflow-024
./build.sh
```

## Running

### With a local reference FASTA (target_fasta mode)

1. Place the reference FASTA in `input_files/` (e.g. `vly_reference.fasta`)
2. Set `target_fasta` to the filename in `global_params.json`
3. Optionally place `exclusion.fasta` in `input_files/`; if absent, close relatives are fetched from NCBI using `organism`
4. Launch Silva and select `workflow-024`

### With NCBI download (online mode)

1. Set `organism` and `email` in `global_params.json`
2. Optionally set `target_gene` for gene-targeted mode (e.g. `cpn60`, `rpoB`)
3. Launch Silva and select `workflow-024`

### Full offline mode

1. Place `target.fasta` and `exclusion.fasta` in `input_files/` Note: the file names must be exactly target.fasta and exclusion.fasta - the pipeline looks for these specific names.
2. Leave `organism` and `target_fasta` empty in `global_params.json`
3. **Set `skip_blast: true` in `global_params.json`** — BLAST requires an organism name to scope its database; without one the pipeline auto-sets `skip_blast=true` and logs a warning anyway, but setting it explicitly avoids the warning and makes the intent clear.
4. Launch Silva and select `workflow-024`


