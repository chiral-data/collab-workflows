---
doc_id: workflow-007
domain: structural-biology
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  Protein binding pocket detection using pocketeer's alpha-sphere algorithm.
  Downloads PDB structures, identifies cavities via Delaunay tessellation and
  alpha-sphere clustering, and generates an interactive 3D HTML report with
  optional rotating GIF.
tags: [pocketeer, pocket-detection, alpha-sphere, protein-structure, fpocket, 3dmol]
---

# Workflow 007: Protein Pocket Analysis

Protein binding pocket detection using [pocketeer](https://github.com/cch1999/pocketeer). This workflow downloads one or more protein structures from RCSB PDB, identifies binding pockets using alpha-sphere geometry, and produces an interactive HTML report with 3Dmol.js visualization and an optional rotating GIF via PyMOL.

## Overview

Pocketeer is a Python reimplementation of the [fpocket](https://github.com/Discngine/fpocket) algorithm for detecting cavities in protein structures. The method works by computing a Delaunay tessellation of the protein's atoms, extracting alpha spheres (circumspheres of the resulting tetrahedra), filtering for buried spheres based on solvent-accessible surface area (SASA), and clustering nearby spheres into pockets. Each pocket is characterized by its volume, score, centroid, and the set of protein residues that line the cavity (Le Guilloux et al., 2009).

The workflow processes multiple PDB structures in a single run, enabling side-by-side comparison of pockets across different proteins or conformations.

## When to use this workflow

Use this workflow when you want to identify potential binding sites on a protein structure before running molecular docking. It is useful as a preliminary step to locate druggable pockets, compare binding sites across homologs, or identify allosteric sites.

For docking small molecules into identified pockets, use workflow-004 (AutoDock Vina) or workflow-002/003 (Smina). For protein-protein docking, use workflow-016 (DiffDock-PP) or workflow-017 (LightDock). For protein function prediction from sequence, use workflow-015 (mDeepFRI).

## Architecture and data flow

```text
[01: Download PDB] ──> [02: Calculate Pockets] ──> [03: Generate Visualization]
        |                        |                            |
   {PDB_ID}.pdb          {PDB_ID}_pockets.json       pocket_visualization.html
                          config.json                 {PDB_ID}_pockets_*.gif
```

Nodes run sequentially: 01 → 02 → 03.

## Input requirements

- **PDB IDs:** A comma-separated list of RCSB PDB identifiers (default: `1FME, 2RH1, 4TOS`). Structures are downloaded automatically from RCSB.
- **Sample proteins:**
  - `1FME` — HIV-1 protease (classic drug target with a clear binding pocket)
  - `2RH1` — Beta-2 adrenergic receptor (GPCR with a well-defined ligand pocket)
  - `4TOS` — Tankyrase 1 (PARP family, Wnt signaling drug target)

## Workflow nodes

### Node 01: Download PDB

**Goal:** Fetch protein structures from the RCSB PDB database.

**Process:** Parses the comma-separated `pdb_ids` parameter, uppercases each ID, and downloads the corresponding PDB file from `https://files.rcsb.org/download/{PDB_ID}.pdb`.

**Outputs:**
- `{PDB_ID}.pdb` — one PDB file per requested structure

### Node 02: Calculate Pockets

**Goal:** Detect binding pockets in each protein structure using the alpha-sphere algorithm.

**Process:** For each PDB file in the inputs directory:
1. Loads the structure with `pocketeer.load_structure()`
2. Runs `pocketeer.find_pockets()` with the configured alpha-sphere parameters
3. Saves detected pockets to `{PDB_ID}_pockets.json` using `pocketeer.write_pockets_json()`
4. Copies the PDB file to outputs for downstream visualization
5. Writes a `config.json` listing all processed PDB IDs

Filenames are validated to contain only alphanumeric characters, hyphens, underscores, and dots to prevent issues in the HTML visualization.

**Scientific notes:** The algorithm performs Delaunay tessellation of protein atoms and extracts alpha spheres — the circumspheres of the resulting tetrahedra. Spheres within the radius range (`r_min` to `r_max`) that are buried (low SASA, determined by `polar_probe_radius`) are retained. Nearby spheres are merged based on `merge_distance`, and clusters with fewer than `min_spheres` spheres are discarded. Each pocket is characterized by its volume (voxel-based estimate in ų), a geometric score, centroid coordinates, and the list of lining residues with chain/residue IDs.

**Outputs:**
- `{PDB_ID}_pockets.json` — per-pocket data: pocket_id, score, volume, centroid, spheres (center, radius, SASA), and lining residues
- `{PDB_ID}.pdb` — copied PDB file
- `config.json` — list of processed PDB IDs

### Node 03: Generate Visualization

**Goal:** Build an interactive HTML report with 3D pocket visualization and optional rotating GIF.

**Process:** Reads `config.json` and the corresponding PDB/pocket JSON files. For each protein, generates a 3Dmol.js interactive viewer using `pocketeer.view_pockets()`. Assembles all viewers into a multi-protein HTML dashboard with:
- Grid layout with per-protein cards showing pocket count and top score
- Interactive 3D viewers with multiple rendering styles (cartoon, stick, sphere, line, surface)
- Focus mode (full-screen single protein) and compare mode (side-by-side selected proteins)
- Keyboard shortcut: Escape to exit focus/compare mode

If `output_format` is `gif` or `both`, additionally renders a rotating 360-degree GIF for each protein using PyMOL, with the top 5 pockets highlighted in distinct colors (red, orange, yellow, green, magenta) as filled sphere surfaces or single centroid spheres.

**Outputs:**
- `pocket_visualization.html` — self-contained interactive HTML dashboard
- `{PDB_ID}_pockets_{representation}.gif` — optional rotating GIF per protein (when `output_format` is `gif` or `both`)

## Parameters

### pdb_ids

- **Type:** string
- **Default:** `"1FME, 2RH1, 4TOS"`
- **Node:** 01
- **Description:** Comma-separated PDB IDs to download and analyze.

### r_min

- **Type:** float
- **Default:** `3.0`
- **Node:** 02
- **Description:** Minimum alpha-sphere radius in angstroms. Spheres smaller than this are discarded. Smaller values detect tighter cavities.

### r_max

- **Type:** float
- **Default:** `6.0`
- **Node:** 02
- **Description:** Maximum alpha-sphere radius in angstroms. Spheres larger than this are discarded (they represent open, solvent-exposed regions rather than pockets).

### polar_probe_radius

- **Type:** float
- **Default:** `1.8`
- **Node:** 02
- **Description:** Probe radius in angstroms used to test atom contact for SASA-based polarity labeling of alpha spheres.

### merge_distance

- **Type:** float
- **Default:** `1.2`
- **Node:** 02
- **Description:** Distance threshold in angstroms for merging nearby alpha-sphere clusters into a single pocket.

### min_spheres

- **Type:** integer
- **Default:** `35`
- **Node:** 02
- **Description:** Minimum number of alpha spheres required to form a pocket cluster. Lower values detect smaller pockets; higher values filter noise.

### ignore_hydrogens

- **Type:** boolean
- **Default:** `true`
- **Node:** 02
- **Description:** Remove hydrogen atoms before pocket detection. Recommended to leave enabled.

### ignore_water

- **Type:** boolean
- **Default:** `true`
- **Node:** 02
- **Description:** Remove water molecules before pocket detection. Recommended to leave enabled.

### ignore_hetero

- **Type:** boolean
- **Default:** `true`
- **Node:** 02
- **Description:** Remove hetero atoms (ligands, cofactors) before pocket detection. Enable for apo (unliganded) pocket detection; disable if you want ligand atoms to influence the tessellation.

### pocket_style

- **Type:** enum
- **Default:** `"filled_surfaces"`
- **Node:** 03
- **Description:** How pockets are rendered in the visualization.

| Value | Description |
|-------|-------------|
| `filled_surfaces` (default) | Each alpha sphere rendered individually, showing pocket shape |
| `single_sphere` | Single sphere at pocket centroid, scale 2.0 |

### render_method

- **Type:** enum
- **Default:** `"draw"`
- **Node:** 03
- **Description:** PyMOL rendering method for GIF output (ignored for HTML-only output).

| Value | Description |
|-------|-------------|
| `draw` (default) | Fast OpenGL rendering |
| `ray` | Ray-traced rendering, higher quality but slower |

### representation

- **Type:** enum
- **Default:** `"surface"`
- **Node:** 03
- **Description:** Protein representation style for PyMOL GIF output. The HTML viewer defaults to cartoon regardless of this setting but offers all styles interactively.

| Value | Description |
|-------|-------------|
| `surface` (default) | Molecular surface with 30% transparency |
| `cartoon` | Secondary structure cartoon with fancy helices |

### output_format

- **Type:** enum
- **Default:** `"html"`
- **Node:** 03
- **Description:** Output format to generate.

| Value | Description |
|-------|-------------|
| `html` (default) | Interactive HTML dashboard only |
| `gif` | Rotating GIF only (requires PyMOL) |
| `both` | Both HTML and GIF |

## Outputs and interpretation

### Pocket score

A geometric score reflecting the size and shape quality of the detected pocket. Higher scores indicate larger, more well-defined cavities. Scores are not directly comparable across different proteins — compare pockets within the same structure to identify the most prominent binding site.

### Pocket volume

The estimated volume of the cavity in cubic angstroms (ų), calculated using a voxel-based approach. Typical drug-binding pockets range from ~200 to ~1500 ų. Very large volumes may indicate channel-like features rather than discrete binding pockets.

### Lining residues

The set of protein residues (chain ID, residue number, residue name) that form the walls of each pocket. These residues define the binding site for downstream docking workflows.

### pocket_visualization.html

Self-contained HTML dashboard with interactive 3Dmol.js viewers for all analyzed proteins. Supports grid view, focus mode, and compare mode for side-by-side analysis.

## Quick start

### Running with Docker

```bash
docker pull ghcr.io/chiral-data/pocketeer:2025_12_08
```

### Running on Silva

1. Select "Protein Pocket Analysis" from the workflow list
2. Enter comma-separated PDB IDs (or use the defaults: 1FME, 2RH1, 4TOS)
3. Adjust pocket detection parameters if needed (defaults work well for most proteins)
4. Click Run

### Test vs production settings

| Setting | Test (default) | Production |
|---------|---------------|------------|
| `pdb_ids` | `1FME, 2RH1, 4TOS` | Your target PDB IDs |
| `min_spheres` | `35` | `20`–`50` depending on pocket size |
| `output_format` | `html` | `both` (for publication figures) |
| `render_method` | `draw` | `ray` (for GIF quality) |

A successful test run with the three default proteins completes in a few minutes and produces an interactive HTML dashboard showing detected pockets for each protein.

## Troubleshooting

### No pockets detected

If pocketeer finds zero pockets for a protein, try lowering `min_spheres` (e.g., to 20) or widening the radius range (`r_min=2.5`, `r_max=7.0`). Some small or highly exposed binding sites may not form large enough alpha-sphere clusters with the default parameters.

### Missing ligand context

By default, `ignore_hetero=true` removes co-crystallized ligands before pocket detection. If you want the pocket geometry to reflect the ligand-bound conformation, set `ignore_hetero=false`.

## References

- Le Guilloux, V., Schmidtke, P. & Tuffery, P. "Fpocket: An open source platform for ligand pocket detection." *BMC Bioinformatics* 10:168, 2009. DOI: https://doi.org/10.1186/1471-2105-10-168
- [pocketeer GitHub](https://github.com/cch1999/pocketeer)
- [pocketeer documentation](https://pocketeer.readthedocs.io/)
