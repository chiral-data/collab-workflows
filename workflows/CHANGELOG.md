# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2026-11-12]

### Added

- workflow-038: BioEmu2 Equilibrium Ensemble Emulation — **scaffold, untested against a GPU**
  (6 nodes: 01-prepare-system, 02-fetch-msa, 03-sample, 04-propagate, 05-analyze, 06-report).
  - Conformational ensemble emulation on BioEmu2 (Microsoft Research AI for Science, MIT,
    public preview 2026-10-20). First workflow in the repo that samples how an existing
    protein *moves* rather than predicting one structure or designing a new one.
  - Covers all three BioEmu2 modes by parameter: `predict-structure`, `sample-equilibrium`,
    and time-lagged `propagate`.
  - `02-fetch-msa` is its own CPU node: the alignment is the only step needing the public
    internet, so fetching it to disk turns a server outage or rate limit into a seconds-long
    CPU failure instead of a failed GPU job, and lets `03-sample` run with `HF_HUB_OFFLINE=1`.
    User-supplied `<chain>.a3m` files skip the server entirely.
  - Nodes 03/04 drive the BioEmu2 Python API rather than its CLI, because the CLI exposes no
    flag for a pre-computed MSA; `--num-steps` is deliberately not exposed, as
    `EmuSession.sample()` has no equivalent.
  - `04-propagate` always writes its declared outputs, passing the ensemble through as an
    identity propagation when disabled, because Silva has no conditionally-skipped node and a
    declared-but-unwritten output is a silent failure.
  - Nodes 04-05 bind the XTC to the System's own topology rather than parsing `topology.cif`,
    sidestepping mdtraj's mmCIF support entirely.
  - Shared image `bioemu2:2026_11_11` (`python:3.12-slim`, CUDA 12.6 via the `[cuda12]` extra
    because the `[cuda]` extra needs a CUDA 13 driver the estate lacks; `[cuda12]` also covers
    compute capability 7.0, so V100 works). HF checkpoint pinned to
    `dd4b0c1b69984b68edcb6200a544213be4d75242` (1.91 GiB). Registered in
    `apps/applications.json`.
  - Sampling frequency is **not** a calibrated Boltzmann probability and sampled frames are
    independent draws, not a trajectory; the README and `report.html` both lead with this.
  - Design doc, decisions and open build-time questions in
    `workflows/workflow-038/bioemu2_pipeline.md`.

## [2026-07-24]

### Added

- workflow-034/035/036: real-pipeline counterparts of workflow-030/031/032, split into their own
  workflow numbers so mock and real runs are separate catalog entries instead of one workflow
  number serving both behaviors via `run.sh`.
  - Each node's `run.sh` delegates to `run_real.sh` by default; carries the real computation
    scripts (`build_cell.py`, `analyze.py`, etc.), `.mdp` files, and `apps/*/Dockerfile` from
    its workflow-030/031/032 counterpart. No `run_mock.sh` or `output_files/`/`sample_outputs/`
    — those stay solely on the mock side to avoid duplicating pre-computed data.
  - `.chiral/workflow.toml` titles drop the `(Mock Run)` prefix.
  - Closes #208.

### Changed

- workflow-030/031/032: reverted `run.sh` back to delegating to `run_mock.sh` by default
  (undoing #206/#207's switch to `run_real.sh`), now that the real pipeline lives in
  workflow-034/035/036. `.chiral/workflow.toml` titles keep the `(Mock Run)` prefix, which is
  accurate again.
  - Removed `run_real.sh`, the real computation scripts, `.mdp` files, and `apps/*/Dockerfile`
    from each node — that code now lives solely in workflow-034/035/036 to avoid duplication.
    README.md points to the real-pipeline counterpart for anyone who needs it.
  - Part of #208.

## [2026-07-12]

### Added

- workflow-032: Heat-Resistant Plastics — Tg / Thermal Stability MD (4 nodes: build-cell, cooling-series, measure-tg, report)
  - AMECC Theme 3: estimates glass-transition temperature, thermal expansion, and high-temperature dimensional stability for battery-case/separator resins (PPS, PA66, PBT, PEEK, PP) via a melt-quench cooling series.
  - Reuses the GAFF2/antechamber/acpype pipeline and `polymer_md`/`gromacs` images from workflow-030/031 — no new Dockerfile needed.
  - Node 01 writes both PDB and SDF from RDKit and feeds the SDF (exact bond orders) to antechamber, avoiding the aromatic-backbone tleap valence failure workflow-031 hit for PET.
  - Node 03's bilinear Tg kink fit flags itself unreliable (`tg_reliable: false`) when the glassy/rubbery segment slopes are nearly parallel and the fitted intersection falls far outside the sampled temperature range, instead of silently reporting a nonsensical extrapolated Tg.
  - Crystallinity (low/medium/high) approximated via initial packing density fraction — documented as trend-only, not a true semi-crystalline lattice.
  - Fixed a degenerate Tg fit (barostat compressibility was water's 4.5e-5 bar⁻¹ instead of a polymer-appropriate 2.0e-6 bar⁻¹, letting the cell drift instead of condense) — tracked for workflow-030/031 in #198.
  - Added mock mode (same pattern as workflow-033): `run.sh` downloads pre-computed `output_files/<node>/` outputs by default; real computation moved to `run_real.sh`. `sample_outputs/` holds a verified reference run's final report.
  - Closes #196.

## [2026-04-08]

### Added

- workflow-016: DiffDock-PP Antibody-Antigen Docking Pipeline (6 nodes: complex-splitting, structure-prep, feature-extraction, inference, analysis, comparison)
  - Comprehensive pipeline for rigid-body docking of antibody-antigen complexes.
  - Aligned with Silva v0.5.1 standards:
    - Implemented `PARAM_` prefix for all runtime environment variables.
    - Leveraged automatic GPU detection (removed manual `use_gpu` flags).
    - Added `params.json` for headless execution support.
  - Performance Optimizations:
    - Pre-downloaded ESM-2 model weights (`esm2_t33_650M_UR50D`) in Docker image to eliminate 2.5-minute runtime download bottleneck.
    - Forced GPU acceleration for inference node.
  - Improved Stability:
    - Fixed fragile `grep` patterns for input resolution in comparison node.
    - Synchronized parameter defaults (10 samples, 20 steps) across all scripts.

## [2026-02-13]

### Added

- workflow-012: Boltz-2 Structure Prediction (4 nodes: sequence-upload, structure-prediction, report, molstar-visualization)
  - Migrated from pre-migration structure_prediction/ directory
  - Reuses boltz_dashboard.py from container-images-for-potter for HTML report generation
  - Includes official prot.yaml example (141 residue protein)
  - GPU-accelerated prediction with Boltz-2
  - Mol* 3D visualization of top-ranked models
- workflow-013: BoltzGen Binder Design (4 nodes: target-upload, binder-design, report, molstar-visualization)
  - New workflow for AI-based protein/peptide binder design
  - Two examples: antibody Fab vs PD-L1 (7uxq), peptide vs BeetleTERT (5cqg)
  - New boltzgen_dashboard.py for design quality assessment
  - GPU-accelerated design with BoltzGen
  - Mol* 3D visualization of top-ranked designs

## [2026-01-06]

### Changed

- workflow-007/03_visualize/generate_output.py: Improved HTML output with full-page responsive layout
  - Wrapped 3Dmol viewer in a complete HTML document structure
  - Full-page viewer container (max-width: 1400px)
  - Added zoom-out script to show protein at 70% for better overview
  - Clean header with title, card-style viewer with shadow
  - Added UI controls for protein representation styles (Cartoon, Stick, Sphere, Line, Surface)
  - Fixed: Surface now properly removed when switching to other representations

## [2025-12-08]

### Added

- workflow-007/example.py: Merged 3 visualization scripts with configurable options (pocket_style, render_method, representation)
- workflow-007/example.py: Tested with pocketeer docker image - successfully generated pocket_visualization.html and rotating GIF
- workflow-007/example.py: Added output_format option ("html", "gif", or "both")
- workflow-007: Split example.py into 3 modular scripts (tested):
  - 01_download/download_pdb.py: Download PDB file
  - 02_pocket/calculate_pockets.py: Calculate pockets and save to JSON
  - 03_visualize/generate_output.py: Generate HTML/GIF visualization
- workflow-007/03_visualize: HTML visualization now uses pocket_style and representation options
  - receptor_cartoon/receptor_surface based on representation setting
  - sphere_scale adjusted based on pocket_style (2.0 for single_sphere, 1.0 for filled_surfaces)
- workflow-007: Added Chiral workflow configuration:
- workflow-003: Migrated from old silva config format to new format:
  - Converted workflow.json to workflow.toml
  - Merged node.json into job.toml for all 6 jobs
  - Updated container.docker_image to container.image
  - Added name, description to job.toml files
  - Removed deprecated node.json and workflow.json files
  - Tested successfully with silva
