# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2026-11-20]

### Added

- workflow-038: BioEmu2 Equilibrium Ensemble Emulation — run end to end on a GPU
  (6 nodes: 01-prepare-system, 02-fetch-msa, 03-sample, 04-propagate, 05-analyze, 06-report).
  - Conformational ensemble emulation on BioEmu2 (Microsoft Research AI for Science, MIT,
    public preview 2026-10-20). First workflow in the repo that samples how an existing
    protein *moves* rather than predicting one structure or designing a new one.
  - Covers all three BioEmu2 modes by parameter: `predict-structure`, `sample-equilibrium`,
    and time-lagged `propagate`.
  - `02-fetch-msa` is its own CPU node: the alignment is the only step needing the public
    internet, so fetching it to disk turns a server outage or rate limit into a seconds-long
    CPU failure instead of a failed GPU job, and lets `03-sample` run with `HF_HUB_OFFLINE=1`.
    User-supplied `<chain>.a3m` files skip the server entirely — today only when the node is
    run directly with Docker, since Silva's one file-typed parameter routes to node 01 (open
    item 7); `msa_depth` does not apply to an uploaded alignment and the node says so rather
    than appearing to subsample it. The alignment moves between
    nodes as the directory `msa`, declared by name rather than as `msa/*.a3m`: Silva matches an
    `inputs` entry by basename, so the pattern staged nothing and nodes 03/04 sampled with no
    alignment at all. Both now refuse to run when node 02 reported protein chains and no A3M
    arrived, since `EmuSession` otherwise fetches its own over the network and the manifest
    would still name node 02's. A heteromer's **positionally paired** alignment is written and
    reattached as well (`msa/<chain>.paired.a3m`): `load_a3m_msa` always returns an unpaired
    MSA, and BioEmu2's pairmixer substitutes zeros for a missing paired feature, so carrying
    only the unpaired half degraded every multi-sequence complex without failing.
  - Nodes 03/04 drive the BioEmu2 Python API rather than its CLI, because the CLI exposes no
    flag for a pre-computed MSA. `num_steps` *is* exposed: `EmuSession.sample()` takes
    `sampler: SamplerPreset | Sampler`, so `sampler_from_preset(preset, mode=…, num_steps=…)`
    reaches the step count. A default run still passes the string preset.
  - Node 01's `system_inventory.json` travels the length of the pipeline, so the report states
    the chains and token count node 01 measured rather than re-reading them from the topology.
    Node 03 prints BioEmu2's own largest-fitting-batch estimate from that token count and warns
    when the batch exceeds it, without applying it — `batch_size` changes which samples are
    drawn, so deriving it from the card would break reproducibility across GPUs.
  - `04-propagate` always writes its declared outputs, passing the ensemble through as an
    identity propagation when disabled, because Silva has no conditionally-skipped node and a
    declared-but-unwritten output is a silent failure.
  - Nodes 04-05 bind the XTC to the System's own topology rather than parsing `topology.cif`,
    sidestepping mdtraj's mmCIF support entirely.
  - Shared image `bioemu2:2026_11_20` (`python:3.12-slim`, 9.58 GB; CUDA 12.6 via the
    `[cuda12]` extra because the `[cuda]` extra needs a CUDA 13 driver the estate lacks and
    pins no torch version at all). Compute-capability floor **measured** rather than inherited:
    `torch.cuda.get_arch_list()` in the built image returns `sm_50 … sm_90`, so `sm_70` is
    present (V100 works) and `sm_100` is absent (not B200). HF checkpoint pinned to
    `dd4b0c1b69984b68edcb6200a544213be4d75242` (1.91 GiB); measured in the image: openmm 8.6.1,
    mdtraj 1.11.2, torch 2.13.0+cu126, numpy 2.4.6, biotite 2.0.0. Registered in
    `apps/applications.json`.
  - The image also pre-warms **126 chemical component definitions** into `BIOEMU2_CCD_CACHE`
    and `BIOEMU2_TEMPLATE_CACHE`. Anything that is not a standard residue — an ion, a cofactor,
    a terminal cap, a modified residue — is built from its CCD entry, which bioemu2 downloads
    from `files.rcsb.org` the first time it meets the code, from inside `system.topology`: a
    call nodes 01, 03, 04 and 05 all make. Five of the seven test entries needed one (1CLL `CA`
    and `EOH`, 1A3N `HEM`, 1AKE `AP5`, 5PTI `PO4`/`UNX`/`DOD`, 2ZTA `ACE`), and node 01 used to
    swallow the failure into a warning and write `n_tokens: null`, which also made the
    batch-size advice and the >10k GPU-fit guard unreachable on every holo system. Node 01 now
    fails there — nodes 03-05 build the same structure — and every node that builds a topology
    names the component and the cache in the error.
  - Node 01 exposes upstream's `include_ligands` and `exclude_ligands`, lists every non-polymer
    component it kept by CCD code in `system_inventory.json`, and warns about two classes worth
    a second look: placeholders (`UNX`, `UNK`, `UNL`, `DUM` — density the depositor could not
    identify) and solvent that upstream does not treat as solvent (`DOD`, the deuterated water
    of a neutron structure; `System.from_cif` omits `HOH` and keeps `DOD`). `n_polymer_residues`
    is labelled for what upstream documents it to be — the entity sequence, including residues
    the input structure does not resolve, so 1CLL is 148 where its crystal resolves 144.
  - Node 03's `--dry-run` is a gate rather than a parameter echo: it verifies the alignment was
    staged, resolves the topology, and compares each staged A3M's query row against the sequence
    the system states for that chain. bioemu2 compares neither — not `get_msa`, not
    `load_a3m_msa` — so an alignment for a different protein was accepted and produced an
    ensemble whose frame count, physicality rate and geometry were all unremarkable.
  - Node 05's collapse warning fired below `1e-3` nm, which an XTC's 0.001 nm coordinate
    quantisation makes unreachable: one conformation written twice and read back measures about
    0.002 nm apart. Now 0.01 nm, with the measurement recorded next to the constant.
  - A single frame is no longer reported as a collapsed ensemble. `mode = structure` with
    `num_samples = 1` is the configuration the README recommends for predicting one folded
    structure; node 05 reported `max_nm: 0.0` over zero pairs for it — the collapse signature —
    and then warned that the ensemble had collapsed and to check the mode. `pairwise_rmsd` now
    carries `n_pairs: 0` with null statistics, node 05 states that a spread needs two frames,
    and node 06's tile reads *n/a, needs 2+ frames*. Measured on a GPU: 1UBQ returns one frame
    at Rg 1.1616 nm and the bundled Trp-cage 1L2Y one at 0.7146 nm against 0.697 expected,
    both physical under `high_quality` checking. One frame from an *equilibrium* run is still
    a wasted run, and the spurious warning was the only thing that had said so, so node 03
    warns before the GPU time is spent and node 05's message is keyed off the mode node 03
    recorded. Two fixtures hold both halves, because a node that said nothing in either case
    would pass a single-mode fixture.
  - `tests/`: a fixture suite over seven real RCSB entries (`tests/README.md`), split into an
    offline driver (`run_fixtures.sh`: unit checks, nodes 01-02 plus node 03 `--dry-run`,
    refusals and their counter-cases, and nodes 04-06 over deterministic ensembles) and a GPU
    campaign (`run_gpu_campaign.sh`) that runs the whole pipeline for real against the MSA
    server. Every check compares what a node *reported* against an independent measurement of
    the input: `probe_entry.py` parses `_atom_site` by hand, so it shares no parser with the
    code under test. 1CLL is there because a calcium ion is an atom named `CA` and mdtraj's
    `name CA` matches atom names, which silently puts 148 atoms in a 144-C-α selection. The
    campaign driver records each node's exit status next to its log, and the checker uses it to
    separate a node that ran and produced nothing from one that was never reached: a system too
    large for the card reports its downstream checks as skipped-with-the-reason, and a separate
    check requires that the node which stopped said why before it did. A `FAIL` therefore always
    means a node that ran and got something wrong. Against `bioemu2:2026_11_20` the offline
    suite reports 70 pass / 0 fail / 0 gap / 2 skipped and the campaign 86 pass / 0 fail /
    42 skipped; two successive campaigns returned every measured number identical, so seed 42
    reproduces a run rather than merely labelling one. Two findings were in the suite itself:
    `make_ensembles.py` seeded from `hash()`, which Python salts per process, so its fixtures
    were redrawn on every run while claiming to be byte-identical; and the medoid check
    required `rmsd_to_medoid_nm` to be zero at the medoid, which mdtraj's float32 RMSD never
    returns (0.0013-0.0018 nm for a frame against itself). Together they made a check that
    passed or failed on the dice. Both are fixed, and the check now asserts the medoid is the
    argmin with its own distance below that floor.
  - Coverage gaps closed in the same pass, each of them a parameter the workflow advertises
    and no test had ever set: `exclude_ligands` and `include_ligands` (two 5PTI fixtures,
    measured — 65 tokens with `PO4`/`UNX`/`DOD`, 63 with the placeholder and the heavy water
    dropped, 58 with every component dropped, which is the polymer residue count exactly),
    the combination of the two that upstream forbids (a refusal naming which parameter does
    what), the raw `sequence` input, `mode = structure`, and `rmsd_cluster_threshold_nm` —
    the parameter behind the report's most quotable number, checked by re-deriving the
    cluster count from the published matrix under node 05's own greedy rule, tie-break
    included. `pdb_id` needs the network and
    was measured by hand: 1UBQ resolves to `RCSB:1UBQ` at 76 residues and 76 tokens.
    `msa_depth` and `use_env` remain unexercised and `tests/README.md` says so.
  - The image also installs `gcc` and `libc6-dev`, which are a **runtime** requirement: triton
    compiles its kernel launchers on first use and shells out to a C compiler, and
    `python:3.12-slim` carries neither. 1UBQ (76 tokens), 5PTI (65) and 2ZTA (72) never reached
    that path; 1CLL (155) failed inside `EmuSession.sample()` with "Failed to find C compiler",
    after the model had loaded, which reads as a size limit or a bad input rather than a missing
    package.
  - Measured on an RTX 3060 Ti (8 GB, sm_86, driver 580.178.04, 7.66 GiB usable), 8 samples per
    system in one batch: 1UBQ 31.2 s with 8/8 frames physical, Rg 1.1597 ± 0.0056 nm against an
    experimental 1.17 nm and max pairwise C-α RMSD 0.2204 nm; 5PTI Rg 1.1214 nm; 2ZTA Rg
    1.5244 nm over 4 clusters (a parallel coiled coil, so 1.39× the globular scaling is the fold
    rather than a failure); 1BNA, which has no protein chain, through all six nodes on the
    all-atom fallback; and 1CLL Rg 1.5677 nm with max pairwise RMSD 1.3682 nm over 8 clusters on
    a 148-C-α selection that excludes its four calcium ions — calmodulin's lobes pivot about a
    flexible central helix, so a broad multi-state ensemble is the literature answer, and the
    contrast with ubiquitin's single 0.22 nm cluster is the nearest thing here to a scientific
    control. 1AKE (542 tokens) and 1A3N (746) run out of memory, with node 03's own OOM message
    — and re-run at batch 1 they fail on the same allocation as at batch 8, 1.06 and 1.46 GiB,
    so on this card the peak is set by the token count and lowering the batch moves nothing.
    Upstream's batch estimate, which node 03 prints, had said 33 and 19 would fit; it is now
    labelled an upper bound rather than a ceiling, and the node's OOM advice says to try a
    lower batch but not to expect much and names the token count instead. The README states
    the ceiling as a size range — 8 GB holds somewhere above 155 and below 542 tokens — because
    that is the variable that binds. Node 02's paired-alignment path is confirmed against the
    live server: 1A3N recorded `paired: true` for its two unique sequences, with unpaired depths
    of 4522 and 5023 rows per chain and a paired alignment of 8794 rows identical in depth across
    all four chains, which is what a positionally paired alignment has to be, and node 03's log
    shows it attached to all four chains before the allocation failed — so no heteromer ensemble
    exists on this hardware, which `tests/README.md` records as the one open gap.
  - Time-lagged propagation run for real for the first time, on 1UBQ at the default 1 ns lag:
    8 endpoints from 8 start points, all 8 physical, mean Rg 1.1628 nm against 1.1597 nm before
    propagation and spread 0.1943 nm against 0.2204 nm. A rigid fold moving slightly and staying
    folded over a nanosecond is the expected answer, so this is strong evidence about the
    plumbing and weak evidence about the dynamics; node 04's **Scientific notes** say so.
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
