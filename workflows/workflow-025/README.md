# Create a silva-runnable workflow for Docking Comparison of Autodock Vina vs gnina from the node structure

Key Files and References: 
- Node structure | '~/dev/collab-workflows/workflows/workflow-025/'
- GitHub Issue | 'https://github.com/chiral-data/collab-workflows/issues/138' 
- Workflow references
    - ‘~/dev/collab-workflows/workflows/workflow-014’
    - '~/dev/collab-workflows/workflows/workflow-018/'
- Silva source code | ‘~/dev/silva’ 
- Silva migration guide | 'https://github.com/chiral-data/collab-workflows/blob/main/SILVA_MIGRATION_GUIDE.md#5-step-3--implementation' (5. Step 3 — Implementation)

Key Tool and Paper References: 
- AutoDock Vina: https://github.com/ccsb-scripps/AutoDock-Vina | pip install vina | Apache 2.0
- GNINA: https://github.com/gnina/gnina | prebuilt binary | Apache 2.0
- Open Babel: https://github.com/openbabel/openbabel | conda install -c conda-forge openbabel | GPL 2.0
- xTB: https://github.com/grimme-lab/xtb | conda install -c conda-forge xtb | LGPL 3.0 
- Buccheri et al., 2025 | 'https://pmc.ncbi.nlm.nih.gov/articles/PMC12388557/'

Docker Files
- AutoDock Vina Dockerfile | 'https://github.com/chiral-data/collab-workflows/blob/main/workflows/workflow-004/Dockerfile'
- Gnina Dockerfile | 'https://github.com/chiral-data/collab-workflows/blob/main/apps/g/gnina_2025_12_04/Dockerfile'

Reference Config File Formats: 
- **workflow.toml**: `workflows/workflow-018/.chiral/workflow.toml`
- **job.toml**: `workflows/workflow-018/00-download/.chiral/job.toml`
- Dependencies are declared centrally in `workflow.toml`, not in individual `job.toml` files
- Parameters use `env = "PARAM_..."` to map to environment variables
- Only nodes with upstream dependencies appear in `[dependencies]` (the first node is omitted)

## Summary of Key Ideas and Goals

### Scientific Context

This workflow compares two molecular docking tools — **AutoDock Vina** and **GNINA** — on a metal-coordinating target, **Carbonic Anhydrase II (PDB: 1OKL)**, which binds Zn²⁺. Both tools share the same MCMC sampling engine (GNINA is a fork of Smina, which forks Vina), so pose generation is effectively held constant. Any difference in ranking, enrichment, or predicted affinity is therefore attributable solely to the **scoring function**: Vina's empirical physics-based energy approximation vs. GNINA's 3D CNN. This makes the comparison scientifically clean and interpretable.

### Central Research Question

*Should a researcher trust the empirical ΔG score or the CNN score for a metal-binding pocket?*

- **Vina's limitation**: Its additive atom-contact terms cannot model coordinate covalent bonds to transition metal ions (Zn²⁺), which skews ΔG predictions and can overscore large, heavy compounds with more pocket contacts.
- **GNINA's advantage/caveat**: The CNN implicitly learns the geometric "picture" of a metal coordination pocket — but can be biased by training data and occasionally score chemically implausible poses highly.

### Key Paper Findings (Buccheri et al., 2025)

For Carbonic Anhydrase II specifically:

| Metric | GNINA | Vina |
|--------|-------|------|
| Pose RMSD | 1.37 Å | 6.78 Å |
| Virtual screening hits | 149 (98% active) | 2 (50% active) |
| Enrichment Factor (EF1%) | 20.75 | 0 |

**Conclusion**: GNINA substantially outperforms Vina for metal-binding pockets; CNN-based rescoring is more effective for drug discovery on targets like CA-II.

### Workflow Goal

Build a reproducible, Silva-runnable pipeline that:
1. Validates inputs and fetches the CA-II receptor from RCSB PDB
2. Prepares structures (receptor → PDBQT, ligands → 3D-optimized via Open Babel + xTB)
3. Runs a redocking QC check with GNINA to confirm the binding pocket is valid before the full screen
4. Runs parallel docking: Vina screening (`04a`) and GNINA screening (`04b`)
5. Generates a side-by-side HTML report with score correlations, ranking histograms, and runtime metrics

### Node Structure Overview

| Node | Tool(s) | Key Output |
|------|---------|------------|
| `01_validate_inputs` | Python, RCSB PDB API | Validated receptor + ligand files, resolution check |
| `02_prepare_structures` | Open Babel, xTB | `receptor.pdbqt`, `optimized_screening_library.pdbqt`, `pocket_config.txt` |
| `03_target_redocking` | GNINA | `qc_validation_results.json` (CNN score + RMSD pass/fail) |
| `04a_run_vina` | AutoDock Vina | `vina_screening_poses.pdbqt`, `vina_runtime_log.txt` |
| `04b_run_gnina` | GNINA | `gnina_screening_poses.pdbqt`, `gnina_runtime_log.txt` |
| `05_generate_report` | Python, Open Babel | `comparative_screening_metrics.csv`, `docking_performance_report.html` |

### Evaluation Criteria

- **Scientific validity**: Redocking QC (node 03) confirms the binding pocket is well-defined before screening
- **Scoring isolation**: Pose generation is identical between tools; observed differences are attributable to the scoring function only
- **Report clarity**: Self-contained HTML report a researcher can interpret without re-running any code
- **End-to-end Silva run**: Full workflow completes locally with all green checkmarks on the 1OKL / Zn²⁺ test case

## Tasks 
- [X] Investigate the GitHub issue and fully the provided read references of files, tools, and papers; create a summary of key ideas and goals. 
- [X] Rewrite the Sample Directory Structure in the README to match the Proposed Node Structure described in the GitHub Issue, keeping a similar format. 
- [] Using existing examples and references, and the sample directory structure, create a workflow by building one node at a time. Make it silva runnable.

**Sample Directory Structure**
```
workflows/workflow-025/
├── .chiral/workflow.toml
├── input_files/
│   └── ligands.smiles
├── 01_validate_inputs/
│   ├── .chiral/job.toml
│   ├── .chiral/test_inputs/sample_ligands.smiles
│   ├── run.sh
│   └── validate_inputs.py
├── 02_prepare_structures/
│   ├── .chiral/job.toml
│   ├── run.sh
│   └── prepare_structures.py
├── 03_target_redocking/
│   ├── .chiral/job.toml
│   ├── run.sh
│   └── redock_target.py
├── 04a_run_vina/
│   ├── .chiral/job.toml
│   ├── run.sh
│   └── run_vina.py
├── 04b_run_gnina/
│   ├── .chiral/job.toml
│   ├── run.sh
│   └── run_gnina.py
├── 05_generate_report/
│   ├── .chiral/job.toml
│   ├── run.sh
│   └── generate_report.py
└── README.md
```

## Outputs