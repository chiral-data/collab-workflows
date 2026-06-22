# Create a silva-runnable workflow for an AI-driven Structure Prediction → Pocket Prediction → Docking pipeline

Key Files and References: 
- Node structure | '~/dev/collab-workflows/workflows/workflow-028/'
- GitHub Issue | 'https://github.com/chiral-data/collab-workflows/issues/172' 
- Workflow references
    - ‘~/dev/collab-workflows/workflows/workflow-025’ (Vina vs GNINA): best reference for a parallel comparison + report workflow
    - '~/dev/collab-workflows/workflows/workflow-014/' (ADMET-AI): reference for SMILES-based property prediction pipeline + Plotly radar charts + Bootstrap 5 layout pattern
    - '~/dev/collab-workflows/workflows/workflow-012/' (Boltz-2): reference for structure prediction workflow + existing workflow-012 image (ghcr.io/chiral-data/boltz:2025_09_05) can be reused for node 01
- Silva migration guide | SILVA_MIGRATION_GUIDE.md in repo root

## Tasks 
- [x] Investigate the GitHub issue and fully the provided read references of files, tools, and papers; create a summary of key ideas and goals. 
- [] Rewrite the Sample Directory Structure in the README to match the Proposed Node Structure described in the GitHub Issue, keeping a similar format. 
- [] Using existing examples and references, and the sample directory structure, create a workflow by building one node at a time. Make it silva runnable.

## Summary of Key Ideas and Goals

This workflow implements **Direction B** from [GitHub Issue #172](https://github.com/chiral-data/collab-workflows/issues/172): an end-to-end AI-driven drug discovery pipeline that goes from a protein sequence directly to docked small molecule poses — no crystal structure required.

### Pipeline Overview

```
protein.yaml (sequence + ligand SMILES + affinity block)
         │
         ▼
01_boltz2_predict ──► structure.cif + plddt_*.npz + pae_*.npz + affinity_*.json
         │
         ▼
02_p2rank_pocket_find (-c alphafold) ──► *_predictions.csv (pocket centers + scores)
         │
         ▼
03_pocket_qc_grid ──► pocket_qc.json (pLDDT distribution per pocket, pass/fail)
                      grid.json (center + computed box size)
                      receptor.pdb (gemmi mmCIF→PDB conversion)
                      ligand.sdf (RDKit SMILES→SDF conversion)
         │
         ▼
04_unimol_docking ──► docked_poses.sdf (3D poses only; no affinity from Uni-Mol)
         │
         ▼
05_generate_report ──► report.html (Molstar 3D viewer + Plotly charts + summary tables)
```

### Tools

**Boltz-2** ([GitHub](https://github.com/jwohlwend/boltz), Wohlwend et al. 2024 / Passaro et al. 2025)
- Diffusion-based 3D structure prediction for proteins, nucleic acids, small molecules, and complexes
- **Input must be YAML** — FASTA mode does not support ligand specification or affinity prediction. Including the ligand SMILES in the YAML enables *holo* prediction, which significantly improves binding site geometry over apo prediction (AlphaFold3 VS paper, biorxiv 2025)
- **Always use mmCIF output** — the `--output_format pdb` flag is broken for protein-ligand complexes ([boltz #298](https://github.com/jwohlwend/boltz/issues/298), IndexError). P2Rank natively accepts mmCIF, so no conversion is needed at this step
- Outputs: `*.cif`, `plddt_*.npz` (per-residue 0–1), `pae_*.npz`, `confidence_*.json`, `affinity_*.json`
- Existing Docker image `ghcr.io/chiral-data/boltz:2025_09_05` (from workflow-012) can be reused

**P2Rank** ([GitHub](https://github.com/rdk/p2rank), Jakubec et al. 2025 *NAR*)
- ML-based ligand binding pocket prediction directly from 3D protein structure
- **Must use `-c alphafold` profile** when input comes from Boltz-2 — the B-factor column contains pLDDT, not thermal displacement, and this profile handles that correctly
- Outputs `*_predictions.csv` with `center_x, center_y, center_z` (SAS point centroid in Å) and `score`, `probability` per pocket
- Does **not** output box size — node 03 must compute `size_x/y/z` from pocket spatial extent (pocket SAS radius + 5–10 Å padding → 15–25 Å box)
- Requires Java runtime in the container; CPU-only

**Uni-Mol Docking V2** ([GitHub](https://github.com/deepmodeling/Uni-Mol), Zhou et al. 2022 / ICLR 2023)
- Universal 3D molecular representation model applied to protein-ligand docking; achieves 77.6% RMSD < 2 Å
- Inputs: receptor **PDB** (not mmCIF) + ligand **SDF** + docking grid JSON with `center_x/y/z` and `size_x/y/z`
- **Outputs only 3D poses (SDF) — no binding affinity or confidence score.** The internal `prmsd_score` used for pose ranking is not written to the output file
- Known issue: RDKit sanitization can silently drop invalid ligands ([Uni-Mol #281](https://github.com/deepmodeling/Uni-Mol/issues/281)) — add a validation check in node 03 before passing to Uni-Mol
- Containerization: requires Uni-Core (custom CUDA C++ extensions); base image `dptechnology/unicore:0.0.1-pytorch1.11.0-cuda11.3`; model weights (464 MB) must be downloaded manually from Dropbox (no pip/HuggingFace distribution)

### Node 03 Conversions (QC + Grid)

Node 03 is the integration hub and handles several format conversions before docking:

| Conversion | Tool | Notes |
|------------|------|-------|
| mmCIF → PDB | `gemmi convert` | Uni-Mol requires PDB; avoid BioPython's `MMCIFIO` (known incompatibility with Boltz-2 CIF's missing `entity.id`) |
| SMILES → SDF | RDKit `AllChem.EmbedMolecule` | Extract SMILES from input YAML; validate before passing to Uni-Mol |
| Pocket center → grid JSON | Python | P2Rank outputs only center; compute box size as pocket SAS radius + padding |
| pLDDT QC | `plddt_*.npz` + P2Rank residue list | Read per-residue values from npz (0–1 scale); report mean, min, std per pocket |

### pLDDT Threshold (≥ 70)

The ≥ 70 threshold is well-supported: AlphaFold's official scale defines ≥ 70 as "confident, backbone correct"; PrankWeb 4 uses it as the default filtering toggle; a Nature Communications 2024 druggable-pocket study required all pocket-centroid residues within 8 Å to have pLDDT > 70.

**Caveat:** Eguida & Rognan 2023 (*JCIM*, PMC9852548) tested 22 drug targets and found 4 of the 5 worst-performing docking targets had binding-site pLDDT ≥ 70 — high pLDDT is necessary but not sufficient. Node 03 should report the **per-pocket pLDDT distribution** (mean, min, std across pocket residues) rather than a binary pass/fail, so the user can make an informed judgment.

### Report Node (05)

The report is a full HTML dashboard — not just a 3D viewer — combining Molstar, Plotly, and Bootstrap 5 (following workflow-014's layout pattern). It includes:

1. **Pipeline summary** — input sequence length, model count, pocket count, QC-passed pockets, docked compounds
2. **Structure confidence panel** — Boltz-2 pLDDT distribution histogram + PAE heatmap (Plotly)
3. **Pocket discovery table** — P2Rank rank, score, probability, mean/min pLDDT per pocket, center coordinates (sortable)
4. **3D interactive viewer** — Molstar: protein colored by pLDDT, P2Rank pocket highlighted, top Uni-Mol docking pose overlaid
5. **Methods & caveats** — auto-generated pipeline description with "pLDDT ≥ 70 is necessary but not sufficient" note

Since Uni-Mol outputs no affinity score, the report uses **Boltz-2's own confidence metrics** (`plddt_*.npz`, `pae_*.npz`, `affinity_*.json`) combined with P2Rank's `score` and `probability` as the confidence signal.

Refs: dockviz (ljmartin, 2024) for PDBe-Molstar embedded HTML; PrankWeb 4 for pocket-colored-by-pLDDT pattern; workflow-014 for Plotly + Bootstrap 5 layout.

### Scientific Motivation

Traditional virtual screening requires a solved crystal structure. This pipeline removes that dependency: Boltz-2 predicts the 3D fold from sequence alone (with the ligand included for holo-mode accuracy), P2Rank finds the druggable pocket using an AlphaFold-aware profile, and Uni-Mol Docking V2 generates docked poses. pLDDT confidence propagates from structure prediction through pocket QC to the report, flagging results from uncertain structural regions.

For validation, the example input should be a target with an **existing holo crystal structure** (≤ 2.5 Å, co-crystallized ligand with known SMILES). This lets us check whether Boltz-2's predicted structure matches the crystal, P2Rank identifies the correct pocket, and Uni-Mol's pose is close to the experimental binding mode.

### Patterns from Reference Workflows

| Reference | Pattern to reuse |
|-----------|-----------------|
| **workflow-012** (Boltz-2) | Reuse `ghcr.io/chiral-data/boltz:2025_09_05`; YAML input format; pLDDT/PAE npz parsing |
| **workflow-025** (Vina vs GNINA) | HTML report structure; Bootstrap 5 + Plotly layout |
| **workflow-014** (ADMET-AI) | Plotly chart patterns; Bootstrap 5 dashboard; SMILES validation via RDKit |

### Key Implementation Challenges

1. **Uni-Core containerization** — requires source compilation against CUDA 11.3; base image `dptechnology/unicore:0.0.1-pytorch1.11.0-cuda11.3`; 464 MB model weights need manual download from Dropbox
2. **Box size computation** — P2Rank outputs only pocket center coordinates; node 03 must derive `size_x/y/z` from pocket spatial extent
3. **Format conversions in node 03** — mmCIF → PDB (`gemmi`), SMILES → SDF (RDKit), with ligand validation before Uni-Mol to avoid silent drops
4. **GPU requirements** — nodes 01 (Boltz-2) and 04 (Uni-Mol) require GPU; node 02 (P2Rank) is CPU-only; `job.toml` container specs must reflect this
5. **No affinity from Uni-Mol** — the report must be built from Boltz-2 metrics and P2Rank scores, not Uni-Mol output

**Sample Directory Structure**
```
workflows/workflow-028/
├── .chiral/workflow.toml
├── input_files/sample_molecules.csv
├── 01_validate_inputs/
│   ├── .chiral/job.toml
│   ├── .chiral/test_inputs/sample_molecules.csv
│   ├── run.sh
│   └── validate.py
├── 02_compute/
│   ├── .chiral/job.toml
│   ├── pre_run.sh          
│   ├── run.sh
│   └── compute_admet.py
├── 03_analyze/
│   ├── .chiral/job.toml
│   ├── run.sh
│   └── analyze.py
├── 04_visualize/
│   ├── .chiral/job.toml
│   ├── run.sh
│   └── generate_report.py
└── README.md
```

## Outputs