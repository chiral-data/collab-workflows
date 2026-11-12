---
doc_id: workflow-038
domain: structural-biology
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  Conformational ensemble emulation with BioEmu2. Takes a sequence, structure or
  RCSB accession, fetches an MSA, and samples an equilibrium ensemble, a folded
  structure, or time-lagged endpoints on one GPU.
tags: [structural-biology, conformational-ensemble, molecular-dynamics, generative-model, gpu]
---

# Workflow 038: BioEmu2 Equilibrium Ensemble Emulation

BioEmu2 ("biomolecular emulator", Microsoft Research AI for Science, public preview 2026-10-20,
MIT) is a ~500M-parameter all-atom diffusion model that emits **many** conformations of a
biomolecular system rather than one. This workflow takes a sequence, a structure or an RCSB
accession, fetches its alignment, samples an ensemble on one GPU, and reports how the molecule
moves.

It answers "what conformations does this protein adopt?" — a question nothing else in this repo
answers. workflow-026 (ColabFold) and workflow-012 (Boltz-2) each predict a single structure;
workflows 030–036 run real MD, but on polymers and at a cost of hours per system.

> **This is an emulator, not a simulator.** Sampling frequency is **not** a calibrated Boltzmann
> probability, so frame counts are not populations and never free energies. Sampled frames are
> **independent draws, not a trajectory**. See [Interpreting the output](#outputs-and-interpretation).

## Overview

BioEmu2 is a diffusion-based all-atom score model conditioned on a multiple sequence alignment,
with optional contact and method conditioning. It was trained for two distinct tasks —
**structure prediction** (one accurate folded structure, sampled with heuristic sharpening at
step scale 1.5) and **equilibrium sampling** (a range of conformations from the modeled
equilibrium ensemble, at step scale 1.0 unsharpened) — and additionally supports **time-lagged
propagation** through implicit transfer operator models, which draws one endpoint a chosen lag
time ahead of each starting conformation.

A "conformation" here is a full all-atom structure: the model handles proteins, nucleic acids,
ligands and multi-chain complexes, though it does no de novo ligand design and binding-affinity
prediction is explicitly out of scope. The ensemble this workflow produces is the model's
estimate of the accessible conformational space, conditioned on the alignment it was given —
which is why alignment depth is reported everywhere and why a shallow MSA is treated as a
caveat rather than a detail.

## When to use this workflow

Use it when you have a sequence or structure and want to know its **conformational
heterogeneity**: whether it has distinct states, which regions are flexible, how compact it is
across the ensemble. It is well suited to screening for conformational change before committing
to an MD campaign, to finding flexible loops, and to generating starting structures for
downstream work. It needs one NVIDIA GPU of compute capability ≥ 7.0 (V100 or newer).

**Do not use it when:** you want a single best structure with per-residue confidence — use
workflow-026 (ColabFold), which reports pLDDT and PAE; this workflow emits no confidence
estimate of any kind. You want free energies, populations or binding affinities — BioEmu2's
frequencies are uncalibrated and affinity is out of scope; use a real MD or FEP pipeline. You
want to design a new protein — use workflow-033, -037 or -013. Your system exceeds ~10k tokens
— upstream notes physically plausible structures become hard above that, and the
context-parallel (≥4 GPU) path is out of scope here. Your protein is an orphan with a shallow
alignment — the ensemble will be poor, and the report will tell you so, but it cannot fix it.

## Architecture and data flow

```text
sequence / structure / system YAML / RCSB accession
   │
   ▼
[01-prepare-system] ──────► system.yaml, system_inventory.json
   │
   ▼
[02-fetch-msa]  (network) ─► msa/*.a3m, msa_summary.json
   │
   ▼
[03-sample]     (GPU) ─────► samples.xtc, topology.cif, physicality.json
   │
   ▼
[04-propagate]  (GPU) ─────► propagated_samples.xtc, propagation_manifest.json
   │
   ▼
[05-analyze] ──────────────► ensemble_metrics.csv, ensemble_summary.json, ensemble_top.pdb
   │
   ▼
[06-report] ───────────────► report.html
```

Nodes run sequentially 01 → 06. Only nodes 03 and 04 need a GPU. Node 02 is the only node that
needs network access; it writes the alignment to disk so node 03 runs fully offline with
`HF_HUB_OFFLINE=1` against the checkpoint baked into the image.

## Input requirements

One system, supplied any of these ways (first match wins):

1. **A file** in the workflow's input directory — a BioEmu2 system YAML (`.yaml`/`.yml`, Boltz
   YAML auto-detected and converted), an mmCIF (`.cif`), a PDB (`.pdb`), or a FASTA
   (`.fasta`/`.fa`/`.faa`).
2. **`pdb_id`** — an RCSB accession, loaded as biological assembly 1 including all non-water
   ligands.
3. **`sequence`** — a raw one-letter protein sequence.

Node 02 also accepts pre-computed `<chain>.a3m` alignment files placed in its input directory,
in which case it contacts no server at all. **This path is currently reachable only when running
node 02 directly with Docker** (mount your A3M files into `/workspace/inputs`); no workflow
parameter routes an uploaded alignment to node 02 yet. See open item 7 in
[bioemu2_pipeline.md](bioemu2_pipeline.md).

See [input_files/README.md](input_files/README.md) for the bundled Trp-cage test system.

## Workflow nodes

### Node 01: Prepare System

**Goal:** One canonical BioEmu2 system YAML, whatever the user supplied.

**Process:** Detects the input by extension and funnels everything through `System` —
`from_yaml`, `from_cif`, `from_description` or `from_pdb_id` — then writes it back out with
`to_yaml()` and loads it again as a validation step. PDB input is converted to mmCIF first;
Boltz YAML is converted with `bioemu2 convert-boltz`.

**Scientific notes:** The inventory records the **token** count, not just residues: BioEmu2's
memory and compute both scale with tokens (one per polymer residue, one per ligand atom), so
that is the number that predicts whether a run fits on the GPU. Over 10k tokens the node warns.

**Outputs:** `system.yaml`, `system_inventory.json`

### Node 02: Fetch MSA

**Goal:** The alignment on disk, so no later node needs the network.

**Process:** Calls `get_msa()` for the system, which retrieves an unpaired A3M per protein chain
from the public ColabFold server on a cache miss, optionally subsamples it to `msa_depth` rows,
and writes `msa/<chain>.a3m`. Pre-supplied A3M files in the node's input directory short-circuit
the whole step (Docker-only today — see open item 7).

**Scientific notes:** Alignment depth is the single strongest predictor of ensemble quality —
upstream names shallow or poor MSAs (orphan proteins, some viral proteins) as the first cause of
degradation. `msa_summary.json` records the depth per chain and the node warns below 32 rows.

**Outputs:** `system.yaml`, `msa_summary.json`, `msa/*.a3m`

### Node 03: Sample Conformations

**Goal:** `num_samples` conformations, with a physicality verdict per frame.

**Process:** Builds an `EmuSession` from the system, the baked-in `bioemu2-preview` checkpoint
and node 02's alignment, samples, runs `check_physicality()`, and writes everything with
`save_frames()`. Uses the **Python API rather than the CLI** because the CLI has no flag for
supplying a pre-computed alignment. `--dry-run` resolves and prints every setting without
touching the GPU.

**Scientific notes:** `mode` selects the task the model was trained for — `structure` sharpens
toward one accurate fold, `equilibrium` samples the ensemble. Physicality is **reported, never
enforced**: rejected frames stay in the XTC by design, as upstream intends, and the pass rate is
surfaced instead.

**Outputs:** `samples.xtc`, `topology.cif`, `physicality.json`, `system.yaml`,
`sampling_manifest.json`, `msa/*.a3m`

### Node 04: Propagate Ensemble

**Goal:** One time-lagged endpoint per sampled conformation.

**Process:** When `propagate` is true, loads the ensemble against the system's own topology and
calls `EmuSession.propagate(lag_time_ns=…)` in equilibrium mode. When false, it copies the
ensemble through unchanged and records `enabled: false` — an identity propagation, so the
declared outputs always exist (a declared-but-unwritten output is a silent failure in Silva).
Nodes 05 and 06 read that flag and never present identity frames as time-lagged.

**Scientific notes:** The result is one endpoint per start point at the chosen lag, **not** a
continuous path. Upstream recommends lags of 1–1000 ns, the range the model was trained for.

**Outputs:** `propagated_samples.xtc`, `propagated_topology.cif`, `propagation_manifest.json`

### Node 05: Analyze Ensemble

**Goal:** The ensemble as numbers.

**Process:** Binds the XTC to the system's own topology (no mmCIF parsing), then computes
per-frame radius of gyration and C-α RMSD to the ensemble medoid, a pairwise RMSD matrix,
per-residue RMSF about the mean structure, and a cluster count by greedy medoid assignment at
`rmsd_cluster_threshold_nm`. Writes a multi-model PDB of cluster centres for the viewer.

**Scientific notes:** Every quantity is geometric. Cluster sizes are **frame counts**, not
populations, and this node deliberately computes no thermodynamic quantity. Above 200 frames the
O(n²) matrix and clustering are skipped and RMSD is measured against frame 0 instead.

**Outputs:** `ensemble_metrics.csv`, `ensemble_summary.json`, `ensemble_top.pdb`

### Node 06: Report

**Goal:** One self-contained HTML report.

**Process:** Embeds all data inline as JavaScript and renders stat tiles, a per-frame Rg line
chart, an Rg histogram, a pairwise-RMSD heatmap on a single sequential hue, a per-residue RMSF
line chart, a per-frame table and a 3Dmol viewer that animates the representative frames.
Bootstrap 5, Plotly and 3Dmol come from CDN; no local assets.

**Scientific notes:** The report **leads** with the interpretation caveats rather than burying
them, and shows the MSA depth beside the results so a narrow ensemble can be attributed.

**Outputs:** `report.html`

## Parameters

### mode

| Value / Range | Description |
|---------------|-------------|
| `equilibrium` (default) | Samples a range of conformations from the modeled equilibrium ensemble (step scale 1.0, no sharpening). The reason to use this workflow. |
| `structure` | One accurate folded structure (step scale 1.5, heuristic sharpening). |

**Trade-off:** `structure` gives a better single answer and a worse picture of flexibility;
`equilibrium` the reverse. They are different trained tasks, not quality settings.

**Test vs production:** Leave at `equilibrium`. Use `structure` with `num_samples = 1` as the
cheapest end-to-end correctness check.

### num_samples

| Value / Range | Description |
|---------------|-------------|
| 1 | Only sensible with `mode = structure`. |
| 8 (default) | Enough to see whether the ensemble is heterogeneous at all. |
| 100–1000+ | Production. Required before any cluster count means anything. |

**Trade-off:** Runtime scales roughly linearly in `num_samples` once the trunk embedding is
paid for; peak memory scales with `batch_size`, not `num_samples`.

**Test vs production:** 8 for a smoke test, 100+ for a real ensemble.

### batch_size

| Value / Range | Description |
|---------------|-------------|
| `0` (default) | Use the BioEmu2 default, which denoises all `num_samples` together. |
| 1–N | Denoise N at a time. Lower it when the GPU runs out of memory. |

**Trade-off:** Larger batches are faster up to a saturation point, then only use more memory.
BioEmu2 runs `ceil(num_samples / batch_size)` batches. Changing it changes the exact samples
drawn, and for `edm_fk`/`edm_fkc` also changes the resampling population.

**Test vs production:** `0` until you hit an out-of-memory error, then halve until it fits.

### sampler

| Value / Range | Description |
|---------------|-------------|
| `edm` (default) | Ordinary, unsteered EDM sampling. |
| `edm_gd` | Gradient guidance toward the default physicality potentials. |
| `edm_fk` | Adds Feynman-Kac resampling. |
| `edm_fkc` | Feynman-Kac corrector. |

**Trade-off:** The steered samplers raise the physicality pass rate at the cost of time, and
they tilt the sampled distribution — which is the point, and also a reason not to read the
result as an unbiased ensemble.

**Test vs production:** `edm` unless the physicality pass rate is poor.

### physicality_filter

| Value / Range | Description |
|---------------|-------------|
| `high_quality` (default) | Upstream's default check set for bond lengths, clashes, chirality and chain breaks. |
| `high_quality_per_residue` | The same checks, applied per residue. |
| `no_violations` | No violation tolerated. |
| `no_filter` | Skip the checks; every frame is recorded as `not_checked`. |

**Trade-off:** None in output content — **frames are never removed**, whatever the setting. This
only changes the verdict recorded per frame.

### propagate / lag_time_ns

| Value / Range | Description |
|---------------|-------------|
| `propagate = false` (default) | Node 04 passes the ensemble through unchanged. |
| `propagate = true` | Draw one endpoint per frame, `lag_time_ns` ahead. |
| `lag_time_ns = 1.0` (default) | Upstream recommends 1–1000 ns, the trained range. |

**Trade-off:** Propagation roughly doubles GPU time. Outside 1–1000 ns the model is being asked
for something it was not trained on, and the node warns.

### msa_depth / use_env

| Value / Range | Description |
|---------------|-------------|
| `msa_depth = 0` (default) | Keep the full alignment. |
| `msa_depth = N` | Subsample to N rows. |
| `use_env = true` (default) | Include the environmental database in the ColabFold search. |

**Trade-off:** Shallower alignments embed faster and use less memory (trunk cost scales with
rows × tokens²) but degrade ensemble quality. Subsampling is a deliberate experiment, not an
optimisation.

### precision / matmul_precision / compile_denoiser

| Value / Range | Description |
|---------------|-------------|
| `precision = bf16-mixed` (default) | Trunk embedding in mixed precision; roughly halves its memory and doubles its speed. |
| `precision = 32-true` | Full precision throughout. |
| `matmul_precision = high` (default) | TF32 tensor cores during denoising. `highest` or `medium` also available. |
| `compile_denoiser = false` (default) | Only worth enabling for many samples on small systems (<128 tokens); compilation itself costs time. |

**Trade-off:** All three trade a little numerical precision for speed and memory. Results change
slightly.

### rmsd_cluster_threshold_nm / top_n

| Value / Range | Description |
|---------------|-------------|
| `rmsd_cluster_threshold_nm = 0.2` (default) | C-α RMSD below which two frames count as the same conformation. |
| `top_n = 10` (default) | Frames written to the multi-model PDB and shown in the viewer. |

**Trade-off:** The cluster count is entirely a function of this threshold; it is a descriptive
choice, not a measurement. Report it alongside any cluster count — the report does.

### seed

`42` by default; batch *i* uses `seed + i`. Set to `-1` for unseeded sampling, which is what
upstream does by default.

## Outputs and interpretation

### samples.xtc + topology.cif

The ensemble. Coordinates live in the XTC; the topology is the CIF. **Bind the XTC to the
system, not to the CIF**, if you load it yourself:

```python
import mdtraj
from bioemu2 import System
system = System.from_yaml("system.yaml")
traj = mdtraj.load_xtc("samples.xtc", top=system.topology)
```

Each frame is an independent all-atom conformation. There is no time ordering, and no frame is
"the" answer.

### physicality.json

One entry per frame: zero-based `frame`, `status` (`passed` / `rejected` / `not_checked`) and
the first failed check in `reason`. **Rejected frames remain in the XTC** — that is upstream's
design, not a bug. Filter with `samples[physicality.passed]` if you need only clean frames. No
nucleic-acid or ligand physicality checks are performed.

### ensemble_metrics.csv

One row per frame per ensemble: `radius_of_gyration_nm` (compactness — a flat series across an
equilibrium run means the ensemble collapsed), `rmsd_to_medoid_nm` (how unusual the frame is
relative to the ensemble centre), and the physicality verdict.

### ensemble_summary.json

Ensemble-level: Rg mean/std/range, pairwise RMSD statistics and the full matrix (up to 200
frames), per-residue RMSF, the medoid frame, and clusters at the chosen threshold. **Cluster
sizes are frame counts.** They are not populations, not probabilities, and not free energies,
because BioEmu2's sampling frequency is not a calibrated Boltzmann probability.

### report.html

Self-contained. Opens in any browser with no local assets.

## Quick start

### Running with Docker

Built from [Dockerfile](Dockerfile) as `bioemu2:2026_11_11`.

```bash
docker build -t bioemu2:2026_11_11 workflows/workflow-038/
```

Smoke-test the image before running any node:

```bash
docker run --rm bioemu2:2026_11_11 bioemu2 --help
```

```bash
docker run --rm --gpus all bioemu2:2026_11_11 nvidia-smi --query-gpu=name,compute_cap --format=csv
```

Compute capability must be 7.0 or newer. Then run a node, node 03 shown:

```bash
docker run --rm --gpus all -v /tmp/wf038/inputs:/workspace/inputs -v /tmp/wf038/outputs:/workspace/outputs -v $(pwd)/workflows/workflow-038/03-sample:/workspace -w /workspace -e PARAM_MODE=equilibrium -e PARAM_NUM_SAMPLES=8 -e PARAM_SEED=42 bioemu2:2026_11_11 bash run.sh --dry-run
```

Drop `--dry-run` to actually sample. Both GPU nodes support it.

### Running on Silva

1. Select this workflow from the workflow list
2. Upload a system YAML, structure or FASTA file (or set `pdb_id`, or `sequence`)
3. Adjust parameters if needed (see Parameters)
4. Click Run

### Test vs production settings

| Setting | Test (default) | Production |
|---------|---------------|------------|
| `num_samples` | 8 | 100–1000+ |
| `batch_size` | 0 (auto) | tuned to the GPU |
| `propagate` | `false` | `true` if you want dynamics |
| input | bundled Trp-cage (20 residues, seconds) | your system, or `pdb_id = "1AKE"` |

The bundled Trp-cage system samples in well under a minute. A successful run ends with
`report.html` showing 8 frames, a non-zero pairwise RMSD spread, and an animating 3D viewer.

## Requirements

- **GPU:** NVIDIA, compute capability ≥ 7.0 (V100, T4, A100, RTX A6000, L40S, H100, B200).
  Below 7.0 is unsupported. Upstream's smallest tested device is an 80 GB A100, so VRAM
  headroom on smaller cards is unmeasured.
- **CUDA:** 12.6 via the `[cuda12]` extra. The `[cuda]` extra needs a CUDA 13 driver, which
  this estate does not have.
- **Network:** node 02 only (public ColabFold MSA server).
- **Disk:** ~1.91 GiB of model weights baked into the image.

## Troubleshooting

**Node 02 fails retrieving the MSA**
The public ColabFold server is unreachable or rate-limiting. Retry later, or run node 02
directly with Docker and mount pre-computed `<chain>.a3m` files into `/workspace/inputs` to skip
the server entirely.

**Node 03 reports out of GPU memory**
Lower `batch_size`. If it already fails at `batch_size = 1`, the system does not fit on one GPU;
context-parallel sampling needs four and is out of scope here.

**Node 03 cannot load the model**
The checkpoint is baked into the image under `$HF_HOME` and `HF_HUB_OFFLINE=1` is set, so this
is an image-build problem, not a network one. Rebuild and check the `snapshot_download` layer.

**The ensemble has collapsed (every frame identical)**
Node 05 warns when the maximum pairwise RMSD is near zero. Check that `mode` is `equilibrium`
rather than `structure`, and check the MSA depth in `msa_summary.json`.

**No frame passed the physicality checks**
Rejected frames are kept by design. Try a steering sampler (`edm_gd` or `edm_fk`), which is
exactly what they are for.

**`System.from_cif()` cannot import my structure**
It reads composition from mmCIF metadata (`_entity_poly` and friends), not from coordinates, so
a coordinates-only CIF will not import. Use an RCSB file, or write a system YAML.

## References

- [microsoft/bioemu2](https://github.com/microsoft/bioemu2) — BioEmu2, MIT licence.
- [microsoft/bioemu2 model card](https://huggingface.co/microsoft/bioemu2) — weights, pinned at
  revision `dd4b0c1b69984b68edcb6200a544213be4d75242`.
- [bioemu2 on PyPI](https://pypi.org/project/bioemu2/) — `0.1.0rc0`, released 2026-10-20.
- [bioemu2_pipeline.md](bioemu2_pipeline.md) — this workflow's design rationale, decisions and
  verification plan.
- BioEmu (v1): *Scalable emulation of protein equilibrium ensembles with generative deep
  learning*, Microsoft Research AI for Science.
