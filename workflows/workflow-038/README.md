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
[02-fetch-msa]  (network) ─► msa/, msa_summary.json, system.yaml, system_inventory.json
   │
   ▼
[03-sample]     (GPU) ─────► samples.xtc, topology.cif, physicality.json, msa/
   │                          sampling_manifest.json, system.yaml
   ├──────────────┐
   │              ▼
   │        [04-propagate] (GPU, identity when disabled)
   │              │   └──► propagated_samples.xtc, propagated_physicality.json,
   │              │        propagation_manifest.json
   │              │
   ├──────────────┼──────────► [05-analyze] ──► ensemble_metrics.csv,
   │              │                 │           ensemble_summary.json, ensemble_top.pdb
   │              │                 │
   └──────────────┴─────────────────┴──────────► [06-report] ──► report.html
```

Nodes run in the order 01 → 02 → 03 → 04 → 05 → 06, but the graph is not a single chain:
node 05 reads from both 03 and 04, and node 06 reads from 03, 04 and 05. That is why the
`enabled` flag in `propagation_manifest.json` has to be plumbed forward, since both later nodes
need to know whether node 04 actually propagated or passed the ensemble through.

Only nodes 03 and 04 need a GPU. Node 02 is the only node that needs network access; it writes
the alignment to disk so node 03 runs offline against the checkpoint baked into the image, with
`HF_HUB_OFFLINE=1` set. Note that `HF_HUB_OFFLINE` constrains Hugging Face only: if the
alignment were ever missing, BioEmu2 would fetch one itself from the public ColabFold server,
so nodes 03 and 04 now fail rather than sample without the alignment they were given.

## Input requirements

One system, supplied any of these ways (first match wins):

1. **A file** in the workflow's input directory — a BioEmu2 system YAML (`.yaml`/`.yml`, Boltz
   YAML auto-detected and converted), an mmCIF (`.cif`), a PDB (`.pdb`), or a FASTA
   (`.fasta`/`.fa`/`.faa`).
2. **`pdb_id`** — an RCSB accession, loaded as biological assembly 1 including all non-water
   ligands.
3. **`sequence`** — a raw one-letter protein sequence.

For a structure or an accession, `include_ligands` (default true) and `exclude_ligands` control
which non-polymer components come along. The default keeps every one of them, which is right —
dropping a component silently is worse than sampling one — but it does mean the cryoprotectant
is sampled too: 1CLL carries an ethanol beside its four structural calcium ions, and all five
become ligand tokens. The node lists every component it kept, by CCD code, so
`exclude_ligands=EOH` is a decision rather than a guess. The two parameters do not apply to a
system YAML, a FASTA or a bare sequence, which declare their own components; the node says so
when they are set anyway.

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
Boltz YAML is converted with `bioemu2 convert-boltz`. It then builds the all-atom structure,
which both counts the tokens and resolves every non-standard component, and **fails** if that
cannot be done: nodes 03, 04 and 05 build the same structure from the same system, so there is
nothing to be gained by discovering it later, on a GPU, after the checkpoint has loaded.

**Scientific notes:** The inventory records the **token** count, not just residues: BioEmu2's
memory and compute both scale with tokens (one per polymer residue, one per ligand atom), so
that is the number that predicts whether a run fits on the GPU. Over 10k tokens the node warns.
The count comes from `AllAtomStructure.num_tokens`, and `system_inventory.json` records which
attribute supplied it in `n_tokens_from`, so a bioemu2 API change shows up as changed provenance
rather than as a token total that is quietly null and a size warning that never fires. Measured:
1UBQ is 76 tokens for 76 residues, and 1CLL is 155 — 148 residues, 4 calcium ions and ethanol's
3 heavy atoms.

`n_polymer_residues` is the length of the **entity sequence**, which is what BioEmu2 models.
Upstream's `System.from_cif` documents it: "polymer sequences include unresolved residues from
the mmCIF entity tables". So the number is not the residue count of the input coordinates —
1CLL is a 148-residue protein whose crystal resolves 144, and the ensemble has all 148. The
inventory records that in `n_polymer_residues_from` and the report labels the row *entity
sequence*, because the two readings differ on most real entries and only one of them is what
was sampled.

Every non-polymer component is listed by CCD code in `ligand_components`, and a **placeholder**
component raises a warning: `UNX` (unknown atom or ion), `UNK` (unknown residue), `UNL`
(unknown ligand) and `DUM` (dummy atom) are what a depositor writes for density that could not
be identified, so sampling one means modelling a molecule whose chemistry nobody has
established. 5PTI carries a `UNX`. The warning names the codes and the `exclude_ligands` string
that drops them.

**Outputs:** `system.yaml`, `system_inventory.json`

### Node 02: Fetch MSA

**Goal:** The alignment on disk, so no later node needs the network.

**Process:** Calls `get_msa()` for the system, which retrieves an unpaired A3M per protein chain
from the public ColabFold server on a cache miss, optionally subsamples it to `msa_depth` rows,
and writes `msa/<chain>.a3m`. Pre-supplied A3M files in the node's input directory short-circuit
the whole step (Docker-only today — see open item 7).

A **heteromer gets two alignments per chain**, and both are written. ColabFold returns the
independent alignment and a positionally paired one whenever the system has more than one unique
protein sequence, and BioEmu2 uses both: its pairmixer is constructed with
`use_paired_feature=True`, and `prepare_structure()` substitutes zeros when the paired side is
absent. Writing only the unpaired half would therefore degrade every multi-sequence complex
without failing. The paired set goes to `msa/<chain>.paired.a3m`, and node 03 reattaches it.

Subsampling happens on the alignment object before the per-chain text is extracted, which
matters: `MSA.__getitem__` applies one set of row indices to every chain and to both alignments
at once, so the paired rows stay positionally matched. Subsampling the files individually
afterwards would silently break the pairing.

A monomer or homomer has one unique sequence and so has no paired alignment. That is correct,
not a missing file, and `msa_summary.json` records which case a run was in under `paired`.

**Scientific notes:** BioEmu2 is conditioned on the alignment, so its depth bounds what the
ensemble can be. Upstream names shallow or poor MSAs — orphan proteins, some viral proteins — as
a cause of degraded quality; it does not rank the causes, and neither does this workflow.
`msa_summary.json` records the depth per chain, and the node warns below 32 rows. **32 is this
workflow's own threshold**, not an upstream number: it is low enough that a real alignment
rarely trips it and high enough to catch a near-empty one. Treat it as a prompt to look at the
depth, not as a pass mark.

**Outputs:** `system.yaml`, `system_inventory.json`, `msa_summary.json`, `msa/` (one
`<chain>.a3m` per protein chain, plus `<chain>.paired.a3m` for a heteromer)

### Node 03: Sample Conformations

**Goal:** `num_samples` conformations, with a physicality verdict per frame.

**Process:** Builds an `EmuSession` from the system, the baked-in `bioemu2-preview` checkpoint
and node 02's alignment, samples, runs `check_physicality()`, and writes everything with
`save_frames()`. Uses the **Python API rather than the CLI** because the CLI has no flag for
supplying a pre-computed alignment. `--dry-run` resolves and prints every setting without
touching the GPU, and it is a real gate rather than a parameter echo: it checks that the
alignment was staged, resolves the system's topology — the step that builds every non-standard
component — and compares each staged A3M's query row against the sequence the system states for
that chain. All three run on a CPU in seconds, and each one is a defect that used to reach the
GPU before being noticed.

The node prints BioEmu2's own estimate of the largest batch that fits the GPU, computed from
the token count node 01 measured, and warns when the batch in use is above it. It does not
apply the estimate: `batch_size` changes which samples are drawn, so deriving it from whichever
card the job landed on would make a run irreproducible while still reporting the same seed.

It reattaches the paired alignment when node 02 wrote one. `load_a3m_msa` returns an unpaired
`MSA` whatever is on disk, so node 03 loads the per-chain files through `get_msa` and then
rebuilds the `MSA` with the paired texts added. A partial paired set is a failure rather than a
fallback: the paired alignment is matched across chains, so using some of it is not a degraded
version of using all of it.

It also checks that the alignment belongs to this system. bioemu2 does **not** compare an A3M's
query row against the sequence in the `System` — neither `get_msa` nor `load_a3m_msa` does — so
an alignment for a different protein is accepted, conditions the model on the wrong evolutionary
signal, and yields an ensemble whose frame count, physicality rate and geometry are all
unremarkable. Node 03 compares them itself, before the model loads, and refuses on a mismatch.
One case is reported rather than refused: where the system's sequence carries `X` for a
non-standard residue or a terminal cap, no alignment can carry a letter at that position, so a
legitimate file differs there by construction.

The node refuses to sample if node 02 reported protein chains and no `.a3m` arrived. That is
not defensiveness about a missing file: `EmuSession` fetches its own alignment when given
none, so the alternative is a node that quietly queries the ColabFold server, exits 0, and
records node 02's alignment in the manifest rather than the one it used. `sampling_manifest.json`
therefore carries `msa_loaded`, `msa_how` and `msa_files` describing what the model was given,
separately from `msa_source` describing what node 02 fetched.

**Scientific notes:** `mode` selects the task the model was trained for — `structure` sharpens
toward one accurate fold, `equilibrium` samples the ensemble. Physicality is **reported, never
enforced**: rejected frames stay in the XTC by design, as upstream intends, and the pass rate is
surfaced instead.

**Outputs:** `samples.xtc`, `topology.cif`, `physicality.json`, `system.yaml`,
`system_inventory.json`, `sampling_manifest.json`, `msa/`

### Node 04: Propagate Ensemble

**Goal:** One time-lagged endpoint per sampled conformation.

**Process:** When `propagate` is true, loads the ensemble against the system's own topology and
calls `EmuSession.propagate(lag_time_ns=…)` in equilibrium mode. When false, it copies the
ensemble through unchanged and records `enabled: false` — an identity propagation, so the
declared outputs always exist (a declared-but-unwritten output is a silent failure in Silva).
Nodes 05 and 06 read that flag and never present identity frames as time-lagged.

It reuses node 03's alignment rather than fetching one, and stops if node 03 recorded an
alignment that did not arrive here. Endpoints drawn against a freshly fetched alignment would
not match the ensemble they started from, and no later node could detect the difference, so
this is a failure and not a warning. `propagation_manifest.json` records which alignment was
used.

It runs `check_physicality()` on the endpoints and writes `propagated_physicality.json`, with
the same contract as node 03: every endpoint keeps its place and the verdict is reported, never
enforced. On a pass-through run that file is node 03's report copied forward, because the frames
are node 03's frames.

`sampler` is restricted here. `edm_fk` and `edm_fkc` turn on resampling, which copies particles
within a group of consecutive samples and so requires every frame in a group to be identical.
The starting frames are node 03's samples, which are all different, so the node rejects those
two up front with that reason rather than letting the failure surface from inside BioEmu2.
`edm` and `edm_gd` are supported; `edm_gd`'s preset sets `resampling_interval=None`.

**Scientific notes:** The result is one endpoint per start point at the chosen lag, **not** a
continuous path. Upstream recommends lags of 1–1000 ns, the range the model was trained for.

Measured once, on ubiquitin at the default 1 ns lag: 8 endpoints from 8 start points, all 8
physical, mean Rg 1.1628 nm against 1.1597 nm before propagation and a spread of 0.1943 nm
against 0.2204 nm. That is what a rigid fold should do over a nanosecond — move slightly and
stay folded — so it is evidence the path is wired correctly, and only weak evidence about the
dynamics. A system where the lag is expected to matter has not been propagated here, because
the ones in the test campaign that would show it do not fit the reference card.

**Outputs:** `propagated_samples.xtc`, `propagated_topology.cif`,
`propagated_physicality.json`, `propagation_manifest.json`

### Node 05: Analyze Ensemble

**Goal:** The ensemble as numbers.

**Process:** Binds the XTC to the system's own topology (no mmCIF parsing), then computes
per-frame radius of gyration and C-α RMSD to the ensemble medoid, a pairwise RMSD matrix,
per-residue RMSF about the mean structure, and a cluster count by greedy medoid assignment at
`rmsd_cluster_threshold_nm`. Writes a multi-model PDB of representative frames for the viewer,
superposed on the medoid so the viewer animates conformational change rather than tumbling.

The C-α set is selected as `protein and name CA`, not `name CA`. mdtraj matches atom names, and
a calcium ion is an atom named CA in a residue named CA, so the plain selection folds structural
calcium into the C-α set and shifts every RMSD, every RMSF and the medoid itself. Chains that
contribute no C-α atom are named in `ensemble_summary.json` under `excluded_chains`; they stay
in the frames and in the radius of gyration, which is computed over every atom.

**Scientific notes:** Every quantity is geometric. Cluster sizes are **frame counts**, not
populations, and this node deliberately computes no thermodynamic quantity. The pairwise pass
streams one row at a time, so the medoid and the ensemble spread are measured over every pair at
any frame count; above 200 frames only the matrix itself is dropped, which costs the heatmap and
the cluster count and nothing else.

The collapse warning fires below **0.01 nm** of maximum pairwise C-α RMSD. Not zero, and not
float epsilon: an XTC quantises coordinates to 0.001 nm, so one conformation written twice and
read back differs by about that much — measured 0.0018-0.0022 nm on a fixture of a single
conformation repeated with random rigid rotations — and the 1e-3 nm test this replaces could never
fire on a real collapsed ensemble, which is the only kind this node ever sees. 0.01 nm is an
order of magnitude above the file format's quantisation and an order below the smallest backbone
difference that means anything. It is a floor for "the frames are not identical", not a
statement that the ensemble has conformational spread; that threshold is roughly 0.1 nm, and it
is a reader's judgement, not a warning.

**Outputs:** `ensemble_metrics.csv`, `ensemble_summary.json`, `ensemble_top.pdb`

### Node 06: Report

**Goal:** One self-contained HTML report.

**Process:** Embeds all data inline as JavaScript and renders stat tiles, a per-frame Rg line
chart, an Rg histogram, a pairwise-RMSD heatmap on a single sequential hue, a per-residue RMSF
line chart, a per-frame table and a 3Dmol viewer that animates the representative frames.
Bootstrap 5, Plotly and 3Dmol come from CDN; no local assets.

The RMSF chart draws one trace per chain, dashed by chain and coloured by ensemble, against
position along the chain rather than author residue number. One trace across every chain strokes
back to the left edge wherever a chain restarts its numbering, which reads as a flexible region
and is an artefact of the axis. The viewer caption says how many of its frames are genuine
cluster centres and how many are filler, and which frame they were superposed on.

This node takes no `top_n`. How many frames are in the viewer is node 05's decision and is read
from `ensemble_summary.json`; a second copy of the parameter here could disagree with the file it
describes.

**Scientific notes:** The report **leads** with the interpretation caveats rather than burying
them. It states the alignment depth beside the results, and separately whether the alignment
reached the model at all, because those two were able to disagree without anything looking
wrong. Provenance records what the model was given, not only what node 02 fetched.

**Outputs:** `report.html`

## Parameters

### include_ligands / exclude_ligands

Which non-polymer components of a structure or RCSB entry come along. `include_ligands` defaults
to true, which is upstream's own default and the right one: every ion, cofactor and modified
residue stays, and dropping one silently would be worse than sampling one. The consequence is
that the crystallisation chemistry stays too, and BioEmu2 samples all of it as ligand tokens.
Measured on real entries: 1CLL carries an ethanol beside its four structural calcium ions
(155 tokens for a 148-residue protein), 5PTI carries a phosphate, a `UNX` placeholder and a
`DOD` — deuterated water, because it is a neutron structure, and upstream omits waters written
as HOH but not those written as DOD.

`exclude_ligands` takes CCD codes, comma separated, and drops those while keeping the rest:
`exclude_ligands=EOH,DOD`. It cannot be combined with `include_ligands=false`, which already
drops everything; the node refuses that combination rather than ignoring one of them. Node 01
lists every component it kept, by code, in the log and in `system_inventory.json`, and warns
about the two classes worth a second look — placeholders (`UNX`, `UNK`, `UNL`, `DUM`: density
nobody identified) and solvent (`DOD`). Neither parameter applies to a system YAML, a FASTA or a
bare sequence, which declare their own components.

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
paid for; peak memory scales with `batch_size`, not `num_samples` — on top of a floor set by
the token count, which for a large system is the whole of it.

**Test vs production:** 8 for a smoke test, 100+ for a real ensemble.

### batch_size

| Value / Range | Description |
|---------------|-------------|
| `0` (default) | Use the BioEmu2 default, which denoises all `num_samples` together. |
| 1–N | Denoise N at a time. Lower it when the GPU runs out of memory. |

**Trade-off:** Larger batches are faster up to a saturation point, then only use more memory.
BioEmu2 runs `ceil(num_samples / batch_size)` batches. Changing it changes the exact samples
drawn, and for `edm_fk`/`edm_fkc` also changes the resampling population.

**Test vs production:** `0` until you hit an out-of-memory error, then halve until it fits —
with the caveat that this works only while the batch is what the memory is going on. For a
system large relative to the card it is not: measured on 8 GB, 542 and 746 tokens failed at
batch 1 on the same allocation as at batch 8 (see **VRAM** under *Requirements*).

### num_steps

| Value / Range | Description |
|---------------|-------------|
| `0` (default) | Use the sampler preset's own step count, which is 200 for every sampler BioEmu2 ships. |
| 20–100 | Faster and coarser. Useful for a smoke test or for sweeping a parameter. |
| 200+ | Diminishing returns; the presets were tuned at 200. |

**Trade-off:** Runtime is close to linear in `num_steps`, and fewer steps means a less converged
denoising trajectory. Reducing it is a speed-for-quality trade, not a free win.

**Implementation note:** the step count belongs to the sampler object, not to
`EmuSession.sample()`. A default run therefore passes the sampler name as a string, which is
BioEmu2's documented path; only a non-zero `num_steps` builds the sampler with
`sampler_from_preset(preset, mode=…, num_steps=…)`, which applies the same step scale the string
path would. If that ever stops working the node fails and names the fallback, instead of
accepting the parameter and ignoring it.

**Test vs production:** `0`. Lower it only when you are deliberately trading quality for speed.

### sampler

| Value / Range | Description |
|---------------|-------------|
| `edm` (default) | Ordinary, unsteered EDM sampling. `EDMSampler(step_scale=1.0)` in equilibrium mode, 1.5 in structure mode. |
| `edm_gd` | Gradient guidance toward the default physicality potentials. `GuidedFKSampler(resampling_interval=None)`. |
| `edm_fk` | The same guidance plus Feynman-Kac resampling every 3 steps. Node 03 only. |
| `edm_fkc` | Feynman-Kac corrector, resampling when the effective sample size falls below 0.9. Node 03 only. |

**Trade-off:** The steered samplers raise the physicality pass rate at the cost of time, and
they tilt the sampled distribution — which is the point, and also a reason not to read the
result as an unbiased ensemble.

Node 04 accepts only `edm` and `edm_gd`. Resampling copies particles within a group of
consecutive samples, so BioEmu2 requires every starting frame in a group to be identical;
propagation starts from node 03's samples, which are all different. Node 04 rejects `edm_fk` and
`edm_fkc` with that reason instead of letting it fail inside the model.

**Test vs production:** `edm` unless the physicality pass rate is poor.

### physicality_filter

| Value / Range | Description |
|---------------|-------------|
| `high_quality` (default) | Upstream's default preset: small violation counts tolerated per frame. |
| `high_quality_per_residue` | The same checks with the counts divided by the residue count, so the tolerance scales with size. |
| `no_violations` | Every threshold at zero. |
| `no_filter` | Run no check at all; every frame is recorded as `not_checked`. |

The checks behind every preset, with upstream's own names, which are what appears in a frame's
`reason` field: `bad_bonds_backbone`, `bad_bonds_sidechain`, `bad_bonds_sulfur` (a bonded S–S
pair outside the disulfide window — no preset tolerates one), `within_residue_clashes`,
`valence_exceeded`, `isolated_atoms`, `same_chain_residue_clashes`,
`cross_chain_residue_clashes`, `duplicate_atoms`, `atoms_mismatch`, `bonds_mismatch`, and
`non_finite_atoms`, which rejects a frame under every preset. There is **no chirality check**,
and no backbone-dihedral or chain-break check either; this table previously claimed one.

**Trade-off:** None in output content — **frames are never removed**, whatever the setting. This
only changes the verdict recorded per frame.

With `no_filter` the manifest reports `n_passed` and `pass_rate` as null rather than as a clean
sweep. BioEmu2 returns an all-true `passed` array in that case, so counting it reports 100% of
frames passing a check that never ran.

### propagate / lag_time_ns

| Value / Range | Description |
|---------------|-------------|
| `propagate = false` (default) | Node 04 passes the ensemble through unchanged. |
| `propagate = true` | Draw one endpoint per frame, `lag_time_ns` ahead. |
| `lag_time_ns = 1.0` (default) | Upstream recommends 1–1000 ns, the trained range. |

**Trade-off:** Propagation roughly doubles GPU time. Outside 1–1000 ns the model is being asked
for something it was not trained on, and the node warns. On ubiquitin at 1 ns the ensemble
barely moves (see node 04's *Scientific notes*), which is the expected answer for a rigid fold
and not a reason to believe nothing happened.

### msa_depth / use_env

| Value / Range | Description |
|---------------|-------------|
| `msa_depth = 0` (default) | Keep the full alignment. |
| `msa_depth = N` | Subsample to N rows. |
| `use_env = true` (default) | Include the environmental database in the ColabFold search. |

**Trade-off:** Shallower alignments give the MSA module less to process, so the trunk embedding
is faster and its activations smaller. The pair representation is quadratic in tokens and is
unaffected by depth, so depth does not buy back the cost of a large system. (An earlier version
of this line asserted a combined `rows × tokens²` cost; that product is not something upstream
states and is not something this workflow has measured.) Subsampling degrades ensemble quality
and is a deliberate experiment, not an optimisation.

`msa_depth` applies only to an alignment fetched from the server. A pre-supplied `<chain>.a3m`
is used exactly as given, and node 02 warns rather than appearing to subsample it.

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

A JSON list, one entry per frame: zero-based `frame`, `status` (`passed` / `rejected` /
`not_checked`) and the first failed check in `reason`. **Rejected frames remain in the XTC** —
that is upstream's design, not a bug. To keep only clean frames, read the list and slice the
trajectory on the frames whose status is `passed`.

Only multi-residue chains containing at least one standard amino acid are checked. A
nucleic-acid chain, a ligand or a lone ion is not, and `not_checked` for the whole file means
either `physicality_filter = no_filter` or a system with no such chain —
`sampling_manifest.json` says which in `physicality_not_checked_because`.

`propagated_physicality.json` is the same format for node 04's endpoints. It carries one extra
status, `unavailable`, used when the check itself could not run: that is deliberately not one of
BioEmu2's three values, so it cannot be read as "not checked".

### ensemble_metrics.csv

One row per frame per ensemble: `radius_of_gyration_nm` (compactness — a flat series across an
equilibrium run means the ensemble collapsed), `rmsd_to_medoid_nm` (how unusual the frame is
relative to the ensemble centre), and the physicality verdict.

### ensemble_summary.json

Ensemble-level: Rg mean/std/range, pairwise RMSD statistics and the full matrix (up to 200
frames), per-residue RMSF, the medoid frame, and clusters at the chosen threshold. **Cluster
sizes are frame counts.** They are not populations, not probabilities, and not free energies,
because BioEmu2's sampling frequency is not a calibrated Boltzmann probability.

`system` carries node 01's findings forward rather than re-deriving them: the source, the chain
ids, the token count, `n_polymer_residues` with its `n_polymer_residues_from` provenance (the
entity sequence, which includes residues the input structure does not resolve), and
`ligand_components`, every non-polymer component by CCD code and asym id. The report shows that
list, because "has ligands: yes" does not distinguish four structural calcium ions from one
ethanol, and BioEmu2 samples both.

`rmsd_atom_selection` records the selection the RMSD and RMSF were computed over, and
`excluded_chains` the chains that contributed no C-α atom. `viewer_frames` are the frames in
`ensemble_top.pdb`, `viewer_frames_from_clusters` says how many of them are genuine cluster
centres rather than consecutive filler, and `viewer_frames_superposed_on_frame` names the frame
they were superposed on. Each per-residue RMSF entry carries both `chain` (the chain id) and
`chain_index`, so a multi-chain trace can be split by chain instead of being drawn as one line
across a discontinuity.

### report.html

Self-contained. Opens in any browser with no local assets.

## Quick start

### Running with Docker

Built from [Dockerfile](Dockerfile) as `bioemu2:2026_11_20`.

```bash
docker build -t bioemu2:2026_11_20 workflows/workflow-038/
```

Smoke-test the image before running any node:

```bash
docker run --rm bioemu2:2026_11_20 bioemu2 --help
```

```bash
docker run --rm --gpus all bioemu2:2026_11_20 nvidia-smi --query-gpu=name,compute_cap --format=csv
```

Compute capability must be 7.0 or newer. Then run a node. Node 03 reads its inputs from
`./inputs`, so the directory has to be populated first — running it against an empty one fails
on the missing `system.yaml`, which is the node working correctly and not a build problem. The
shortest path is to run nodes 01 and 02 first and hand their outputs on:

```bash
mkdir -p /tmp/wf038/{inputs,outputs}
```

```bash
docker run --rm -v /tmp/wf038/outputs:/workspace/outputs -v $(pwd)/workflows/workflow-038/01-prepare-system:/workspace -w /workspace -e PARAM_PDB_ID=1UBQ bioemu2:2026_11_20 bash run.sh
```

Then copy that node's `outputs/` into the next node's `inputs/` and repeat for node 02, whose
outputs are node 03's inputs:

```bash
docker run --rm --gpus all -v /tmp/wf038/inputs:/workspace/inputs -v /tmp/wf038/outputs:/workspace/outputs -v $(pwd)/workflows/workflow-038/03-sample:/workspace -w /workspace -e PARAM_MODE=equilibrium -e PARAM_NUM_SAMPLES=8 -e PARAM_SEED=42 bioemu2:2026_11_20 bash run.sh --dry-run
```

Drop `--dry-run` to actually sample. Both GPU nodes support it, and node 03's dry run checks
the alignment was staged, resolves the system's topology and compares each A3M's query row
against the sequence the system states — a real test of the plumbing, the chemistry and the
input's identity rather than of the parameters alone. It needs no GPU.

The whole suite in [tests/](tests/README.md) is built on that: `tests/run_fixtures.sh` runs
nodes 01, 02 and the node 03 dry run over seven RCSB entries plus nodes 04-06 over
deterministic ensembles, all offline, and `tests/run_gpu_campaign.sh` runs the real thing on a
card.

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

Measured on an RTX 3060 Ti (8 GB, compute capability 8.6, driver 580.178.04): ubiquitin (1UBQ,
76 tokens), 8 samples in `equilibrium` mode, 31.2 s in node 03, 8 of 8 frames passing
physicality. The bundled Trp-cage system is smaller again. A successful run ends with
`report.html` showing 8 frames and an animating 3D viewer.

Do not read "non-zero pairwise RMSD" as a pass: two independently sampled frames always differ,
so that condition is satisfied by any output at all, including a collapsed one. For a folded
protein in `equilibrium` mode the maximum pairwise C-α RMSD should clear roughly 0.1 nm, the
scale of ordinary backbone fluctuation; below that the ensemble is one structure plus noise.
That figure is a floor for spotting collapse, not a validation of the ensemble. Node 05 warns
below 0.01 nm, which is a different and much weaker claim — it is the point at which the frames
stop being distinguishable from one conformation written twice, given that an XTC quantises
coordinates to 0.001 nm. Clearing 0.01 nm means the frames differ; clearing 0.1 nm means the
ensemble has conformational spread. Measured: 1UBQ, 8 samples, max pairwise C-α RMSD 0.2204 nm,
Rg 1.1597 ± 0.0056 nm against an experimental 1.17 nm.

## Requirements

- **GPU:** NVIDIA, compute capability ≥ 7.0 (V100, T4, A100, RTX A6000, L40S, H100).
  Not B200: sm_100 needs CUDA 12.8 or newer and this image is built on cu126. **Measured** in
  the built image — `torch.cuda.get_arch_list()` returns
  `['sm_50', 'sm_60', 'sm_70', 'sm_75', 'sm_80', 'sm_86', 'sm_90']` — so `sm_70` is present and
  the 7.0 floor is real rather than inherited, and `sm_100` is absent, which is the B200
  exclusion.
- **VRAM:** the limit is the **token count**, not the sample count. Measured on an 8 GB RTX
  3060 Ti (sm_86, 7.66 GiB usable): 1UBQ (76 tokens, 31.2 s), 5PTI (65), 2ZTA (72) and 1BNA
  (69) sample 8 in one batch in well under a minute, and 1CLL (155) samples. 1AKE (542 tokens)
  and 1A3N (746) **run out of memory at batch 8 and again at batch 1** — each time on the same
  allocation, 1.06 GiB and 1.46 GiB respectively — so on this card lowering the batch does not
  move the ceiling at all. Upstream's own estimate, which node 03 prints, said 33 and 19 would
  fit; it is linear in the token count and calibrated on an 80 GB A100, and on a small card it
  is an upper bound to come down from rather than a ceiling. Plan by size: 8 GB holds a few
  hundred tokens, somewhere above 155 and below 542, and a system past that needs a larger
  card regardless of `batch_size`.
- **CUDA:** 12.6 via the `[cuda12]` extra, which pins `torch==2.13.0+cu126`. The `[cuda]` extra
  pins no torch at all — it would resolve anywhere in BioEmu2's `>=2.11,<2.14` range, so the
  image would differ between builds — and it needs the CUDA 13 cuEquivariance ops, which raise
  the driver floor. The cu13 kernels are a speed optimisation, not a capability.
- **Network:** node 02 for the public ColabFold MSA server, and nodes 01, 03, 04 and 05 **only**
  for a non-standard component whose CCD definition is not already in the image. Anything that
  is not a standard residue — an ion, a cofactor, a terminal cap, a modified residue — is built
  from its CCD entry, which bioemu2 fetches from `files.rcsb.org` the first time it meets the
  code, from inside `system.topology`. The image pre-warms 126 codes (ions, caps, cofactors,
  sugars, crystallisation additives, the frequent modified residues) into `BIOEMU2_CCD_CACHE`
  and `BIOEMU2_TEMPLATE_CACHE`, so the common cases need no network at all; an unlisted code
  does, and every node that builds a topology says exactly that, naming the component.
- **Disk:** 9.58 GB image, of which ~1.91 GiB is model weights and ~2 MB the CCD cache.

## Troubleshooting

**Node 02 fails retrieving the MSA**
The public ColabFold server is unreachable or rate-limiting. Retry later, or run node 02
directly with Docker and mount pre-computed `<chain>.a3m` files into `/workspace/inputs` to skip
the server entirely. bioemu2 also caches each chain's result under `BIOEMU2_MSA_CACHE`, which
the image points at `/models/colabfold_msa` inside the container: mount a host directory there
and a repeated run skips the server. `tests/run_gpu_campaign.sh` does exactly that, because the
MSA fetches rather than the GPU are most of its wall time.

**Node 03 reports out of GPU memory**
Try a lower `batch_size`, but check the token count in the same log before spending long on it.
Measured on an 8 GB card, a 542-token and a 746-token system each failed at `batch_size = 1` on
the same allocation that failed at batch 8: the peak was set by the size of the system, not by
how many samples were denoised together. If batch 1 fails, the system needs a larger card —
context-parallel sampling needs four GPUs and is out of scope here. Do not read the batch
estimate node 03 prints as a ceiling either: it is upstream's formula, and it said 33 would fit
for the system that fitted none.

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

**A node fails with a connection error for `files.rcsb.org/ligands/download/XXX.cif`**
`XXX` is a chemical component the image does not already carry. bioemu2 builds every
non-standard component from its CCD definition and downloads it the first time it meets the
code, from inside `system.topology` — which nodes 01, 03, 04 and 05 all call. Either give that
node network access, or add the code to the pre-warm list in the `Dockerfile` and rebuild, or
drop the component with `exclude_ligands=XXX`. The node names the code in the error for exactly
this reason.

**A component cannot be modelled at all**
`UNL`, the unknown-ligand placeholder, raises inside bioemu2's own template builder even with
the CCD in hand. Node 01 fails, lists the system's non-polymer components and names the
`exclude_ligands` parameter. `UNX`, `UNK` and `DUM` do build, and node 01 warns about them
instead: they are a depositor's codes for density that could not be identified, so sampling one
means modelling a molecule of unknown chemistry.

**The log is full of `libcue_ops.so: libcuda.so.1`, pynvml failures and "NVIDIA RTX A6000"**
Expected on the CPU nodes (01, 02, 05, 06), which run without `--gpus`. bioemu2 probes for the
CUDA equivariance ops and for NVML, finds neither, and falls back to a built-in default device
description — the "RTX A6000" in that line is a placeholder, not the card in the machine. It is
noise. On a GPU node the probe succeeds and the lines do not appear; `nvidia-smi` in node 03's
log is the authoritative statement of which device ran.

## References

- [microsoft/bioemu2](https://github.com/microsoft/bioemu2) — BioEmu2, MIT licence.
- [microsoft/bioemu2 model card](https://huggingface.co/microsoft/bioemu2) — weights, pinned at
  revision `dd4b0c1b69984b68edcb6200a544213be4d75242`.
- [bioemu2 on PyPI](https://pypi.org/project/bioemu2/) — `0.1.0rc0`, released 2026-10-20.
- [bioemu2_pipeline.md](bioemu2_pipeline.md) — this workflow's design rationale, decisions and
  verification plan.
- BioEmu (v1): *Scalable emulation of protein equilibrium ensembles with generative deep
  learning*, Microsoft Research AI for Science.
