---
doc_id: workflow-036
domain: materials-science
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  Estimates the glass-transition temperature (Tg), thermal expansion, and
  high-temperature dimensional stability of heat-resistant battery-case and
  separator resins (PPS, PA66, PBT, PEEK, PP) using all-atom GAFF2 MD
  (GROMACS) and a stepwise melt-quench cooling series, delivering an HTML
  report with a density/specific-volume Tg kink fit and a Mol* structure view.
tags: [molecular-dynamics, polymer, gromacs, materials, glass-transition, thermal-expansion, gaff2, tg]
---

# Workflow 036: Polymer MD — Heat-Resistant Plastics for EV Battery Protection

All-atom molecular-dynamics pipeline that builds an amorphous polymer cell for
one of five candidate resins, runs a melt-quench cooling series in GROMACS
with the GAFF2 force field, fits a bilinear (two-segment) kink to the
density/specific-volume vs. temperature curve to estimate Tg, and reports
thermal expansion (CTE) and high-temperature dimensional stability (volume
change) — the properties that determine whether a resin protects an EV
battery case or separator from heat.

---

## When to use this workflow

Use this workflow to get a fast, first-principles ranking of candidate
heat-resistant resins — PPS, PA66, PBT, PEEK, or PP (reference) — by Tg,
thermal expansion, and dimensional stability at a chosen service temperature,
without running an experimental DMA/DSC test programme.

The default run times are short (150 ps melt, 100 ps per quench stage) for
rapid screening; increase `stage_time_ps` to 500–1000 ps for production-
quality density convergence and a more reliable Tg. Do not use this workflow
for true semi-crystalline morphology prediction — crystallinity is
approximated via initial packing density only (see Node 01 below), not a
real crystal lattice.

---

## Workflow structure

```
01-build-cell → 02-cooling-series → 03-measure-tg → 04-report
```

| Node | Image | Input | Output |
|------|-------|-------|--------|
| `01-build-cell` | `polymer_md:2026_06_26` | `global_params.json` | `system.gro`, `topol.top`, `cell.pdb`, `build_report.json` |
| `02-cooling-series` | `gromacs:2025_11_12` | `system.gro`, `topol.top`, `build_report.json` | `cooling_series.json`, `equilibrated.gro` |
| `03-measure-tg` | `polymer_md:2026_06_26` | `cooling_series.json`, `build_report.json` | `tg_report.json` |
| `04-report` | `polymer_md:2026_06_26` | `tg_report.json`, `build_report.json`, `cell.pdb` | `report.html`, `summary.json` |

No new Docker images were built for this workflow — it reuses
`polymer_md:2026_06_26` and `gromacs:2025_11_12`, already published from
workflow-030/031, since the toolchain (RDKit, AmberTools, packmol, acpype,
GROMACS) is identical.

---

## Node details

### Node 01 — Build Cell

Builds a GAFF2-parametrised amorphous polymer cell using AmberTools and GROMACS:

1. **RDKit** — generates 3-D coordinates for a short end-capped oligomer (ETKDG + MMFF), writing both a PDB (coordinates, for packmol) and an SDF (exact Kekulized bond orders, for antechamber)
2. **antechamber** — computes BCC partial charges from the SDF
3. **parmchk2** — generates missing GAFF2 torsion/vdW parameters (`.frcmod`)
4. **packmol** — packs `n_chains` oligomer chains into a cubic box at a crystallinity-scaled fraction of target density
5. **tleap** — builds AMBER topology (`.prmtop`, `.inpcrd`)
6. **acpype** — converts AMBER topology to GROMACS format (`.top`, `.gro`)

Feeding antechamber the SDF instead of a bond-order-less PDB avoids a known
failure mode: for aromatic/conjugated backbones (PPS, PBT, PEEK), antechamber
can mis-guess hybridisation from geometry alone, which `tleap` then rejects
with `... has force field coordination 4 but only 3 bonded neighbors`. This is
the same class of bug workflow-031 hit and fixed for PET.

**Crystallinity** (`low` / `medium` / `high`) is approximated by varying the
initial packing density fraction (0.55 / 0.65 / 0.75 of the resin's target
amorphous density) — a **trend-only proxy**, not a real semi-crystalline
lattice. Building a true crystalline unit cell per resin was out of scope;
`build_report.json` documents this explicitly.

### Node 02 — Cooling Series

Chained GROMACS stages (each continuing from the previous stage's coordinates
and velocities), using 8 OpenMP threads + GPU offload (GAFF2 / PME):

| Stage | Ensemble | Duration |
|-------|----------|----------|
| Energy minimisation | — | `em_steps` (steepest descent) |
| Melt | NPT at resin `melt_temp_c` | `melt_time_ps` |
| Quench 1 | NPT at 200 °C | `stage_time_ps` |
| Quench 2 | NPT at 150 °C | `stage_time_ps` |
| Quench 3 | NPT at 80 °C | `stage_time_ps` |
| Quench 4 | NPT at 25 °C | `stage_time_ps` |

This fixed 5-point ladder (melt + 200/150/80/25 °C) always runs, independent
of the student's chosen `temperature` parameter — it's what makes the Tg
kink-fit possible, and the four quench points already coincide with the four
selectable `temperature` levels, so one MD run answers both "what is Tg?" and
"what are this resin's properties at my chosen temperature?" Density is
averaged over the last third of each stage's trajectory into
`cooling_series.json`.

**Barostat compressibility**: `melt.mdp`/`quench.mdp` use `compressibility =
2.0e-6 bar⁻¹` (polymer-melt order of magnitude), not the `4.5e-5 bar⁻¹`
(water's compressibility) inherited from workflow-030/031's mdp files. With
the water-like value, PPS's cell drifted/expanded instead of condensing —
density bounced noisily around 520–610 kg/m³ (well under half of the 1350
kg/m³ target) across every stage, even at 10x longer stage times, and the
resulting Tg fit was degenerate. Lowering compressibility to a value
appropriate for a stiff polymer restored a monotonic density-vs-temperature
trend and a numerically sane (`tg_reliable: true`) fit on the very same
short default run length — see
[issue #196](https://github.com/chiral-data/collab-workflows/issues/196) for
the diagnosis. workflow-030/031 shared this mdp template and had the same
latent issue — confirmed and fixed there too, see
[issue #198](https://github.com/chiral-data/collab-workflows/issues/198).

### Node 03 — Measure Tg

Pure-Python analysis (no numpy):

- **Tg**: fits a two-segment (bilinear) line to specific volume vs.
  temperature. With 5 ladder points there are exactly two candidate
  breakpoints that leave ≥2 points on each side — the minimum viable case for
  a kink fit, chosen by minimising total least-squares error; Tg is the
  intersection of the two fitted lines.
- **Thermal expansion**: CTE (glassy and rubbery) from each segment's slope
  of specific volume vs. temperature.
- **Dimensional stability**: % volume change between the 25 °C reference
  point and the student's selected `temperature`, read directly off the same
  series (no extra simulation).
- Literature Tg is carried through from Node 01 as a labeled reference — not
  computed by this pipeline.

### Node 04 — Report

Self-contained `report.html`:

- KPI strip: estimated Tg vs. literature, CTE (glassy/rubbery), volume change % at the selected temperature
- Density/specific-volume vs. temperature chart with the two fitted lines and the Tg kink marked
- Cross-resin Tg bar chart vs. literature values
- Dimensional-stability table (density and specific volume at every ladder stage)
- **Mol\*** 3-D viewer (`molstar@latest` via CDN) of the representative packed cell (`cell.pdb`)
- Callouts documenting the two known simplifications (GAFF2 vs. ideal OPLS-AA; crystallinity as a packing-density proxy)

---

## Parameters (`global_params.json`)

| Parameter | Default | Description |
|-----------|---------|--------------|
| `resin_type` | `"PPS"` | Resin identifier: `PPS`, `PA66`, `PBT`, `PEEK`, or `PP` (reference) |
| `temperature` | `150.0` | Measurement temperature (°C): 25, 80, 150, or 200 — selects which ladder point is reported for dimensional stability |
| `crystallinity` | `"medium"` | `low`, `medium`, or `high` — approximated via initial packing density (see Node 01) |
| `n_chains` | `20` | Number of oligomer chains in the amorphous cell |
| `em_steps` | `5000` | Steepest-descent EM steps |
| `melt_time_ps` | `150.0` | NPT melt-stage duration (ps) |
| `stage_time_ps` | `100.0` | Duration (ps) of each of the four NPT quench stages |

---

## Resin library

| Resin | Oligomer proxy | Target density (g/cc) | Melt stage temp | Literature Tg |
|---|---|---|---|---|
| PPS | 5-ring phenylene sulfide | 1.35 | 285 °C | 88 °C |
| PA66 | 3-mer hexamethylenediamine/adipic acid | 1.14 | 260 °C | 50 °C |
| PBT | 3-mer butanediol/terephthalate | 1.31 | 225 °C | 45 °C |
| PEEK | 3-mer ether-ether-ketone | 1.30 | 343 °C | 143 °C |
| PP (ref.) | 5-mer propylene | 0.855 | 200 °C | −10 °C |

PA66's repeat unit is heavier than the other resins' proxies, so it uses a
3-mer (like PBT/PEEK) rather than a 5-mer — a 5-mer PA66 oligomer made
antechamber's AM1-BCC charge SCF (`sqm`) run for 15+ minutes, well past what's
reasonable for the "rapid screening" default.

---

## Known simplifications

- **Force field**: aromatic backbones (PPS, PEEK, PBT) are chemically better
  suited to OPLS-AA (per the original brief), but this workflow reuses the
  GAFF2/antechamber/acpype pipeline shared with workflow-030/031, since it
  auto-parametrises arbitrary SMILES-derived oligomers with no manual
  atom-typing — there is no existing template in this repo for hand-authoring
  OPLS-AA residues per resin. Tracked in
  [issue #196](https://github.com/chiral-data/collab-workflows/issues/196).
- **Crystallinity** is a packing-density proxy at cell-build time, not a real
  semi-crystalline lattice — trend comparison only.
- **Tg fit** uses the minimum viable 5-point ladder for a bilinear kink fit;
  increase `stage_time_ps` for a better-converged curve.

---

## Real pipeline

Runs the actual RDKit/AmberTools/GROMACS pipeline by default: each node's
`run.sh` (what `job.toml` executes) delegates to `run_real.sh`. For a fast,
deterministic mock run of this same pipeline (downloads pre-computed outputs
instead of computing them, no Docker/GPU/wait required), see
[workflow-032](../workflow-032/README.md).

---

## Docker images

| Image | Used by | Contents |
|-------|---------|----------|
| `polymer_md:2026_06_26` | 01, 03, 04 | AmberTools 22 (antechamber, parmchk2, tleap), acpype, RDKit, Python 3 |
| `gromacs:2025_11_12` | 02 | GROMACS 2023.2 (AVX2_256, GPU-offload), 8 OpenMP threads |
