---
doc_id: workflow-030
domain: materials-science
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  Predicts Young's modulus and specific modulus of a neat or glass-fibre-filled
  polypropylene resin using all-atom GAFF2 MD (GROMACS) and the Halpin–Tsai GF
  correction model, delivering an HTML report with mechanical and density results.
tags: [molecular-dynamics, polymer, gromacs, materials, specific-modulus, gaff2, halpin-tsai]
---

# Workflow 030: Polymer MD — Specific Modulus Prediction (PP)

All-atom molecular-dynamics pipeline that builds a polypropylene (PP) simulation cell,
runs melt-quench equilibration and NPT production in GROMACS with the GAFF2 force field,
applies a Halpin–Tsai glass-fibre (GF) correction for filled grades, and reports Young's
modulus, density, and specific modulus (E / ρ) — the key metric for EV lightweighting.

---

## When to use this workflow

Use this workflow when you need a fast, first-principles estimate of the specific modulus
of neat PP or PP + GF compounds at a given temperature and fibre loading, without running
a full experimental test programme.  The default run times are intentionally short
(500 ps equil, 1 000 ps production) for rapid screening; increase `equil_time_ps` /
`prod_time_ps` to 5 000–10 000 ps for production-quality density convergence.

Do not use this workflow for non-PP resins (the current force-field parametrisation
targets PP only), or for crystalline or semi-crystalline morphology prediction (the
pipeline builds an amorphous cell from scratch).

---

## Workflow structure

```
01-build-cell → 02-equilibrate → 03-measure-properties → 04-apply-gf-correction → 05-report
```

| Node | Image | Input | Output | Wall time (2026-06-28, GPU) |
|------|-------|-------|--------|-----------------------------|
| `01-build-cell` | `polymer_md:2026_06_26` | `global_params.json` | `system.gro`, `topol.top`, `build_report.json` | ~30 s |
| `02-equilibrate` | `gromacs:2025_11_12` | `system.gro`, `topol.top` | `equilibrated.gro`, `density.xvg`, `equil_report.json` | 75.6 s |
| `03-measure-properties` | `gromacs:2025_11_12` | `equilibrated.gro`, `topol.top` | `properties.json`, `density_prod.xvg`, `stress_strain.xvg` | 46.2 s |
| `04-apply-gf-correction` | `polymer_md:2026_06_26` | `properties.json`, `build_report.json` | `corrected_properties.json` | ~2 s |
| `05-report` | `polymer_md:2026_06_26` | `corrected_properties.json`, `build_report.json` | `report.html`, `summary.json` | ~2 s |

**Total wall time (default params, GPU):** ~2.5–3 min

---

## Node details

### Node 01 — Build Cell

Builds a GAFF2-parametrised amorphous PP simulation cell using AmberTools and GROMACS.

1. **antechamber** — computes BCC partial charges on a single PP oligomer (50 atoms, 5-mer)
2. **parmchk2** — generates missing GAFF2 torsion/vdW parameters (`.frcmod`)
3. **packmol** — packs 20 oligomer chains into a cubic box at 60 % of target density (~0.27 s)
4. **tleap** — builds AMBER topology (`.prmtop`, `.inpcrd`)
5. **acpype** — converts AMBER topology to GROMACS format (`.top`, `.gro`)

Output `build_report.json` records box size, MW per chain, and fibre loading for downstream nodes.

### Node 02 — Equilibrate

Three-stage GROMACS NPT equilibration using 8 OpenMP threads + GPU offload (GAFF2 / PME):

| Stage | Ensemble | Duration | Wall time |
|-------|----------|----------|-----------|
| Energy minimisation | — | 5 000 steps (steepest descent) | 8.9 s |
| Melt NPT | NPT at 473 K | 200 ps (200 000 steps, 1 fs dt) | 22.0 s |
| Equilibrate NPT | NPT at target T | 500 ps (500 000 steps, 1 fs dt) | 43.7 s |

Density is averaged over the last third of the equilibration trajectory;
the result is stored in `equil_report.json`.

### Node 03 — Measure Properties

Two GROMACS runs to extract Young's modulus and density:

| Stage | Ensemble | Duration | Wall time |
|-------|----------|----------|-----------|
| Production NPT | NPT at target T | 1 000 ps (1 000 000 steps, 1 fs dt) | 43.7 s |
| Deformation NVT | NVT at target T | 50 ps (50 000 steps, 1 fs dt, 2 % strain) | 2.5 s |

Young's modulus is derived from the stress–strain curve (`stress_strain.xvg`).
Specific modulus is E / ρ in kN·m/kg.

### Node 04 — Apply GF Correction

Pure-Python Halpin–Tsai composite model.  When `fiber_loading > 0`, it blends the
neat-resin modulus (from Node 03) with an assumed glass-fibre modulus (72 GPa, aspect
ratio 20) to give `corrected_youngs_modulus_gpa` and `corrected_specific_modulus_kNm_kg`.
At 0 wt% loading the node is a pass-through (no correction applied).

### Node 05 — Report

Generates `report.html` (self-contained HTML table + EV lightweighting context) and
`summary.json` from the corrected properties.

---

## Parameters (`global_params.json`)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `resin_type` | `"PP"` | Resin identifier (PP only in v1) |
| `temperature` | `23.0` | Measurement temperature (°C) |
| `fiber_loading` | `0.0` | Glass-fibre loading (wt%) |
| `n_chains` | `20` | Number of oligomer chains in the cell |
| `em_steps` | `5000` | Steepest-descent EM steps |
| `melt_time_ps` | `200.0` | Melt NPT duration (ps) |
| `equil_time_ps` | `500.0` | Equilibration NPT duration (ps) |
| `prod_time_ps` | `1000.0` | Production NPT duration (ps) |
| `deform_time_ps` | `50.0` | Deformation NVT duration (ps) |

---

## Sample output (2026-06-28, default params)

```json
{
  "resin_type": "PP",
  "temperature_c": 23.0,
  "fiber_loading_wt_pct": 0.0,
  "youngs_modulus_gpa": -0.026,
  "density_kg_m3": 383.2,
  "specific_modulus_kNm_kg": -67.8,
  "gf_correction_applied": false,
  "simulation_note": "fiber_loading=0 — neat resin properties"
}
```

> **Note on negative modulus:** The default run times are intentionally short for speed.
> A negative Young's modulus indicates the cell has not fully converged — increase
> `equil_time_ps` to ≥ 5 000 ps and `prod_time_ps` to ≥ 5 000 ps for physically
> meaningful results.  The `sample_outputs/` folder contains a converged reference run.

---

## Docker images

| Image | Used by | Contents |
|-------|---------|----------|
| `polymer_md:2026_06_26` | 01, 04, 05 | AmberTools 22 (antechamber, parmchk2, tleap), acpype, Python 3 |
| `gromacs:2025_11_12` | 02, 03 | GROMACS 2023.2 (AVX2_256, GPU-offload), 8 OpenMP threads |
