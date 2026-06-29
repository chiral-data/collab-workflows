---
doc_id: workflow-031
domain: materials-science
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  Measures the O₂ or H₂O diffusion coefficient through a polymer film (LDPE, PP,
  EVOH, PA6, or PET) using all-atom GAFF2 MD (GROMACS) and the Einstein relation,
  delivering an HTML report with MSD curve and barrier ranking vs. literature values.
tags: [molecular-dynamics, polymer, gromacs, diffusion, barrier-film, gaff2, spce, trappe]
---

# Workflow 031: Barrier Films — Plastics That Protect Food

All-atom molecular-dynamics pipeline that builds an amorphous polymer cell,
equilibrates it via melt-quench NPT, inserts small-molecule penetrants (O₂ or H₂O),
runs an NVT diffusion trajectory, and reports the diffusion coefficient D via the
Einstein relation — the key quantity for food-packaging barrier performance.

---

## When to use this workflow

Use this workflow when you need a first-principles estimate of gas or moisture barrier
performance for a candidate resin, or to rank five common packaging polymers (LDPE, PP,
EVOH, PA6, PET) by their diffusivity without running permeability lab tests.

The default run times are short for rapid screening. Increase `diff_time_ps` to
≥ 20 000 ps for high-barrier resins (EVOH, PET) where slow diffusion requires long
simulations to reach the true Einstein (log-log slope ≈ 1) regime.

---

## Parameters

| Parameter | Default | Description |
|---|---|---|
| `resin_type` | `LDPE` | Polymer matrix: `LDPE`, `PP`, `EVOH`, `PA6`, or `PET` (see note below) |
| `penetrant` | `O2` | Small molecule: `O2` (TraPPE 2-site) or `H2O` (SPC/E) |
| `temperature` | `23.0` | Test temperature in °C |
| `n_chains` | `20` | Oligomer chains in the simulation cell |
| `em_steps` | `5000` | Energy minimisation steps |
| `melt_time_ps` | `200.0` | NPT melt duration (ps) |
| `equil_time_ps` | `500.0` | NPT equilibration at target temperature (ps) |
| `n_penetrant` | `5` | Penetrant molecules inserted |
| `diff_time_ps` | `5000.0` | NVT diffusion run length (ps) |

---

## Resin support and known limitations

| Resin | Default SMILES proxy | Status |
|---|---|---|
| `LDPE` | n-decane (C₁₀) | ✅ Supported (default) |
| `PP` | propylene 5-mer | ✅ Supported |
| `EVOH` | E/VOH 5-mer | ✅ Supported |
| `PA6` | caprolactam 5-mer | ✅ Supported |
| `PET` | terephthalate 5-mer | ⚠️ **Known issue — see below** |

### PET: tleap connectivity failure

PET contains aromatic rings. When packmol assembles the multi-chain cell as a PDB file,
bond-order information is lost. tleap then misidentifies the hybridisation of some ring
carbons and exits with:

```
Atom .R<UNL 1>.A<C36 47> has force field coordination 4 but only 3 bonded neighbors.
```

**Workaround (not yet automated):** switch packmol to mol2 I/O so bond types are
preserved through the full pipeline. Until that is fixed, use `LDPE`, `PP`, `EVOH`, or
`PA6` for automated runs. PET manual workaround tracked in the PR.

---

## Node pipeline

```
01-build-cell  →  02-equilibrate  →  03-insert-penetrant  →  04-diffusion-md  →  05-report
(barrier_films)    (gromacs)          (gromacs)               (gromacs)           (barrier_films)
```

### 01-build-cell
RDKit 3-D conformer → antechamber GAFF2 + AM1-BCC charges → parmchk2 → packmol
(mol2 I/O) → tleap → acpype → `system.gro`, `topol.top`, `build_report.json`.

### 02-equilibrate
Energy minimisation → NPT melt at `melt_temp_c` (Berendsen barostat) →
NPT quench to `temperature` → `equilibrated.gro`, `density.xvg`, `equil_report.json`.

### 03-insert-penetrant
Writes force-field ITP (`penetrant.itp`) for the chosen penetrant, calls
`gmx insert-molecules`, patches the topology, outputs `system_with_penetrant.gro` and
`topol_penetrant.top`.

**Force fields:**
- **O₂** — TraPPE 2-site (Potoff & Siepmann 2001): σ = 3.02 Å, ε/k_B = 49 K, bond = 1.21 Å
- **H₂O** — SPC/E (Berendsen 1987): r_OH = 1.0 Å, ∠HOH = 109.47°, q_O = −0.8476 e

### 04-diffusion-md
Brief EM to relax inserted molecules → NVT diffusion run → `gmx msd` on penetrant group →
Einstein relation D = slope / 6 (nm²/ps → cm²/s). Reports log-log slope as a diffusive-
regime diagnostic (slope ≈ 1 required for reliable D).

### 05-report
HTML report with:
- MSD log-log curve with slope = 1 reference line
- Computed D vs literature D (colour-coded pass/fail)
- Barrier ranking bar chart (log scale, longer = better barrier)

---

## Literature D reference values (cm²/s, ~23 °C)

| Resin | O₂ | H₂O |
|---|---|---|
| PET  | 3.4 × 10⁻¹⁰ | 5.0 × 10⁻¹² |
| LDPE | 4.5 × 10⁻⁷  | 1.5 × 10⁻⁸  |
| PP   | 2.0 × 10⁻⁸  | 3.0 × 10⁻¹⁰ |
| EVOH | 2.0 × 10⁻¹³ | 1.5 × 10⁻¹¹ |
| PA6  | 1.5 × 10⁻⁹  | 5.0 × 10⁻¹¹ |
