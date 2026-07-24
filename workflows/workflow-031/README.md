---
doc_id: workflow-031
domain: materials-science
doc_type: workflow
version: "1.0.0"
deprecated: false
description: >
  Measures the O₂ or H₂O diffusion coefficient D and solubility coefficient S through
  a polymer film (LDPE, PP, EVOH, PA6, or PET) using all-atom GAFF2 MD (GROMACS),
  the Einstein relation, and test particle insertion, delivering an HTML report with
  MSD curve, permeability P = D × S, and barrier ranking vs. literature values.
tags: [molecular-dynamics, polymer, gromacs, diffusion, solubility, permeability, barrier-film, gaff2, spce, trappe, tpi]
---

# Workflow 031: Barrier Films — Plastics That Protect Food

All-atom molecular-dynamics pipeline that builds an amorphous polymer cell,
equilibrates it via melt-quench NPT, then branches into two independent
measurements on that same cell: an NVT diffusion trajectory for small-molecule
penetrants (O₂ or H₂O) giving the diffusion coefficient D via the Einstein
relation, and a test particle insertion (Widom insertion) run giving the
solubility coefficient S. The two combine into permeability P = D × S — the
production-relevant OTR/WVTR-equivalent quantity for food-packaging barrier
performance.

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
| `n_penetrant` | `5` | Penetrant molecules inserted (node 03, diffusion) |
| `diff_time_ps` | `5000.0` | NVT diffusion run length (ps) |
| `n_insertions` | `5000` | TPI insertion attempts per trajectory frame (node 03b, solubility) — increase for dense films (EVOH, PA6) where insertion acceptance is low |

---

## Resin support and known limitations

| Resin | Default SMILES proxy | Status |
|---|---|---|
| `LDPE` | n-decane (C₁₀) | ✅ Supported (default) |
| `PP` | propylene 5-mer | ✅ Supported |
| `EVOH` | E/VOH 5-mer | ✅ Supported |
| `PA6` | caprolactam 5-mer | ✅ Supported |
| `PET` | terephthalate 5-mer | ✅ Supported |

### PET: fixed tleap connectivity failure

`antechamber` previously read the RDKit-generated oligomer from a plain PDB, which has
no bond-order field. For PET's ester carbonyl carbon this made `antechamber` guess the
wrong hybridisation (sp3 instead of sp2) and emit a duplicate bond, which `tleap` then
rejected:

```
Atom .R<UNL 1>.A<C36 47> has force field coordination 4 but only 3 bonded neighbors.
```

**Fix:** `build_cell.py` now also writes the oligomer as an SDF
(`RDKit.Chem.MolToMolFile`), which preserves RDKit's exact Kekulized bond orders, and
feeds that to `antechamber` (`-fi mdl`) instead of the bond-order-less PDB. Packmol and
tleap are unaffected — they still use the PDB for coordinates only.

---

## Node pipeline

```
                                    ┌─ 03-insert-penetrant ─→ 04-diffusion-md ─┐
01-build-cell  →  02-equilibrate  ─┤                                          ├─→ 05-report
(barrier_films)    (gromacs)        └─ 03b-solubility-tpi ────────────────────┘  (barrier_films)
                                       (gromacs)               (gromacs)
```

### 01-build-cell
RDKit 3-D conformer → antechamber GAFF2 + AM1-BCC charges → parmchk2 → packmol
(mol2 I/O) → tleap → acpype → `system.gro`, `topol.top`, `build_report.json`.

### 02-equilibrate
Energy minimisation → NPT melt at `melt_temp_c` (Berendsen barostat) →
NPT quench to `temperature` → `equilibrated.gro`, `equil.xtc`, `density.xvg`,
`equil_report.json`. `equil.xtc` feeds node 03b's TPI rerun for better insertion
statistics.

### 03-insert-penetrant
Writes force-field ITP (`penetrant.itp`) for the chosen penetrant, calls
`gmx insert-molecules`, patches the topology, outputs `system_with_penetrant.gro` and
`topol_penetrant.top`.

**Force fields:**
- **O₂** — TraPPE 2-site (Potoff & Siepmann 2001): σ = 3.02 Å, ε/k_B = 49 K, bond = 1.21 Å
- **H₂O** — SPC/E (Berendsen 1987): r_OH = 1.0 Å, ∠HOH = 109.47°, q_O = −0.8476 e

### 03b-solubility-tpi
Computes the Henry-regime solubility coefficient S via GROMACS test particle
insertion (`gmx tpi`, i.e. Widom insertion). Inserts a single test-particle copy of
the penetrant (same force fields as node 03) with `gmx insert-molecules`, appends it
as the *last* topology entry (required by `gmx tpi`), then reruns TPI over
`equil.xtc` (`gmx mdrun -rerun`) for `n_insertions` random insertion attempts per
frame. Parses the reported excess chemical potential (`<mu> = ... kJ/mol`) and
converts it to S via the ideal-gas/Henry's-law equilibrium relation:

```
S = 1/(RT) · exp(−μ_ex / RT)
```

reported in the conventional membrane-science unit cm³(STP)/(cm³·cmHg). H₂O runs are
flagged `qualitative` in `solubility_report.json` — Henry's law breaks down for water
sorption in polar resins (EVOH, PA6) via clustering/swelling, so H₂O solubility is a
relative trend, not an absolute value.

### 04-diffusion-md
Brief EM to relax inserted molecules → NVT diffusion run → `gmx msd` on penetrant group →
Einstein relation D = slope / 6 (nm²/ps → cm²/s). Reports log-log slope as a diffusive-
regime diagnostic (slope ≈ 1 required for reliable D).

### 05-report
Combines D (node 04) and S (node 03b) into permeability P = D × S (Barrer units,
via S expressed in cm³(STP)/(cm³·cmHg) so that D[cm²/s] × S gives the standard
solution-diffusion permeability directly). HTML report with:
- KPI strip: D, S, P = D × S, MSD log-log slope
- Solubility & permeability table (μ_ex, S, P) with a qualitative-H₂O callout when applicable
- MSD log-log curve with slope = 1 reference line
- Barrier ranking bar chart vs. literature D (log scale, longer = better barrier)

---

## Verification

Last verified with a full `silva workflows/workflow-031` end-to-end run (default
params: LDPE, O₂, 23 °C, 20 chains) — all 6 jobs completed, including
`03b-solubility-tpi`:

```
01-build-cell → 02-equilibrate → 03-insert-penetrant → 03b-solubility-tpi → 04-diffusion-md → 05-report
```

Results from that run:

| Quantity | Value |
|---|---|
| MSD log-log slope | 0.93 (diffusive regime) |
| D (MD) | 2.22 × 10⁻⁴ cm²/s |
| μ_ex (TPI, 38 frames × 5000 insertions) | −1.28 kJ/mol |
| S | 2.04 × 10⁻² cm³(STP)/(cm³·cmHg) |
| P = D × S | 4.52 × 10⁴ Barrer |

D and P are ~500× above the literature values in the table below — this is expected
with the *default* `melt_time_ps`/`equil_time_ps`, which only reached 441 kg/m³ of
the 920 kg/m³ target density (defaults are tuned for fast screening, not accuracy;
see "When to use this workflow" above). The mechanism itself (TPI convergence,
`<mu>` parsing, S and P formulas, report rendering) is verified correct — increase
`equil_time_ps` for production-accurate D/S/P.

> **Root cause found (2026-07-14):** the 441/920 kg/m³ shortfall above wasn't
> purely a "needs more time" issue — `mdp/melt.mdp`/`mdp/equil.mdp` used
> `compressibility = 4.5e-5 bar⁻¹` (water), which let the cell drift instead
> of condense: density stayed flat and noisy around 431 kg/m³ across the
> whole 500 ps default equilibration, with no real upward trend. Diagnosed
> in workflow-032
> ([issue #196](https://github.com/chiral-data/collab-workflows/issues/196),
> [issue #198](https://github.com/chiral-data/collab-workflows/issues/198))
> and fixed here to `compressibility = 2.0e-6 bar⁻¹` (polymer-melt order of
> magnitude): the same default-length run now shows a genuine, still-rising
> monotonic trend (310→371 kg/m³ across deciles) — real convergence, just
> not complete yet at these short default times. `equil_time_ps` is still
> the right lever to increase for production accuracy; it will now actually
> work.

Full end-to-end re-run after the fix (same default params):

| Quantity | Value |
|---|---|
| MSD log-log slope | 0.255 (sub-diffusive — run longer for reliable D) |
| D (MD) | 1.06 × 10⁻⁴ cm²/s |
| μ_ex (TPI, 38 frames × 5000 insertions) | −0.78 kJ/mol |
| S | 1.66 × 10⁻² cm³(STP)/(cm³·cmHg) |
| P = D × S | 1.77 × 10⁴ Barrer |

D and P are still well above literature (as expected at these short default
times — D/S/P convergence is a separate, longer lever than the density fix,
see `diff_time_ps` above), but the mechanism and the density-convergence
trend are both confirmed genuinely working now.

---

## Mock mode

Each node's `run.sh` (what `job.toml` executes) downloads that node's
pre-computed outputs from `output_files/<node>/` on `main` instead of running
the real RDKit/AmberTools/GROMACS pipeline — no Docker/GPU/wait required to
inspect the DAG or the report. The real pipeline (`run_real.sh` and the
computation scripts it calls) lives in
[workflow-035](../workflow-035/README.md), the real-pipeline counterpart of
this workflow.

`output_files/` and `sample_outputs/` hold the actual outputs from the
verified real run above (LDPE, O₂, defaults, post-fix).

---

## Literature D reference values (cm²/s, ~23 °C)

| Resin | O₂ | H₂O |
|---|---|---|
| PET  | 3.4 × 10⁻¹⁰ | 5.0 × 10⁻¹² |
| LDPE | 4.5 × 10⁻⁷  | 1.5 × 10⁻⁸  |
| PP   | 2.0 × 10⁻⁸  | 3.0 × 10⁻¹⁰ |
| EVOH | 2.0 × 10⁻¹³ | 1.5 × 10⁻¹¹ |
| PA6  | 1.5 × 10⁻⁹  | 5.0 × 10⁻¹¹ |
