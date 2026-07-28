# Node 01 — Build Amorphous Cell

Generates a GROMACS-ready amorphous polymer cell from a resin type and chain count.

**Pipeline:**
```
RDKit (oligomer 3D) → antechamber (GAFF2 + AM1-BCC) → parmchk2 (frcmod)
→ packmol (amorphous cell) → tleap (AmberTop) → acpype (GROMACS topology)
```

---

## Inputs

No upstream node inputs. The cell is generated entirely from parameters.

## Outputs

| File | Description |
|------|-------------|
| `system.gro` | GROMACS coordinate file for the packed cell |
| `topol.top` | GROMACS topology with GAFF2 force field parameters |
| `cell.pdb` | Packed amorphous cell PDB (from packmol) |
| `build_report.json` | Build metadata: resin, MW, box size, chain count |

---

## Parameters

| Parameter | Default | Options | Description |
|-----------|---------|---------|-------------|
| `PARAM_RESIN_TYPE` | `PP` | PP, PA6, PC, PET, ABS | Polymer resin |
| `PARAM_FIBER_LOADING` | `0` | 0, 15, 30 | Glass-fiber loading (wt%) — recorded in report; rule-of-mixtures correction applied in node 04 |
| `PARAM_N_CHAINS` | `20` | integer | Number of oligomer chains in the cell |

---

## Resin library

| Resin | Oligomer | MW/chain (g/mol) | Amorphous density (g/cc) | Melt temp (°C) |
|-------|----------|-----------------|--------------------------|----------------|
| PP | 5-mer isotactic | 226 | 0.855 | 200 |
| PA6 | 5-mer | 571 | 1.084 | 270 |
| PC | 3-mer (BPA) | 763 | 1.20 | 250 |
| PET | 3-mer | 576 | 1.335 | 280 |
| ABS | 4-mer (SAN approx.) | 741 | 1.05 | 240 |

The cell is packed at **60% of target density** so molecules fit without clashes. Node 02 (NPT melt-quench) compresses to full amorphous density.

---

## Docker image

```bash
docker build -t polymer_md:2026_07_28 apps/polymer_md_2026_07_28/
```

Tools: `rdkit`, `ambertools` (antechamber, parmchk2, tleap), `packmol`, `acpype`.

---

## Local test

```bash
mkdir -p /tmp/test-01/outputs
docker run --rm \
  -v /tmp/test-01/outputs:/workspace/outputs \
  -v $(pwd):/workspace \
  -w /workspace \
  -e PARAM_RESIN_TYPE=PP \
  -e PARAM_FIBER_LOADING=0 \
  -e PARAM_N_CHAINS=20 \
  polymer_md:2026_07_28 bash run.sh
```

Expected output: `system.gro` (1000 atoms for PP 20-chain), `topol.top`, `cell.pdb`, `build_report.json`.
