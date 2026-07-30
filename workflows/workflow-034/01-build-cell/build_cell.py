#!/usr/bin/env python3
"""Node 01 — Build amorphous polymer cell for GROMACS MD.

Pipeline:
  RDKit (oligomer 3D) → antechamber (GAFF2 + AM1-BCC) → parmchk2 (frcmod)
  → packmol (amorphous cell) → tleap (AmberTop) → acpype (GROMACS topology)
"""

import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

OUTDIR = Path("outputs")
OUTDIR.mkdir(exist_ok=True)

WORKDIR = Path("/tmp/build_cell")
WORKDIR.mkdir(exist_ok=True)

# ── parameters ────────────────────────────────────────────────────────────────
RESIN_TYPE    = os.environ.get("PARAM_RESIN_TYPE",    "PP").upper()
FIBER_LOADING = float(os.environ.get("PARAM_FIBER_LOADING", "0"))
N_CHAINS      = int(os.environ.get("PARAM_N_CHAINS",  "20"))
# Pack at 60% of target density so molecules fit comfortably; NPT in node 02 compresses to full density.
PACK_DENSITY_FRAC = 0.60
# Mass/density sizing alone can size a box smaller than a single extended chain
# (e.g. PA6 at low n_chains) — geometrically impossible for packmol to solve.
# Floor the box side at this multiple of the oligomer's own max pairwise atomic
# distance so the box can always fit at least one chain with clearance.
EXTENT_FLOOR_MULT = 1.5

# ── resin library ─────────────────────────────────────────────────────────────
# Oligomer SMILES: 5-mer end-capped chains — compact enough for packmol to solve
# quickly while preserving backbone connectivity for GAFF2 parameterisation.
# PP  — 5-mer:  CH3-[CH(CH3)-CH2]5-CH3
# PA6 — 5-mer:  H2N-[(CH2)5-CO-NH]4-(CH2)5-COOH
# PC  — 3-mer:  phenol-capped BPA-PC
# PET — 3-mer:  diol-capped
# ABS — 4-mer SAN approximation (acrylonitrile/styrene, no butadiene)
RESINS = {
    "PP": {
        "name":        "Polypropylene",
        "smiles":      "CC(C)" + "CC(C)" * 4 + "C",   # 5-mer
        "density_gcc": 0.855,
        "melt_temp_c": 200,
    },
    "PA6": {
        "name":        "Nylon-6 (PA6)",
        "smiles":      "NCCCCCC(=O)" * 5 + "O",        # 5-mer
        "density_gcc": 1.084,
        "melt_temp_c": 270,
    },
    "PC": {
        "name":        "Polycarbonate (BPA-PC)",
        "smiles":      (
            "OC(=O)Oc1ccc(C(C)(C)c2ccc(OC(=O)Oc3ccc"
            "(C(C)(C)c4ccc(O)cc4)cc3)cc2)cc1"
        ),
        "density_gcc": 1.20,
        "melt_temp_c": 250,
    },
    "PET": {
        "name":        "Polyethylene Terephthalate (PET)",
        "smiles":      (
            "OCCOC(=O)c1ccc(C(=O)OCCOC(=O)c2ccc"
            "(C(=O)OCCOC(=O)c3ccc(C(=O)O)cc3)cc2)cc1"
        ),
        "density_gcc": 1.335,
        "melt_temp_c": 280,
    },
    "ABS": {
        "name":        "ABS (SAN matrix approximation)",
        "smiles":      (
            "CC(C#N)CC(c1ccccc1)CC(C#N)CC(c1ccccc1)CC(C#N)C"
        ),
        "density_gcc": 1.05,
        "melt_temp_c": 240,
    },
}


def run(cmd, **kw):
    print(f"  $ {cmd}", flush=True)
    subprocess.run(cmd, shell=True, check=True, **kw)


def oligomer_extent_angstrom(mol) -> float:
    """Max pairwise atomic distance in the embedded 3D conformer (Å)."""
    conf = mol.GetConformer()
    positions = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    return max(
        positions[i].Distance(positions[j])
        for i in range(len(positions))
        for j in range(i + 1, len(positions))
    )


def build_oligomer(smiles: str, out_pdb: Path, out_sdf: Path) -> tuple[int, float]:
    """Generate 3D coordinates with ETKDG + MMFF, write PDB + SDF.
    Returns (atom count, max pairwise atomic distance in Å)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        sys.exit(f"ERROR: invalid SMILES for {RESIN_TYPE}")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
        sys.exit("ERROR: 3D embedding failed — reduce n_repeat or check SMILES")
    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    # SDF keeps RDKit's exact (Kekulized) bond orders for antechamber; PDB (no
    # bond-order field) is only used for packmol, which just needs coordinates.
    Chem.MolToMolFile(mol, str(out_sdf))
    Chem.MolToPDBFile(mol, str(out_pdb))
    n_atoms  = mol.GetNumAtoms()
    extent_a = oligomer_extent_angstrom(mol)
    print(f"  Oligomer: {n_atoms} atoms, {extent_a:.1f} Å extent → {out_pdb}")
    return n_atoms, extent_a


def gaff2_params(sdf: Path) -> tuple[Path, Path]:
    """antechamber + parmchk2 → GAFF2 mol2 and frcmod."""
    mol2   = WORKDIR / "mol.mol2"
    frcmod = WORKDIR / "mol.frcmod"
    # Use residue name UNL to match RDKit's default PDB output so tleap can
    # map atom types from the mol2 onto the packmol-packed cell.pdb. Feed the
    # SDF (not the PDB) so antechamber sees exact bond orders instead of
    # guessing them from geometry alone.
    run(
        f"antechamber -i {sdf} -fi mdl -o {mol2} -fo mol2 "
        f"-c bcc -s 2 -rn UNL -at gaff2",
        cwd=WORKDIR,
    )
    run(f"parmchk2 -i {mol2} -f mol2 -o {frcmod} -s gaff2", cwd=WORKDIR)
    return mol2, frcmod


def build_packed_cell(single_pdb: Path, n_chains: int, box_a: float) -> Path:
    """packmol → amorphous cell PDB."""
    packed_pdb = WORKDIR / "cell.pdb"
    inp = (
        f"tolerance 2.0\n"
        f"filetype pdb\n"
        f"output {packed_pdb}\n\n"
        f"structure {single_pdb}\n"
        f"  number {n_chains}\n"
        f"  inside box 0. 0. 0. {box_a:.2f} {box_a:.2f} {box_a:.2f}\n"
        f"end structure\n"
    )
    inp_file = WORKDIR / "packmol.inp"
    inp_file.write_text(inp)
    cmd = f"packmol < {inp_file}"
    print(f"  $ {cmd}", flush=True)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        forced_pdb = WORKDIR / "cell.pdb_FORCED"
        if "ENDED WITHOUT PERFECT PACKING" in result.stdout and forced_pdb.exists():
            # packmol hit its iteration budget without perfect packing but still
            # wrote a best-effort solution; packmol itself recommends using it
            # as an MD starting configuration, so fall back to it here rather
            # than treating the nonzero exit as fatal.
            print("  packmol: no perfect packing within tolerance — using forced solution", flush=True)
            shutil.copy(forced_pdb, packed_pdb)
        else:
            print(result.stdout[-2000:], file=sys.stderr)
            print(result.stderr[-2000:], file=sys.stderr)
            sys.exit(f"ERROR: packmol failed: {cmd}")
    return packed_pdb


def amber_to_gromacs(mol2: Path, frcmod: Path, packed_pdb: Path) -> None:
    """tleap (AmberTop) → acpype → GROMACS .gro / .top."""
    leap_in = WORKDIR / "tleap.in"
    leap_in.write_text(
        f"source leaprc.gaff2\n"
        f"UNL = loadmol2 {mol2}\n"
        f"loadamberparams {frcmod}\n"
        f"system = loadpdb {packed_pdb}\n"
        f"setbox system vdw 0\n"
        f"saveamberparm system {WORKDIR}/system.prmtop {WORKDIR}/system.inpcrd\n"
        f"quit\n"
    )
    run(f"tleap -f {leap_in}", cwd=WORKDIR)
    run(
        f"acpype -p {WORKDIR}/system.prmtop -x {WORKDIR}/system.inpcrd "
        f"-b system -o gmx",
        cwd=WORKDIR,
    )
    amb2gmx = WORKDIR / "system.amb2gmx"
    shutil.copy(amb2gmx / "system_GMX.gro", OUTDIR / "system.gro")
    shutil.copy(amb2gmx / "system_GMX.top", OUTDIR / "topol.top")
    shutil.copy(packed_pdb,                  OUTDIR / "cell.pdb")


def box_side_angstrom(n_chains: int, mw_g_per_mol: float, density_gcc: float,
                      pack_frac: float = 1.0) -> float:
    """Cubic box side in Å. pack_frac < 1 gives a looser initial cell so
    molecules fit without clashes; NPT in node 02 compresses to full density."""
    AVOGADRO = 6.02214076e23
    mass_g   = (n_chains * mw_g_per_mol) / AVOGADRO
    vol_cc   = mass_g / (density_gcc * pack_frac)
    return (vol_cc * 1e24) ** (1.0 / 3.0)   # cm³ → Å³ → Å side


def main():
    if RESIN_TYPE not in RESINS:
        sys.exit(f"ERROR: unknown resin '{RESIN_TYPE}'. Choose from {list(RESINS)}")

    cfg = RESINS[RESIN_TYPE]
    print(f"\n=== Build Cell: {cfg['name']} ===")
    print(f"  chains={N_CHAINS}  fiber={FIBER_LOADING} wt%  target density={cfg['density_gcc']} g/cc")
    if FIBER_LOADING > 0:
        print(f"  Note: glass-fiber correction ({FIBER_LOADING:.0f} wt%) applied in node 04")

    # 1. Oligomer 3D structure (PDB for packmol, SDF for antechamber)
    pdb_single = WORKDIR / "oligomer.pdb"
    sdf_single = WORKDIR / "oligomer.sdf"
    _, extent_a = build_oligomer(cfg["smiles"], pdb_single, sdf_single)

    # 2. Molecular weight
    mol = Chem.AddHs(Chem.MolFromSmiles(cfg["smiles"]))
    mw  = Descriptors.ExactMolWt(mol)

    # 3. Box size — mass/density sizing, floored against the oligomer's own
    # extent so the box can never end up geometrically smaller than a single
    # chain (silently unsolvable for packmol at low n_chains).
    box_a_density = box_side_angstrom(N_CHAINS, mw, cfg["density_gcc"], pack_frac=PACK_DENSITY_FRAC)
    box_a_floor   = EXTENT_FLOOR_MULT * extent_a
    box_a = max(box_a_density, box_a_floor)
    print(f"  MW/chain={mw:.1f} g/mol  box={box_a:.1f} Å  "
          f"(density-based={box_a_density:.1f} Å, extent-floor={box_a_floor:.1f} Å, "
          f"packed at {PACK_DENSITY_FRAC*100:.0f}% density)")

    # 4. GAFF2 parameters
    mol2, frcmod = gaff2_params(sdf_single)

    # 5. Pack amorphous cell
    packed_pdb = build_packed_cell(pdb_single, N_CHAINS, box_a)

    # 6. Convert to GROMACS
    amber_to_gromacs(mol2, frcmod, packed_pdb)

    # 7. Build report
    report = {
        "resin_type":          RESIN_TYPE,
        "resin_name":          cfg["name"],
        "n_chains":            N_CHAINS,
        "mw_chain_g_per_mol":  round(mw, 2),
        "density_gcc":         cfg["density_gcc"],
        "box_angstrom":        round(box_a, 2),
        "pack_density_frac":   PACK_DENSITY_FRAC,
        "melt_temp_c":         cfg["melt_temp_c"],
        "fiber_loading_wt_pct": FIBER_LOADING,
        "fiber_correction_applied": False,
    }
    (OUTDIR / "build_report.json").write_text(json.dumps(report, indent=2))
    print(f"\nDone — outputs: {[f.name for f in OUTDIR.iterdir()]}")


if __name__ == "__main__":
    main()
