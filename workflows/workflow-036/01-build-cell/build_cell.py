#!/usr/bin/env python3
"""Node 01 — Build amorphous polymer cell for GROMACS MD.

Pipeline:
  RDKit (oligomer 3D) -> antechamber (GAFF2 + AM1-BCC) -> parmchk2 (frcmod)
  -> packmol (amorphous cell) -> tleap (AmberTop) -> acpype (GROMACS topology)

Crystallinity (low/medium/high) is approximated by varying the initial
packing density fraction — a trend-only proxy, not a real semi-crystalline
lattice (building a true crystalline unit cell is out of scope; see README).
"""

import json
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
RESIN_TYPE    = os.environ.get("PARAM_RESIN_TYPE",    "PPS").upper()
CRYSTALLINITY = os.environ.get("PARAM_CRYSTALLINITY", "medium").lower()
N_CHAINS      = int(os.environ.get("PARAM_N_CHAINS",  "20"))

# Packing density fraction of each resin's target amorphous density.
# Higher crystallinity -> tighter initial packing (denser, more ordered
# proxy). NPT in node 02 relaxes/compresses further from this starting point.
CRYSTALLINITY_PACK_FRAC = {
    "low":    0.55,
    "medium": 0.65,
    "high":   0.75,
}

# ── resin library ─────────────────────────────────────────────────────────────
# Oligomer SMILES: short end-capped chains — compact enough for packmol to
# solve quickly while preserving backbone connectivity for GAFF2
# parameterisation. density_gcc is the amorphous-phase target used to size
# the packing box; melt_temp_c is the NPT melt-stage target in node 02 (a
# randomisation temperature, not necessarily the resin's true Tm).
# Literature Tg values are carried through for the node 03/04 comparison —
# they are reference numbers, not computed by this pipeline.
RESINS = {
    "PPS": {
        "name":          "Polyphenylene Sulfide (PPS)",
        "smiles":        "Sc1ccc(Sc2ccc(Sc3ccc(Sc4ccc(Sc5ccccc5)cc4)cc3)cc2)cc1",  # 5-ring
        "density_gcc":   1.35,
        "melt_temp_c":   285,
        "literature_tg_c": 88,
    },
    "PA66": {
        "name":          "Nylon-6,6 (PA66)",
        "smiles":        ("NCCCCCCNC(=O)CCCCC(=O)" * 3) + "O",  # 3-mer — heavier repeat unit than PA6, kept short so antechamber's AM1-BCC SCF stays fast
        "density_gcc":   1.14,
        "melt_temp_c":   260,
        "literature_tg_c": 50,
    },
    "PBT": {
        "name":          "Polybutylene Terephthalate (PBT)",
        "smiles": (
            "OCCCCOC(=O)c1ccc(C(=O)OCCCCOC(=O)c2ccc"
            "(C(=O)OCCCCOC(=O)c3ccc(C(=O)O)cc3)cc2)cc1"
        ),  # 3-mer
        "density_gcc":   1.31,
        "melt_temp_c":   225,
        "literature_tg_c": 45,
    },
    "PEEK": {
        "name":          "Polyether Ether Ketone (PEEK)",
        "smiles": (
            "Oc1ccc(Oc2ccc(C(=O)c3ccc(Oc4ccc(Oc5ccc"
            "(C(=O)c6ccccc6)cc5)cc4)cc3)cc2)cc1"
        ),  # 3-mer ether-ether-ketone
        "density_gcc":   1.30,
        "melt_temp_c":   343,
        "literature_tg_c": 143,
    },
    "PP": {
        "name":          "Polypropylene (PP, reference)",
        "smiles":        "CC(C)" + "CC(C)" * 4 + "C",   # 5-mer, same proxy as workflow-030
        "density_gcc":   0.855,
        "melt_temp_c":   200,
        "literature_tg_c": -10,
    },
}


def run(cmd, **kw):
    print(f"  $ {cmd}", flush=True)
    subprocess.run(cmd, shell=True, check=True, **kw)


def build_oligomer(smiles: str, out_pdb: Path, out_sdf: Path) -> int:
    """Generate 3D coordinates with ETKDG + MMFF, write PDB + SDF. Returns atom count.

    The SDF keeps RDKit's exact Kekulized bond orders for antechamber; a plain
    PDB has no bond-order field, which makes antechamber guess hybridisation
    from geometry alone. For aromatic/conjugated backbones (PPS, PBT, PEEK)
    that guess can be wrong (e.g. an aromatic-ring carbon read as sp3), which
    tleap then rejects with a "coordination 4 but only 3 bonded neighbors"
    error — the same failure workflow-031 hit and fixed for PET. The PDB is
    still used for packmol/tleap, which only need atomic coordinates.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        sys.exit(f"ERROR: invalid SMILES for {RESIN_TYPE}")
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.useRandomCoords = True
    embedded = False
    for seed in range(10):
        params.randomSeed = seed
        if AllChem.EmbedMolecule(mol, params) == 0:
            embedded = True
            break
    if not embedded:
        sys.exit("ERROR: 3D embedding failed after 10 attempts — check SMILES")

    AllChem.MMFFOptimizeMolecule(mol, maxIters=2000)
    Chem.MolToMolFile(mol, str(out_sdf))
    Chem.MolToPDBFile(mol, str(out_pdb))
    n_atoms = mol.GetNumAtoms()
    print(f"  Oligomer: {n_atoms} atoms -> {out_pdb}")
    return n_atoms


def gaff2_params(sdf: Path) -> tuple[Path, Path]:
    """antechamber + parmchk2 -> GAFF2 mol2 and frcmod."""
    mol2   = WORKDIR / "mol.mol2"
    frcmod = WORKDIR / "mol.frcmod"
    # Use residue name UNL to match RDKit's default PDB output so tleap can
    # map atom types from the mol2 onto the packmol-packed cell.pdb. Feed the
    # SDF (not the PDB) so antechamber sees exact bond orders — see
    # build_oligomer() docstring.
    run(
        f"antechamber -i {sdf} -fi mdl -o {mol2} -fo mol2 "
        f"-c bcc -s 2 -rn UNL -at gaff2",
        cwd=WORKDIR,
    )
    run(f"parmchk2 -i {mol2} -f mol2 -o {frcmod} -s gaff2", cwd=WORKDIR)
    return mol2, frcmod


def build_packed_cell(single_pdb: Path, n_chains: int, box_a: float) -> Path:
    """packmol -> amorphous cell PDB."""
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
    run(f"packmol < {inp_file}")
    return packed_pdb


def amber_to_gromacs(mol2: Path, frcmod: Path, packed_pdb: Path) -> None:
    """tleap (AmberTop) -> acpype -> GROMACS .gro / .top."""
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
    molecules fit without clashes; NPT in node 02 compresses from there."""
    AVOGADRO = 6.02214076e23
    mass_g   = (n_chains * mw_g_per_mol) / AVOGADRO
    vol_cc   = mass_g / (density_gcc * pack_frac)
    return (vol_cc * 1e24) ** (1.0 / 3.0)   # cm³ -> Å³ -> Å side


def main():
    if RESIN_TYPE not in RESINS:
        sys.exit(f"ERROR: unknown resin '{RESIN_TYPE}'. Choose from {list(RESINS)}")
    if CRYSTALLINITY not in CRYSTALLINITY_PACK_FRAC:
        sys.exit(f"ERROR: unknown crystallinity '{CRYSTALLINITY}'. "
                  f"Choose from {list(CRYSTALLINITY_PACK_FRAC)}")

    cfg = RESINS[RESIN_TYPE]
    pack_frac = CRYSTALLINITY_PACK_FRAC[CRYSTALLINITY]
    print(f"\n=== Build Cell: {cfg['name']} ===")
    print(f"  chains={N_CHAINS}  crystallinity={CRYSTALLINITY} (pack_frac={pack_frac})"
          f"  target density={cfg['density_gcc']} g/cc")

    # 1. Oligomer 3D structure (PDB for packmol, SDF for antechamber)
    pdb_single = WORKDIR / "oligomer.pdb"
    sdf_single = WORKDIR / "oligomer.sdf"
    build_oligomer(cfg["smiles"], pdb_single, sdf_single)

    # 2. Molecular weight
    mol = Chem.AddHs(Chem.MolFromSmiles(cfg["smiles"]))
    mw  = Descriptors.ExactMolWt(mol)

    # 3. Box size
    box_a = box_side_angstrom(N_CHAINS, mw, cfg["density_gcc"], pack_frac=pack_frac)
    print(f"  MW/chain={mw:.1f} g/mol  box={box_a:.1f} Å  (packed at {pack_frac*100:.0f}% density)")

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
        "crystallinity":       CRYSTALLINITY,
        "pack_density_frac":   pack_frac,
        "crystallinity_note":  ("Approximated via initial packing density fraction — "
                                 "trend comparison only, not a real semi-crystalline lattice."),
        "melt_temp_c":         cfg["melt_temp_c"],
        "literature_tg_c":     cfg["literature_tg_c"],
    }
    (OUTDIR / "build_report.json").write_text(json.dumps(report, indent=2))
    print(f"\nDone — outputs: {[f.name for f in OUTDIR.iterdir()]}")


if __name__ == "__main__":
    main()
