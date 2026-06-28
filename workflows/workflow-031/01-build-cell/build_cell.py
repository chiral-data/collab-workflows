"""
Build an amorphous polymer cell for barrier-film MD.
Pipeline: RDKit 5-mer → antechamber GAFF2 + AM1-BCC → parmchk2 → packmol → tleap → acpype
"""
import json, math, os, pathlib, shutil, subprocess, sys, tempfile

RESIN_TYPE    = os.environ.get("PARAM_RESIN_TYPE",    "PET")
FIBER_LOADING = 0.0   # no GF in barrier films
N_CHAINS      = int(os.environ.get("PARAM_N_CHAINS", 20))
PACK_DENSITY_FRAC = 0.60

RESIN_LIBRARY = {
    # smiles (5-mer), density g/cc, melt_temp_c, full name
    "PET":  {
        "smiles":   "OCC(=O)c1ccc(cc1)C(=O)OCC(=O)c1ccc(cc1)C(=O)OCC(=O)c1ccc(cc1)C(=O)OCC(=O)c1ccc(cc1)C(=O)OCC",
        "density":  1.335,
        "melt_c":   280,
        "name":     "Polyethylene terephthalate",
    },
    "LDPE": {
        # low-density PE — linear 5-mer of ethylene + short branch via 1-butene co-unit
        "smiles":   "CCCCCCCCCC",   # n-decane as tractable LDPE proxy (C10)
        "density":  0.920,
        "melt_c":   115,
        "name":     "Low-density polyethylene",
    },
    "PP":   {
        "smiles":   "CC(C)CC(C)CC(C)CC(C)CC(C)C",
        "density":  0.855,
        "melt_c":   200,
        "name":     "Polypropylene",
    },
    "EVOH": {
        # ethylene-vinyl alcohol — alternating E/VOH 5-mer proxy
        "smiles":   "CCOC(O)CC(O)CC(O)CC(O)CC(O)C",
        "density":  1.19,
        "melt_c":   190,
        "name":     "Ethylene-vinyl alcohol copolymer",
    },
    "PA6":  {
        "smiles":   "NCCCCCC(=O)NCCCCCC(=O)NCCCCCC(=O)NCCCCCC(=O)NCCCCCC(=O)O",
        "density":  1.084,
        "melt_c":   270,
        "name":     "Polyamide 6 (Nylon-6)",
    },
}

if RESIN_TYPE not in RESIN_LIBRARY:
    sys.exit(f"Unknown resin '{RESIN_TYPE}'. Choose from: {list(RESIN_LIBRARY)}")

resin    = RESIN_LIBRARY[RESIN_TYPE]
SMILES   = resin["smiles"]
DENSITY  = resin["density"]
MELT_C   = resin["melt_c"]

print(f"=== Build Cell: {resin['name']} ===")
print(f"  chains={N_CHAINS}  target density={DENSITY} g/cc")


def run(cmd, cwd=None):
    r = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=cwd)
    if r.returncode != 0:
        print(r.stdout[-2000:], file=sys.stderr)
        print(r.stderr[-2000:], file=sys.stderr)
        sys.exit(f"FAILED: {cmd}")
    return r.stdout


def box_side_angstrom(n_chains, smiles, density_gcc, pack_frac=PACK_DENSITY_FRAC):
    from rdkit.Chem import Descriptors, MolFromSmiles, AddHs
    mol  = AddHs(MolFromSmiles(smiles))
    mw   = Descriptors.MolWt(mol)
    mass_g   = n_chains * mw / 6.022e23
    vol_cc   = mass_g / (density_gcc * pack_frac)
    side_cm  = vol_cc ** (1/3)
    return side_cm * 1e8   # cm → Å


def build_oligomer_pdb(smiles, out_pdb):
    from rdkit.Chem import MolFromSmiles, AddHs
    from rdkit.Chem.AllChem import EmbedMolecule, MMFFOptimizeMolecule, ETKDGv3
    from rdkit.Chem import MolToPDBFile
    mol = AddHs(MolFromSmiles(smiles))
    ps  = ETKDGv3(); ps.randomSeed = 42
    EmbedMolecule(mol, ps)
    MMFFOptimizeMolecule(mol)
    MolToPDBFile(mol, str(out_pdb))
    n_atoms = mol.GetNumAtoms()
    print(f"  Oligomer: {n_atoms} atoms → {out_pdb}")
    return n_atoms


WORKDIR = pathlib.Path(tempfile.mkdtemp(prefix="build_barrier_"))

# 1. RDKit → PDB
oligomer_pdb = WORKDIR / "oligomer.pdb"
build_oligomer_pdb(SMILES, oligomer_pdb)

# 2. antechamber GAFF2 + AM1-BCC charges
mol2 = WORKDIR / "mol.mol2"
print(f"  $ antechamber ...")
run(f"antechamber -i {oligomer_pdb} -fi pdb -o {mol2} -fo mol2 "
    f"-c bcc -s 2 -rn UNL -at gaff2", cwd=WORKDIR)

# 3. parmchk2
frcmod = WORKDIR / "mol.frcmod"
run(f"parmchk2 -i {mol2} -f mol2 -o {frcmod} -s gaff2", cwd=WORKDIR)

# 4. packmol
box_a = box_side_angstrom(N_CHAINS, SMILES, DENSITY)
box_a = max(box_a, 20.0)   # minimum sensible box

cell_pdb = WORKDIR / "cell.pdb"
packmol_inp = WORKDIR / "packmol.inp"
packmol_inp.write_text(f"""
tolerance 2.0
output {cell_pdb}
filetype pdb
seed 1234567

structure {oligomer_pdb}
  number {N_CHAINS}
  inside box 0. 0. 0. {box_a:.3f} {box_a:.3f} {box_a:.3f}
end structure
""")
print(f"  $ packmol  (box={box_a:.1f} Å, {N_CHAINS} chains @ {PACK_DENSITY_FRAC*100:.0f}% density)")
run(f"packmol < {packmol_inp}", cwd=WORKDIR)

# 5. tleap → Amber topology
prmtop   = WORKDIR / "system.prmtop"
inpcrd   = WORKDIR / "system.inpcrd"
tleap_in = WORKDIR / "tleap.in"
tleap_in.write_text(f"""
source leaprc.gaff2
UNL = loadmol2 {mol2}
loadamberparams {frcmod}
sys = loadpdb {cell_pdb}
setbox sys vdw
saveamberparm sys {prmtop} {inpcrd}
quit
""")
print("  $ tleap ...")
out = run(f"tleap -f {tleap_in}", cwd=WORKDIR)
if "Errors = 0" not in out:
    for line in out.splitlines():
        if "Error" in line or "error" in line:
            print(f"  tleap: {line}", file=sys.stderr)

# 6. acpype → GROMACS
print("  $ acpype ...")
run(f"acpype -p {prmtop} -x {inpcrd} -b system -o gmx", cwd=WORKDIR)

# 7. Collect outputs
out_dir  = pathlib.Path("outputs")
out_dir.mkdir(exist_ok=True)
amb2gmx  = WORKDIR / "system.amb2gmx"
shutil.copy(amb2gmx / "system_GMX.gro", out_dir / "system.gro")
shutil.copy(amb2gmx / "system_GMX.top", out_dir / "topol.top")
shutil.copy(cell_pdb,                   out_dir / "cell.pdb")

report = {
    "resin_type":       RESIN_TYPE,
    "resin_name":       resin["name"],
    "n_chains":         N_CHAINS,
    "density_gcc":      DENSITY,
    "melt_temp_c":      MELT_C,
    "box_angstrom":     round(box_a, 2),
    "pack_density_frac": PACK_DENSITY_FRAC,
}
(out_dir / "build_report.json").write_text(json.dumps(report, indent=2))

print(f"Done — outputs: system.gro  topol.top  cell.pdb  build_report.json")
