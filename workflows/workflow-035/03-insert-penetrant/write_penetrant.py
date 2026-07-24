"""
Write GRO + ITP files for O2 (TraPPE 2-site) or H2O (SPC/E).
Called by run.sh before gmx insert-molecules.
"""
import sys

PENETRANT = sys.argv[1] if len(sys.argv) > 1 else "O2"

# ── O₂ — TraPPE 2-site model ─────────────────────────────────────────────────
# Potoff & Siepmann (2001) J. AIChE 47:1676
# σ=3.02 Å, ε/kB=49 K → ε=0.4074 kJ/mol, bond=1.21 Å, partial charges ±0
O2_GRO = """\
O2 penetrant — TraPPE 2-site
    2
    1O2    O1    1   0.000   0.000   0.000
    1O2    O2    2   0.000   0.000   0.121
   0.50000   0.50000   0.50000
"""

O2_ITP = """\
; O2 TraPPE 2-site (Potoff & Siepmann 2001)
[ atomtypes ]
; name  at.num  mass      charge  ptype  sigma(nm)  epsilon(kJ/mol)
OT      8       15.9994   0.000   A      0.30200    0.40740

[ moleculetype ]
O2  3

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge   mass
  1   OT    1      O2       O1    1     0.0      15.9994
  2   OT    1      O2       O2    2     0.0      15.9994

[ bonds ]
; i  j  funct  r0(nm)  kb(kJ/mol/nm²)
  1  2  1      0.1210  40000.0
"""

# ── H₂O — SPC/E model ────────────────────────────────────────────────────────
# Berendsen et al. (1987); q_O=-0.8476, q_H=+0.4238, r_OH=0.1 nm, HOH=109.47°
import math
cos_half = math.cos(math.radians(109.47 / 2))
sin_half = math.sin(math.radians(109.47 / 2))
r_OH = 0.1  # nm
Hx = r_OH * sin_half
Hz = r_OH * cos_half

H2O_GRO = f"""\
H2O penetrant — SPC/E
    3
    1SOL    OW    1   0.000   0.000   0.000
    1SOL   HW1    2  {Hx:.4f}   0.000  {Hz:.4f}
    1SOL   HW2    3  {-Hx:.4f}  0.000  {Hz:.4f}
   0.50000   0.50000   0.50000
"""

H2O_ITP = """\
; H2O SPC/E model (Berendsen 1987)
[ atomtypes ]
; name  at.num  mass      charge   ptype  sigma(nm)  epsilon(kJ/mol)
OW      8       15.9994  -0.8476   A      0.31660    0.65017
HW      1        1.0080   0.4238   A      0.00000    0.00000

[ moleculetype ]
SOL  2

[ atoms ]
; nr  type  resnr  residue  atom  cgnr  charge    mass
  1   OW    1      SOL      OW    1     -0.8476   15.9994
  2   HW    1      SOL     HW1    1      0.4238    1.0080
  3   HW    1      SOL     HW2    1      0.4238    1.0080

[ bonds ]
  1  2  1  0.10000  345000.0
  1  3  1  0.10000  345000.0

[ angles ]
  2  1  3  1  109.47  383.0
"""

if PENETRANT == "O2":
    open("penetrant.gro", "w").write(O2_GRO)
    open("penetrant.itp", "w").write(O2_ITP)
    print(f"  Wrote O2 TraPPE GRO + ITP")
elif PENETRANT in ("H2O", "H2O"):
    open("penetrant.gro", "w").write(H2O_GRO)
    open("penetrant.itp", "w").write(H2O_ITP)
    print(f"  Wrote H2O SPC/E GRO + ITP")
else:
    sys.exit(f"Unknown penetrant '{PENETRANT}'. Choose O2 or H2O.")
