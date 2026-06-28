#!/bin/bash
set -euo pipefail

TARGET_TEMP_C=${PARAM_TEMPERATURE:-23.0}
DIFF_TIME_PS=${PARAM_DIFF_TIME_PS:-5000.0}
PENETRANT=${PARAM_PENETRANT:-O2}

GMX=/usr/local/gromacs/avx2_256/bin/gmx

TARGET_TEMP_K=$(python3 -c "print(${TARGET_TEMP_C} + 273.15)")
DIFF_NSTEPS=$(python3 -c "print(int(${DIFF_TIME_PS} / 0.002))")
MOL_NAME=$(python3 -c "print('O2' if '${PENETRANT}' == 'O2' else 'SOL')")

echo "=== Diffusion MD ==="
echo "  Penetrant  : ${PENETRANT} (${MOL_NAME})"
echo "  Temperature: ${TARGET_TEMP_K} K"
echo "  Run length : ${DIFF_TIME_PS} ps (${DIFF_NSTEPS} steps, dt=2 fs)"
echo ""

mkdir -p outputs
cp inputs/system_with_penetrant.gro .
cp inputs/topol_penetrant.top       .

# Need penetrant.itp in working dir (referenced by topology)
python3 - <<'PYEOF'
import os, sys
sys.path.insert(0, '/workspace')
penetrant = os.environ.get('PARAM_PENETRANT', 'O2')
exec(open('/workspace/03-insert-penetrant/write_penetrant.py').read() if False else '')

import math
if penetrant == 'O2':
    itp = """
[ atomtypes ]
OT  8  15.9994  0.000  A  0.30200  0.40740

[ moleculetype ]
O2  3

[ atoms ]
1  OT  1  O2  O1  1  0.0  15.9994
2  OT  1  O2  O2  2  0.0  15.9994

[ bonds ]
1  2  1  0.1210  40000.0
"""
else:
    cos_half = math.cos(math.radians(109.47/2))
    sin_half = math.sin(math.radians(109.47/2))
    itp = """
[ atomtypes ]
OW  8  15.9994  -0.8476  A  0.31660  0.65017
HW  1   1.0080   0.4238  A  0.00000  0.00000

[ moleculetype ]
SOL  2

[ atoms ]
1  OW  1  SOL  OW   1  -0.8476  15.9994
2  HW  1  SOL  HW1  1   0.4238   1.0080
3  HW  1  SOL  HW2  1   0.4238   1.0080

[ bonds ]
1  2  1  0.10000  345000.0
1  3  1  0.10000  345000.0

[ angles ]
2  1  3  1  109.47  383.0
"""
open('penetrant.itp','w').write(itp)
print(f'  Wrote penetrant.itp for {penetrant}')
PYEOF

# ── Energy minimisation of penetrant-inserted cell ────────────────────────────
echo "=== Step 1: Brief EM to relax inserted molecules ==="
cat > em_fast.mdp << 'EOF'
integrator  = steep
emtol       = 1000.0
emstep      = 0.001
nsteps      = 2000
cutoff-scheme = Verlet
coulombtype = PME
rcoulomb    = 1.2
rvdw        = 1.2
pbc         = xyz
constraints = none
EOF

$GMX grompp -f em_fast.mdp -c system_with_penetrant.gro -p topol_penetrant.top \
            -o em.tpr -maxwarn 10
$GMX mdrun -v -deffnm em -ntmpi 1

# ── Diffusion NVT ─────────────────────────────────────────────────────────────
echo ""
echo "=== Step 2: NVT diffusion run (${DIFF_TIME_PS} ps) ==="
sed -e "s/DIFF_NSTEPS/${DIFF_NSTEPS}/g" \
    -e "s/TARGET_TEMP/${TARGET_TEMP_K}/g" \
    mdp/nvt_diffusion.mdp > nvt_diffusion.mdp

$GMX grompp -f nvt_diffusion.mdp -c em.gro -p topol_penetrant.top \
            -o diffusion.tpr -maxwarn 10
$GMX mdrun -v -deffnm diffusion -ntmpi 1

# ── MSD analysis ──────────────────────────────────────────────────────────────
echo ""
echo "=== Step 3: MSD analysis ==="
echo "${MOL_NAME}" | $GMX msd \
  -f diffusion.xtc \
  -s diffusion.tpr \
  -o outputs/msd.xvg \
  -mol 2>/dev/null

# ── Compute D via Einstein relation ───────────────────────────────────────────
python3 compute_diffusion.py \
  --msd    outputs/msd.xvg \
  --report outputs/diffusion_report.json \
  --penetrant "${PENETRANT}" \
  --resin "$(python3 -c "import json; print(json.load(open('inputs/build_report.json'))['resin_type'])")" \
  --temp   "${TARGET_TEMP_C}"

echo ""
echo "Done — outputs: msd.xvg  diffusion_report.json"
