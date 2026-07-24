#!/bin/bash
set -euo pipefail

PENETRANT=${PARAM_PENETRANT:-O2}
TARGET_TEMP_C=${PARAM_TEMPERATURE:-23.0}
N_INSERTIONS=${PARAM_N_INSERTIONS:-5000}
GMX=/usr/local/gromacs/avx2_256/bin/gmx

TARGET_TEMP_K=$(python3 -c "print(${TARGET_TEMP_C} + 273.15)")

echo "=== Solubility via Test Particle Insertion: ${PENETRANT} ==="
echo "  Temperature : ${TARGET_TEMP_K} K"
echo "  Insertions  : ${N_INSERTIONS} per trajectory frame"
echo ""

mkdir -p outputs
cp inputs/equilibrated.gro .
cp inputs/equil.xtc        .
cp inputs/topol.top        .

# ── Write the test-particle GRO + ITP ─────────────────────────────────────────
python3 write_test_particle.py "${PENETRANT}"

# ── Insert exactly one test particle (TPI relocates it randomly each trial;
#    this only gives grompp valid starting coordinates) ──────────────────────
$GMX insert-molecules \
  -f equilibrated.gro \
  -ci test_particle.gro \
  -nmol 1 \
  -o system_with_tp.gro \
  -try 1000 2>&1 | tail -5

# ── Append the test particle as the LAST molecule in the topology (required
#    for gmx tpi — it always test-inserts the final moleculetype) ────────────
python3 -c "
import re, pathlib
top = pathlib.Path('topol.top').read_text()
mol_name = 'O2' if '${PENETRANT}' == 'O2' else 'SOL'
itp_line = '#include \"test_particle.itp\"\n'
if itp_line not in top:
    top = re.sub(r'(\[ moleculetype \])', itp_line + r'\1', top, count=1)
top = top.rstrip() + f'\n{mol_name}                1\n'
pathlib.Path('topol_tpi.top').write_text(top)
"

# ── TPI mdp ────────────────────────────────────────────────────────────────────
sed -e "s/TPI_NSTEPS/${N_INSERTIONS}/g" \
    -e "s/TARGET_TEMP/${TARGET_TEMP_K}/g" \
    mdp/tpi.mdp > tpi.mdp

$GMX grompp -f tpi.mdp -c system_with_tp.gro -p topol_tpi.top -o tpi.tpr -maxwarn 10

# ── Rerun over the equilibration trajectory for better insertion statistics ──
echo ""
echo "=== Running TPI over equilibration trajectory ==="
$GMX mdrun -deffnm tpi -rerun equil.xtc -tpi outputs/tpi.xvg -ntmpi 1 2>&1 | tee tpi_stdout.log

N_FRAMES=$(grep -c "Reading frame" tpi_stdout.log || true)

# ── Compute solubility S from the TPI log ─────────────────────────────────────
python3 compute_solubility.py \
  --log      tpi.log \
  --report   outputs/solubility_report.json \
  --penetrant "${PENETRANT}" \
  --resin    "$(python3 -c "import json; print(json.load(open('inputs/build_report.json'))['resin_type'])")" \
  --temp     "${TARGET_TEMP_C}" \
  --n_frames "${N_FRAMES}" \
  --n_insertions_per_frame "${N_INSERTIONS}"

echo ""
echo "Done — outputs: solubility_report.json  tpi.xvg"
