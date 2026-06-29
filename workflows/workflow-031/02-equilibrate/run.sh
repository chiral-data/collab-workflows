#!/bin/bash
set -euo pipefail

# ── parameters ────────────────────────────────────────────────────────────────
TARGET_TEMP_C=${PARAM_TEMPERATURE:-23.0}          # °C → convert to K below
EM_STEPS=${PARAM_EM_STEPS:-5000}
MELT_TIME_PS=${PARAM_MELT_TIME_PS:-200.0}
EQUIL_TIME_PS=${PARAM_EQUIL_TIME_PS:-500.0}

GMX=/usr/local/gromacs/avx2_256/bin/gmx

# Convert °C to K
TARGET_TEMP_K=$(python3 -c "print(${TARGET_TEMP_C} + 273.15)")

# Read melt temperature from node 01 build report
MELT_TEMP_K=$(python3 -c "
import json, pathlib
r = json.loads(pathlib.Path('inputs/build_report.json').read_text())
print(r['melt_temp_c'] + 273.15)
")

# nsteps from time (dt = 0.001 ps)
MELT_NSTEPS=$(python3  -c "print(int(${MELT_TIME_PS}  / 0.001))")
EQUIL_NSTEPS=$(python3 -c "print(int(${EQUIL_TIME_PS} / 0.001))")

echo "=== Polymer Melt-Quench Equilibration ==="
echo "  Melt temp  : ${MELT_TEMP_K} K"
echo "  Target temp: ${TARGET_TEMP_K} K"
echo "  EM steps   : ${EM_STEPS}"
echo "  Melt NPT   : ${MELT_TIME_PS} ps (${MELT_NSTEPS} steps)"
echo "  Equil NPT  : ${EQUIL_TIME_PS} ps (${EQUIL_NSTEPS} steps)"
echo ""

mkdir -p outputs

cp inputs/system.gro .
cp inputs/topol.top  .

# ── Step 1: Energy minimisation ───────────────────────────────────────────────
echo "=== Step 1: Energy minimisation ==="
sed "s/EM_STEPS/${EM_STEPS}/g" mdp/em.mdp > em.mdp
$GMX grompp -f em.mdp -c system.gro -p topol.top -o em.tpr -maxwarn 5
$GMX mdrun -v -deffnm em -ntmpi 1

# ── Step 2: NPT melt ──────────────────────────────────────────────────────────
echo ""
echo "=== Step 2: NPT at melt temperature (${MELT_TEMP_K} K) ==="
sed -e "s/MELT_NSTEPS/${MELT_NSTEPS}/g" \
    -e "s/MELT_TEMP/${MELT_TEMP_K}/g" \
    mdp/melt.mdp > melt.mdp
$GMX grompp -f melt.mdp -c em.gro -p topol.top -o melt.tpr -maxwarn 5
$GMX mdrun -v -deffnm melt -ntmpi 1

# ── Step 3: NPT equilibration at target temperature ───────────────────────────
echo ""
echo "=== Step 3: NPT equilibration at target temperature (${TARGET_TEMP_K} K) ==="
sed -e "s/EQUIL_NSTEPS/${EQUIL_NSTEPS}/g" \
    -e "s/TARGET_TEMP/${TARGET_TEMP_K}/g" \
    mdp/equil.mdp > equil.mdp
$GMX grompp -f equil.mdp -c melt.gro -t melt.cpt -p topol.top -o equil.tpr -maxwarn 5
$GMX mdrun -v -deffnm equil -ntmpi 1

# ── Extract density time series ───────────────────────────────────────────────
echo ""
echo "=== Extracting density ==="
echo "Density" | $GMX energy -f equil.edr -o outputs/density.xvg 2>/dev/null

# ── Write report ──────────────────────────────────────────────────────────────
python3 -c "
import pathlib, json
xvg = pathlib.Path('outputs/density.xvg').read_text()
values = [float(l.split()[1]) for l in xvg.splitlines() if l and not l.startswith(('#','@'))]
avg = sum(values[-len(values)//3:]) / max(len(values)//3, 1)
build = json.loads(pathlib.Path('inputs/build_report.json').read_text())
report = {
    'target_temp_k':     ${TARGET_TEMP_K},
    'melt_temp_k':       ${MELT_TEMP_K},
    'equil_time_ps':     ${EQUIL_TIME_PS},
    'avg_density_kg_m3': round(avg, 1),
    'resin_type':        build['resin_type'],
}
pathlib.Path('outputs/equil_report.json').write_text(json.dumps(report, indent=2))
print(f'  Average density (last third): {avg:.1f} kg/m3')
"

cp equil.gro outputs/equilibrated.gro

echo ""
echo "Done — outputs: equilibrated.gro  density.xvg  equil_report.json"
