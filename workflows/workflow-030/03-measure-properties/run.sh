#!/bin/bash
set -euo pipefail

# ── parameters ────────────────────────────────────────────────────────────────
TARGET_TEMP_C=${PARAM_TEMPERATURE:-23.0}
PROD_TIME_PS=${PARAM_PROD_TIME_PS:-1000.0}
DEFORM_TIME_PS=${PARAM_DEFORM_TIME_PS:-50.0}

GMX=/usr/local/gromacs/avx2_256/bin/gmx

TARGET_TEMP_K=$(python3 -c "print(${TARGET_TEMP_C} + 273.15)")
PROD_NSTEPS=$(python3  -c "print(int(${PROD_TIME_PS}  / 0.001))")
DEFORM_NSTEPS=$(python3 -c "print(int(${DEFORM_TIME_PS} / 0.001))")

# Read initial box x-length (nm) from the equilibrated structure (last line of GRO)
L0_NM=$(python3 -c "
gro = open('inputs/equilibrated.gro').readlines()
print(gro[-1].split()[0])
")

# Deform rate targets 2% total strain: rate = 0.02 * L0 / DEFORM_TIME_PS  (nm/ps)
DEFORM_RATE=$(python3 -c "print(round(0.02 * ${L0_NM} / ${DEFORM_TIME_PS}, 6))")

echo "=== Measure Properties ==="
echo "  Target temp : ${TARGET_TEMP_K} K"
echo "  Production  : ${PROD_TIME_PS} ps (${PROD_NSTEPS} steps)"
echo "  Deformation : ${DEFORM_TIME_PS} ps (${DEFORM_NSTEPS} steps)"
echo "  Box L0      : ${L0_NM} nm"
echo "  Deform rate : ${DEFORM_RATE} nm/ps  → 2% total strain"
echo ""

mkdir -p outputs

cp inputs/equilibrated.gro .
cp inputs/topol.top        .

# ── Step 1: Production NPT ────────────────────────────────────────────────────
echo "=== Step 1: Production NPT (${PROD_TIME_PS} ps) ==="
sed -e "s/PROD_NSTEPS/${PROD_NSTEPS}/g" \
    -e "s/TARGET_TEMP/${TARGET_TEMP_K}/g" \
    mdp/prod.mdp > prod.mdp
$GMX grompp -f prod.mdp -c equilibrated.gro -p topol.top -o prod.tpr -maxwarn 5
$GMX mdrun -v -deffnm prod -ntmpi 1

echo "Density" | $GMX energy -f prod.edr -o outputs/density_prod.xvg 2>/dev/null

# ── Step 2: Uniaxial deformation NVT ─────────────────────────────────────────
echo ""
echo "=== Step 2: Uniaxial deformation NVT (${DEFORM_TIME_PS} ps, rate=${DEFORM_RATE} nm/ps) ==="
sed -e "s/DEFORM_NSTEPS/${DEFORM_NSTEPS}/g" \
    -e "s/TARGET_TEMP/${TARGET_TEMP_K}/g" \
    -e "s/DEFORM_RATE/${DEFORM_RATE}/g" \
    mdp/deform.mdp > deform.mdp
$GMX grompp -f deform.mdp -c prod.gro -t prod.cpt -p topol.top -o deform.tpr -maxwarn 5
$GMX mdrun -v -deffnm deform -ntmpi 1

echo "Pres-XX" | $GMX energy -f deform.edr -o outputs/stress_strain.xvg 2>/dev/null

# ── Step 3: Analysis ──────────────────────────────────────────────────────────
echo ""
echo "=== Step 3: Analysis ==="
TARGET_TEMP_K=${TARGET_TEMP_K} \
L0_NM=${L0_NM} \
DEFORM_RATE=${DEFORM_RATE} \
DEFORM_TIME_PS=${DEFORM_TIME_PS} \
PROD_TIME_PS=${PROD_TIME_PS} \
python3 analyze.py

echo ""
echo "Done — outputs: properties.json  density_prod.xvg  stress_strain.xvg"
