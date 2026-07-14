#!/bin/bash
set -euo pipefail

# ── parameters ────────────────────────────────────────────────────────────────
EM_STEPS=${PARAM_EM_STEPS:-5000}
MELT_TIME_PS=${PARAM_MELT_TIME_PS:-150.0}
STAGE_TIME_PS=${PARAM_STAGE_TIME_PS:-100.0}

GMX=/usr/local/gromacs/avx2_256/bin/gmx

# Read melt temperature from node 01 build report
MELT_TEMP_C=$(python3 -c "
import json, pathlib
r = json.loads(pathlib.Path('inputs/build_report.json').read_text())
print(r['melt_temp_c'])
")
MELT_TEMP_K=$(python3 -c "print(${MELT_TEMP_C} + 273.15)")

EM_NSTEPS=$(python3 -c "print(int(${EM_STEPS}))")
MELT_NSTEPS=$(python3  -c "print(int(${MELT_TIME_PS}  / 0.001))")
STAGE_NSTEPS=$(python3 -c "print(int(${STAGE_TIME_PS} / 0.001))")

echo "=== Melt-Quench Cooling Series ==="
echo "  Melt temp   : ${MELT_TEMP_K} K (${MELT_TEMP_C} C)"
echo "  EM steps    : ${EM_NSTEPS}"
echo "  Melt NPT    : ${MELT_TIME_PS} ps (${MELT_NSTEPS} steps)"
echo "  Quench ladder: 200 / 150 / 80 / 25 C, ${STAGE_TIME_PS} ps each (${STAGE_NSTEPS} steps)"
echo ""

mkdir -p outputs

cp inputs/system.gro .
cp inputs/topol.top  .

# ── Step 1: Energy minimisation ───────────────────────────────────────────────
echo "=== Step 1: Energy minimisation ==="
sed "s/EM_STEPS/${EM_NSTEPS}/g" mdp/em.mdp > em.mdp
$GMX grompp -f em.mdp -c system.gro -p topol.top -o em.tpr -maxwarn 5
$GMX mdrun -v -deffnm em -ntmpi 1

# ── Step 2: NPT melt ──────────────────────────────────────────────────────────
echo ""
echo "=== Step 2: NPT melt (${MELT_TEMP_K} K) ==="
sed -e "s/MELT_NSTEPS/${MELT_NSTEPS}/g" \
    -e "s/MELT_TEMP/${MELT_TEMP_K}/g" \
    mdp/melt.mdp > melt.mdp
$GMX grompp -f melt.mdp -c em.gro -p topol.top -o melt.tpr -maxwarn 5
$GMX mdrun -v -deffnm melt -ntmpi 1
echo "Density" | $GMX energy -f melt.edr -o outputs/density_melt.xvg 2>/dev/null

# ── Step 3: chained NPT quench ladder ─────────────────────────────────────────
QUENCH_TEMPS_C=(200 150 80 25)
QUENCH_NAMES=(q200 q150 q80 q25)

PREV_GRO=melt.gro
PREV_CPT=melt.cpt

for i in "${!QUENCH_TEMPS_C[@]}"; do
    TEMP_C=${QUENCH_TEMPS_C[$i]}
    NAME=${QUENCH_NAMES[$i]}
    TEMP_K=$(python3 -c "print(${TEMP_C} + 273.15)")

    echo ""
    echo "=== Step 3.$((i+1)): NPT quench to ${TEMP_C} C (${TEMP_K} K) — ${NAME} ==="
    sed -e "s/QUENCH_NSTEPS/${STAGE_NSTEPS}/g" \
        -e "s/QUENCH_TEMP/${TEMP_K}/g" \
        mdp/quench.mdp > "${NAME}.mdp"
    $GMX grompp -f "${NAME}.mdp" -c "${PREV_GRO}" -t "${PREV_CPT}" -p topol.top -o "${NAME}.tpr" -maxwarn 5
    $GMX mdrun -v -deffnm "${NAME}" -ntmpi 1
    echo "Density" | $GMX energy -f "${NAME}.edr" -o "outputs/density_${NAME}.xvg" 2>/dev/null

    PREV_GRO="${NAME}.gro"
    PREV_CPT="${NAME}.cpt"
done

cp "${PREV_GRO}" outputs/equilibrated.gro

# ── Build cooling_series.json ─────────────────────────────────────────────────
echo ""
echo "=== Extracting density-vs-temperature series ==="
python3 -c "
import json, pathlib

def last_third_avg(xvg_path):
    vals = [float(l.split()[1]) for l in pathlib.Path(xvg_path).read_text().splitlines()
            if l and not l.startswith(('#', '@'))]
    n = max(len(vals) // 3, 1)
    return sum(vals[-n:]) / n

stages = [('melt', ${MELT_TEMP_C})] + list(zip(
    ['${QUENCH_NAMES[0]}', '${QUENCH_NAMES[1]}', '${QUENCH_NAMES[2]}', '${QUENCH_NAMES[3]}'],
    [${QUENCH_TEMPS_C[0]}, ${QUENCH_TEMPS_C[1]}, ${QUENCH_TEMPS_C[2]}, ${QUENCH_TEMPS_C[3]}],
))

series = []
for name, temp_c in stages:
    avg_density = last_third_avg(f'outputs/density_{name}.xvg')
    series.append({
        'stage':                      name,
        'temp_c':                     float(temp_c),
        'avg_density_kg_m3':          round(avg_density, 2),
        'avg_specific_volume_cm3_g':  round(1000.0 / avg_density, 6) if avg_density > 0 else None,
    })

build = json.loads(pathlib.Path('inputs/build_report.json').read_text())
report = {
    'resin_type':     build['resin_type'],
    'crystallinity':  build['crystallinity'],
    'stage_time_ps':  ${STAGE_TIME_PS},
    'series':         series,
}
pathlib.Path('outputs/cooling_series.json').write_text(json.dumps(report, indent=2))

for s in series:
    print(f\"  {s['stage']:>6}  {s['temp_c']:6.1f} C  {s['avg_density_kg_m3']:8.1f} kg/m3  \"
          f\"{s['avg_specific_volume_cm3_g']:.4f} cm3/g\")
"

echo ""
echo "Done — outputs: cooling_series.json  equilibrated.gro"
