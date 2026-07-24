#!/bin/bash
set -euo pipefail

PENETRANT=${PARAM_PENETRANT:-O2}
N_PENETRANT=${PARAM_N_PENETRANT:-5}
GMX=/usr/local/gromacs/avx2_256/bin/gmx

echo "=== Insert Penetrant: ${PENETRANT} (${N_PENETRANT} molecules) ==="

mkdir -p outputs
cp inputs/equilibrated.gro .
cp inputs/topol.top        .

# ── Write penetrant GRO and ITP ───────────────────────────────────────────────
python3 write_penetrant.py "${PENETRANT}"

# ── gmx insert-molecules ──────────────────────────────────────────────────────
$GMX insert-molecules \
  -f equilibrated.gro \
  -ci penetrant.gro \
  -nmol ${N_PENETRANT} \
  -o system_with_penetrant.gro \
  -try 1000 2>&1 | tail -5

ACTUAL=$(grep -c "^${PENETRANT}" system_with_penetrant.gro || true)
echo "  Inserted: $(grep 'Added' /dev/stdin <<< "$(cat /dev/null)" || true)"
echo "  Checking insertion..."
python3 -c "
import subprocess, sys
out = subprocess.run(['grep', '-c', '${PENETRANT}', 'system_with_penetrant.gro'],
                    capture_output=True, text=True)
n = int(out.stdout.strip()) if out.returncode == 0 else 0
print(f'  Lines with ${PENETRANT} in GRO: {n}')
"

# ── Update topology ───────────────────────────────────────────────────────────
python3 update_topology.py "${PENETRANT}" "${N_PENETRANT}"

# ── Write report ──────────────────────────────────────────────────────────────
python3 -c "
import json, pathlib
build = json.loads(pathlib.Path('inputs/build_report.json').read_text())
report = {
    'resin_type':   build['resin_type'],
    'penetrant':    '${PENETRANT}',
    'n_penetrant':  ${N_PENETRANT},
    'forcefield':   'TraPPE' if '${PENETRANT}' == 'O2' else 'SPC/E',
}
pathlib.Path('outputs/penetrant_report.json').write_text(json.dumps(report, indent=2))
"

cp system_with_penetrant.gro outputs/system_with_penetrant.gro
cp topol_penetrant.top       outputs/topol_penetrant.top

echo ""
echo "Done — outputs: system_with_penetrant.gro  topol_penetrant.top  penetrant_report.json"
