#!/usr/bin/env python3
"""
PRODIGY affinity scoring for filter-passing binder designs.

Chain convention: chain A = binder, chain B = receptor (RFdiffusion convention).
PRODIGY --selection A B: A = binder, B = receptor.

PRODIGY accuracy note: r=0.73, RMSE=1.89 kcal/mol on natural crystal structures.
Performance on predicted structures is lower and unbenchmarked. Treat ΔG values as
relative rankings within this campaign, not absolute affinity predictions.

Zero-passing-designs: if filter_report.json has no passing designs, write empty
output files and exit 0 — the workflow should not fail with an exception.
"""
import csv
import glob
import json
import os
import shutil
import subprocess
import sys

os.makedirs('./outputs/results/top10', exist_ok=True)

binder_chain = os.environ.get('PARAM_BINDER_CHAIN', 'A')
receptor_chain = os.environ.get('PARAM_RECEPTOR_CHAIN', 'B')
temperature = float(os.environ.get('PARAM_TEMPERATURE', '25.0'))
distance_cutoff = float(os.environ.get('PARAM_DISTANCE_CUTOFF', '5.5'))
top_n = int(os.environ.get('PARAM_TOP_N', '10'))
dg_cutoff = float(os.environ.get('PARAM_DG_CUTOFF', '-8.0'))
pymol_selection = os.environ.get('PARAM_PYMOL_SELECTION', 'true').lower() == 'true'

# Load filter report and collect passing designs
with open('./inputs/folded/filter_report.json') as f:
    filter_report = json.load(f)

filter_map = {r['design_id']: r for r in filter_report}
passing_ids = [r['design_id'] for r in filter_report if r['pass']]

print(f'Passing designs from filter step: {len(passing_ids)}/{len(filter_report)}', flush=True)

if not passing_ids:
    print('\nNo designs passed the fold-and-filter step.', flush=True)
    print('Suggestions to recover more designs:', flush=True)
    print('  - Increase max_bb_rmsd (currently set in node 04)', flush=True)
    print('  - Decrease min_iptm (currently set in node 04)', flush=True)
    print('  - Decrease min_plddt_binder (currently set in node 04)', flush=True)
    print('  - Increase max_pae_interaction (currently set in node 04)', flush=True)
    print('  - Generate more backbones (increase num_designs in node 02)', flush=True)
    with open('./outputs/results/prodigy_all_designs.json', 'w') as f:
        json.dump([], f)
    with open('./outputs/results/ranked_designs.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'design_id', 'dg', 'kd', 'iptm', 'plddt_binder',
            'bb_rmsd', 'pae_interaction', 'weak_binder',
        ])
        writer.writeheader()
    sys.exit(0)

# Locate passing PDB files
pdb_dir = './inputs/folded'
prodigy_results = []

for design_id in passing_ids:
    pdb_path = os.path.join(pdb_dir, f'{design_id}_rank001.pdb')
    if not os.path.exists(pdb_path):
        print(f'WARNING: PDB not found for passing design {design_id}: {pdb_path}', flush=True)
        continue

    cmd = [
        'prodigy',
        '--selection', binder_chain, receptor_chain,
        '--temperature', str(temperature),
        '--distance-cutoff', str(distance_cutoff),
        '-q',
        pdb_path,
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        stdout = result.stdout.strip()
    except subprocess.CalledProcessError as e:
        print(f'WARNING: PRODIGY failed for {design_id}: {e.stderr.strip()}', flush=True)
        continue

    # Parse tab-separated stdout: filename\tΔG\tKd
    parts = stdout.split()
    if len(parts) < 2:
        print(f'WARNING: Unexpected PRODIGY output for {design_id}: "{stdout}"', flush=True)
        continue

    try:
        dg = float(parts[-1])
        kd = None
    except ValueError:
        print(f'WARNING: Could not parse PRODIGY values for {design_id}: "{stdout}"', flush=True)
        continue

    filter_rec = filter_map.get(design_id, {})
    record = {
        'design_id': design_id,
        'dg': round(dg, 3),
        'kd': kd,
        'iptm': filter_rec.get('iptm'),
        'plddt_binder': filter_rec.get('plddt_binder'),
        'bb_rmsd': filter_rec.get('bb_rmsd'),
        'pae_interaction': filter_rec.get('pae_interaction'),
        'pdb_path': pdb_path,
    }
    prodigy_results.append(record)
    print(f'  {design_id}: ΔG={dg:.2f} kcal/mol', flush=True)

if not prodigy_results:
    print('WARNING: PRODIGY produced no results. Writing empty output files.', flush=True)
    with open('./outputs/results/prodigy_all_designs.json', 'w') as f:
        json.dump([], f)
    with open('./outputs/results/ranked_designs.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'design_id', 'dg', 'kd', 'iptm', 'plddt_binder',
            'bb_rmsd', 'pae_interaction', 'weak_binder',
        ])
        writer.writeheader()
    sys.exit(0)

# Sort by ΔG ascending (more negative = stronger predicted binder)
prodigy_results.sort(key=lambda r: r['dg'])

# Flag weak binders
for r in prodigy_results:
    r['weak_binder'] = r['dg'] > dg_cutoff

# Write full JSON
json_out = []
for r in prodigy_results:
    rec = dict(r)
    rec.pop('pdb_path', None)
    json_out.append(rec)
with open('./outputs/results/prodigy_all_designs.json', 'w') as f:
    json.dump(json_out, f, indent=2)

# Write CSV
csv_fields = ['design_id', 'dg', 'kd', 'iptm', 'plddt_binder', 'bb_rmsd', 'pae_interaction', 'weak_binder']
with open('./outputs/results/ranked_designs.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=csv_fields, extrasaction='ignore')
    writer.writeheader()
    writer.writerows(json_out)

# Copy top N PDBs
top_designs = prodigy_results[:top_n]
for r in top_designs:
    src = r['pdb_path']
    dest = f'./outputs/results/top10/{r["design_id"]}_rank001.pdb'
    shutil.copy(src, dest)

    if pymol_selection:
        pml_path = f'./outputs/results/top10/{r["design_id"]}.pml'
        pdb_basename = os.path.basename(dest)
        with open(pml_path, 'w') as pml:
            pml.write('# Run from the directory containing this file.\n')
            pml.write(f'load {pdb_basename}\n')
            pml.write(f'select binder, chain {binder_chain}\n')
            pml.write(f'select receptor, chain {receptor_chain}\n')
            pml.write('select interface, (binder within 5.0 of receptor) or (receptor within 5.0 of binder)\n')
            pml.write('show sticks, interface\n')
            pml.write('color cyan, binder\n')
            pml.write('color salmon, receptor\n')
            pml.write('zoom interface\n')

print(f'\nRanked {len(prodigy_results)} designs by ΔG', flush=True)
print(f'Top {len(top_designs)} exported to ./outputs/results/top10/', flush=True)
strong = sum(1 for r in prodigy_results if not r['weak_binder'])
print(f'Designs below ΔG cutoff ({dg_cutoff} kcal/mol): {strong}/{len(prodigy_results)}', flush=True)
