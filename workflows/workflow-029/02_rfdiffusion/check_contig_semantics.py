#!/usr/bin/env python3.9
"""Assert what RFdiffusion's ContigMap does with a segmented receptor contig.

Node 02 joins the receptor's segments with '/', relying on RFdiffusion advancing
its positional index across each junction by the real size of the gap, and by 200
across a chain break. That arithmetic lives in RFdiffusion, not here, so nothing
in this repo would notice if a new image changed it.

Deliberately not part of the default test run: it needs the pinned RFdiffusion
image (no GPU). Run it inside that image when the image tag is bumped:

    python3.9 check_contig_semantics.py
"""
import sys

sys.path.insert(0, '/app/RFdiffusion')

from rfdiffusion.contigs import ContigMap

# A real node 01 record: 2B5I chain A, two breaks of 4 missing residues each.
SEGMENTS = ['A6-75', 'A80-98', 'A103-133']
BINDER_LENGTH = 80

# The receptor block ends in '/0', which is what marks it as the receptor rather
# than a chain to design; the binder follows as a free length range.
contig = f"{'/'.join(SEGMENTS)}/0 {BINDER_LENGTH}-{BINDER_LENGTH}"
print(f'Contig: [{contig}]', flush=True)

# ContigMap reads only parsed_pdb['pdb_idx'] on the contig-string path, to check
# every named residue exists, so the segments themselves stand in for the target.
pdb_idx = []
for segment in SEGMENTS:
    first, last = (int(number) for number in segment[1:].split('-'))
    pdb_idx += [(segment[0], number) for number in range(first, last + 1)]

contig_map = ContigMap({'pdb_idx': pdb_idx}, contigs=[contig])
ref = contig_map.ref
rf = contig_map.rf

binder = [i for i, entry in enumerate(ref) if entry == ('_', '_')]
receptor = [i for i, entry in enumerate(ref) if entry != ('_', '_')]
assert len(binder) == BINDER_LENGTH, f'binder is {len(binder)} residues, not {BINDER_LENGTH}'
print(f'Binder: {len(binder)} residues; receptor: {len(receptor)} residues', flush=True)

failures = []

# Across the binder/receptor boundary RFdiffusion inserts its 200-index chain
# break, so the step is 201: one position plus the break.
boundary_step = rf[receptor[0]] - rf[binder[-1]]
print(f'Binder -> receptor boundary step: {boundary_step}', flush=True)
if boundary_step != 201:
    failures.append(f'binder boundary step is {boundary_step}, expected 201')

# Inside the receptor every step must be 1 + the residues the gap swallowed, so
# a junction stays a gap rather than becoming a bond or a second chain.
for first, second in zip(receptor, receptor[1:]):
    missing = ref[second][1] - ref[first][1] - 1
    step = rf[second] - rf[first]
    expected = 1 + missing
    if missing:
        print(f'Junction {ref[first][0]}{ref[first][1]}|{ref[second][0]}{ref[second][1]}: '
              f'{missing} missing, index step {step} (expected {expected})', flush=True)
    if step != expected:
        failures.append(f'step {ref[first][1]}|{ref[second][1]} is {step}, expected {expected}')

if failures:
    for failure in failures:
        print(f'FAIL: {failure}', flush=True)
    sys.exit(1)

print('ContigMap semantics hold for this image', flush=True)
