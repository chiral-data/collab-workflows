#!/usr/bin/env python3
import glob
import os
import subprocess
import sys


def read_receptor_segments(path):
    """The receptor's contiguous segments as node 01 found them, or None."""
    with open(path) as f:
        for line in f:
            if line.startswith('REMARK 999 SEGMENTS'):
                return line[len('REMARK 999 SEGMENTS'):].strip()
            if not line.startswith('REMARK'):
                break
    return None


def derive_receptor_segments(path):
    """Fallback for a target carrying no segment record: split chain A on every
    gap in its numbering. RFdiffusion indexes CA records only, so they define the
    ranges a contig may name."""
    numbers = []
    seen = set()
    with open(path) as f:
        for line in f:
            if line.startswith('ATOM') and line[21] == 'A' and line[12:16].strip() == 'CA':
                if line[26] != ' ':
                    continue
                number = int(line[22:26])
                if number not in seen:
                    seen.add(number)
                    numbers.append(number)
    if not numbers:
        return None
    numbers.sort()
    segments = [[numbers[0], numbers[0]]]
    for number in numbers[1:]:
        if number == segments[-1][1] + 1:
            segments[-1][1] = number
        else:
            segments.append([number, number])
    return '/'.join(f'A{first}-{last}' for first, last in segments)


def read_hotspot_mapping(path):
    """Hotspot IDs as node 01 resolved them, or None if it wrote no mapping."""
    with open(path) as f:
        for line in f:
            if line.startswith('REMARK 999 HOTSPOTS'):
                return line[len('REMARK 999 HOTSPOTS'):].strip()
            if not line.startswith('REMARK'):
                break
    return None


def chain_residue_counts(lines):
    """Residues per chain letter, counting CA records, in first-seen order."""
    counts = {}
    for line in lines:
        if line.startswith('ATOM') and len(line) > 26 and line[12:16].strip() == 'CA':
            counts.setdefault(line[21], set()).add(line[22:27])
    return {chain: len(residues) for chain, residues in counts.items()}


def swap_binder_to_chain_a(path, binder_length):
    """Make the generated binder chain A and the receptor chain B.

    The contig puts the receptor first, so RFdiffusion emits it as chain A and
    the binder as chain B — the reverse of what nodes 03-05 assume. Each line is
    rewritten from the original, so the two letters cannot collide mid-swap.

    The swap is blind, so it is guarded first. Two ways it can be handed
    something it would corrupt:

    Re-running this node over a non-empty backbones/ directory. The swap is not
    idempotent, so a second pass over an already-swapped file puts the receptor
    back on chain A. This happens in practice — it is what the guard caught first.

    RFdiffusion emitting the receptor as more than one chain, which a junction
    wide enough to read as a chain break would do, since '/' advances its index
    by the size of the gap. Whether any real numbering gap reaches that width is
    untested; check_contig_semantics.py has a fixture for it.

    Either way chain B is not the binder, and flipping the letters hands node 03
    a piece of the receptor to redesign. Nothing downstream would notice: node 03
    never checks the binder's length, and node 04 compares the folded chain A
    against backbone chain A, which are equal by construction in that case. So
    the mislabel has to be caught here or not at all.
    """
    with open(path) as f:
        lines = f.readlines()

    counts = chain_residue_counts(lines)
    if sorted(counts) != ['A', 'B']:
        print(f'ERROR: {os.path.basename(path)} has chains {sorted(counts)}, expected exactly '
              "['A', 'B'] (receptor then binder). The receptor was probably emitted as more "
              'than one chain, which would make the binder relabelling assign the wrong '
              'chain. Residue counts per chain: '
              + ', '.join(f'{c}={n}' for c, n in sorted(counts.items())), flush=True)
        sys.exit(1)
    if counts['B'] != binder_length:
        print(f'ERROR: {os.path.basename(path)} chain B has {counts["B"]} residues, expected '
              f'{binder_length}. Chain B should be the generated binder; a mismatch means the '
              'contig order and the emitted chain order disagree, so relabelling would mark '
              'the receptor as the binder.', flush=True)
        sys.exit(1)

    with open(path, 'w') as f:
        for line in lines:
            if line.startswith(('ATOM', 'HETATM', 'TER')) and len(line) > 21 and line[21] in 'AB':
                line = line[:21] + ('B' if line[21] == 'A' else 'A') + line[22:]
            f.write(line)


os.makedirs('./outputs/backbones', exist_ok=True)

num_designs = int(os.environ.get('PARAM_NUM_DESIGNS', '100'))
binder_length = int(os.environ.get('PARAM_BINDER_LENGTH', '80'))
hotspot_residues = os.environ.get('PARAM_HOTSPOT_RESIDUES', '').strip()
noise_scale = float(os.environ.get('PARAM_NOISE_SCALE', '1.0'))
num_diffusion_steps = int(os.environ.get('PARAM_NUM_DIFFUSION_STEPS', '50'))

result = subprocess.run(['nvidia-smi'], capture_output=True)
if result.returncode != 0:
    print('ERROR: GPU requested but nvidia-smi failed. Ensure a CUDA-capable GPU is available '
          'and the NVIDIA Container Toolkit is installed.', flush=True)
    sys.exit(1)
print('GPU check passed', flush=True)

# Count receptor residues from chain A of validated_target.pdb (ATOM lines only)
receptor_len = 0
seen_residues = set()
with open('./inputs/validated_target.pdb') as f:
    for line in f:
        if line.startswith('ATOM') and line[21] == 'A':
            res_id = line[22:26].strip()
            if res_id not in seen_residues:
                seen_residues.add(res_id)
                receptor_len += 1

if receptor_len == 0:
    print('ERROR: No ATOM records found for chain A in validated_target.pdb', flush=True)
    sys.exit(1)

mapped_hotspots = read_hotspot_mapping('./inputs/validated_target.pdb')

segments = read_receptor_segments('./inputs/validated_target.pdb')
if segments is not None:
    segment_source = 'recorded by node 01'
else:
    segments = derive_receptor_segments('./inputs/validated_target.pdb')
    segment_source = 'derived from the target'

if not segments:
    print('ERROR: Could not determine the receptor segments of validated_target.pdb', flush=True)
    sys.exit(1)

print(f'Receptor length (chain A residues): {receptor_len}', flush=True)
print(f'Receptor segments ({segment_source}): {segments}', flush=True)
print(f'Binder length: {binder_length}', flush=True)
print(f'Number of designs: {num_designs}', flush=True)

# Chain A is hardcoded — node 01 guarantees chain A output. Segments are joined
# with '/', which advances RFdiffusion's positional index by the real size of
# each gap; '/0' would instead declare a 200-index break, i.e. a separate chain.
# check_contig_semantics.py asserts that arithmetic against the pinned image;
# re-run it if the image tag changes, since nothing here would notice a drift.
contig = f'[{segments}/0 {binder_length}-{binder_length}]'

cmd = [
    'python3.9', '/app/RFdiffusion/scripts/run_inference.py',
    f'inference.input_pdb=./inputs/validated_target.pdb',
    f'contigmap.contigs={contig}',
    f'inference.output_prefix=./outputs/backbones/design',
    f'inference.num_designs={num_designs}',
    f'inference.model_directory_path=/app/RFdiffusion/models',
    f'inference.schedule_directory_path=/tmp/schedules',
    f'denoiser.noise_scale_ca={noise_scale}',
    f'diffuser.T={num_diffusion_steps}',
]

if mapped_hotspots is not None:
    # Node 01 resolved these against the coordinates it wrote; RFdiffusion matches
    # hotspots against that same numbering
    source = 'mapped by node 01'
    hotspots = mapped_hotspots
else:
    source = 'raw parameter — target carries no hotspot mapping'
    hotspots = hotspot_residues.replace(' ', '')

if hotspots:
    hotspot_arg = '[' + hotspots + ']'
    cmd.append(f'ppi.hotspot_res={hotspot_arg}')
    print(f'Hotspot residues ({source}): {hotspot_arg}', flush=True)
else:
    print('No hotspot residues specified — RFdiffusion will select automatically', flush=True)

print(f'Running: {" ".join(cmd)}', flush=True)
subprocess.run(cmd, check=True)

pdb_files = glob.glob('./outputs/backbones/design_*.pdb')

for pdb_file in sorted(pdb_files):
    swap_binder_to_chain_a(pdb_file, binder_length)

print(f'Generated {len(pdb_files)} backbone PDB(s) in ./outputs/backbones/ '
      f'(binder relabelled to chain A, receptor to chain B)', flush=True)
