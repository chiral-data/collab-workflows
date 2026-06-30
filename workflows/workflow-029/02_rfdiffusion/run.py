#!/usr/bin/env python3
import glob
import os
import subprocess
import sys

os.makedirs('./outputs/backbones', exist_ok=True)

num_designs = int(os.environ.get('PARAM_NUM_DESIGNS', '100'))
binder_length = int(os.environ.get('PARAM_BINDER_LENGTH', '80'))
hotspot_residues = os.environ.get('PARAM_HOTSPOT_RESIDUES', '').strip()
noise_scale = float(os.environ.get('PARAM_NOISE_SCALE', '1.0'))
num_diffusion_steps = int(os.environ.get('PARAM_NUM_DIFFUSION_STEPS', '50'))
use_gpu = os.environ.get('PARAM_USE_GPU', 'true').lower() == 'true'

if use_gpu:
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

print(f'Receptor length (chain A residues): {receptor_len}', flush=True)
print(f'Binder length: {binder_length}', flush=True)
print(f'Number of designs: {num_designs}', flush=True)

# Contig: chain A is always hardcoded — node 01 guarantees chain A output
contig = f'[A1-{receptor_len}/0 {binder_length}-{binder_length}]'

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

if hotspot_residues:
    # Format: "A55,A58,A102" → "[A55,A58,A102]"
    hotspot_arg = '[' + hotspot_residues.replace(' ', '') + ']'
    cmd.append(f'ppi.hotspot_res={hotspot_arg}')
    print(f'Hotspot residues: {hotspot_arg}', flush=True)
else:
    print('No hotspot residues specified — RFdiffusion will select automatically', flush=True)

print(f'Running: {" ".join(cmd)}', flush=True)
subprocess.run(cmd, check=True)

pdb_files = glob.glob('./outputs/backbones/design_*.pdb')
print(f'Generated {len(pdb_files)} backbone PDB(s) in ./outputs/backbones/', flush=True)
