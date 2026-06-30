#!/usr/bin/env python3
"""
Run ProteinMPNN on each RFdiffusion backbone.

RFdiffusion convention: chain A = binder, chain B = receptor.
ProteinMPNN designs chain A while holding chain B fixed.

Design ID naming: RFdiffusion outputs design_0.pdb, design_1.pdb, etc.
This script renames outputs to zero-padded IDs: design_000.fasta, design_001.fasta, etc.
This canonical padded ID is used in all downstream nodes.
"""
import glob
import os
import shutil
import subprocess
import sys
import tempfile

os.makedirs('./outputs/sequences', exist_ok=True)

num_seq_per_target = int(os.environ.get('PARAM_NUM_SEQ_PER_TARGET', '8'))
sampling_temp = float(os.environ.get('PARAM_SAMPLING_TEMP', '0.1'))
chain_to_design = os.environ.get('PARAM_CHAIN_TO_DESIGN', 'A')
use_soluble_model = os.environ.get('PARAM_USE_SOLUBLE_MODEL', 'true').lower() == 'true'
use_gpu = os.environ.get('PARAM_USE_GPU', 'true').lower() == 'true'

if use_gpu:
    result = subprocess.run(['nvidia-smi'], capture_output=True)
    if result.returncode != 0:
        print('ERROR: GPU requested but nvidia-smi failed.', flush=True)
        sys.exit(1)
    print('GPU check passed', flush=True)

# Glob and sort backbone PDBs by integer index for stable ordering
backbone_paths = sorted(
    glob.glob('./inputs/backbones/design_*.pdb'),
    key=lambda p: int(os.path.splitext(os.path.basename(p))[0].split('_')[1])
)

if not backbone_paths:
    print('ERROR: No backbone PDB files found in ./inputs/backbones/', flush=True)
    sys.exit(1)

print(f'Found {len(backbone_paths)} backbone PDB(s)', flush=True)
print(f'Sequences per backbone: {num_seq_per_target}', flush=True)
print(f'Sampling temperature: {sampling_temp}', flush=True)
print(f'SolubleMPNN: {use_soluble_model}', flush=True)

fasta_count = 0

for padded_idx, backbone_path in enumerate(backbone_paths):
    design_id = f'design_{padded_idx:03d}'
    print(f'Processing {os.path.basename(backbone_path)} → {design_id}', flush=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        out_dir = os.path.join(tmpdir, 'mpnn_out')
        os.makedirs(out_dir)

        cmd = [
            'python', '/opt/ProteinMPNN/protein_mpnn_run.py',
            '--pdb_path', backbone_path,
            '--pdb_path_chains', chain_to_design,
            '--out_folder', out_dir,
            '--num_seq_per_target', str(num_seq_per_target),
            '--sampling_temp', str(sampling_temp),
            '--batch_size', '1',
        ]
        if use_soluble_model:
            cmd.append('--use_soluble_model')

        print(f'  Running: {" ".join(cmd)}', flush=True)
        subprocess.run(cmd, check=True)

        # ProteinMPNN writes FASTA to {out_dir}/seqs/
        seqs_dir = os.path.join(out_dir, 'seqs')
        fa_files = glob.glob(os.path.join(seqs_dir, '*.fa'))
        if not fa_files:
            print(f'  WARNING: No .fa files produced for {backbone_path}', flush=True)
            continue

        dest = f'./outputs/sequences/{design_id}.fasta'
        shutil.copy(fa_files[0], dest)
        fasta_count += 1
        print(f'  Wrote {dest}', flush=True)

print(f'Total FASTA files written: {fasta_count}', flush=True)
if fasta_count == 0:
    print('ERROR: No FASTA files were produced', flush=True)
    sys.exit(1)
