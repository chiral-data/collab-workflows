#!/usr/bin/env python3
"""
ColabFold fold-and-filter for binder designs.

Key conventions:
- Chain A = binder, chain B = receptor (relabelled by node 02)
- ProteinMPNN FASTA files contain num_seq_per_target records each.
  We iterate every record, creating one ColabFold job per sequence.
  Design IDs: design_001_seq_0, design_001_seq_1, ...
- PAE and pLDDT extracted from *_scores_rank_001_*.json (ColabFold 1.6.1).
  No separate PAE file exists in this version.
- RMSD alignment uses sequence position, not residue numbers.
  Backbone reference = chain A of the RFdiffusion output PDB.
"""
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from Bio import PDB
from Bio.PDB import PDBParser, PPBuilder, Superimposer

os.makedirs('./outputs/folded', exist_ok=True)

num_recycle = int(os.environ.get('PARAM_NUM_RECYCLE', '6'))
msa_mode = os.environ.get('PARAM_MSA_MODE', 'mmseqs2_uniref_env')
min_iptm = float(os.environ.get('PARAM_MIN_IPTM', '0.6'))
min_plddt_binder = float(os.environ.get('PARAM_MIN_PLDDT_BINDER', '80.0'))
max_bb_rmsd = float(os.environ.get('PARAM_MAX_BB_RMSD', '1.5'))
max_pae_interaction = float(os.environ.get('PARAM_MAX_PAE_INTERACTION', '10.0'))

if not os.path.exists('/dev/nvidia0'):
    print('ERROR: GPU requested but /dev/nvidia0 not found.', flush=True)
    sys.exit(1)
print('GPU check passed', flush=True)

print(f'Filters: iPTM>={min_iptm}, pLDDT_binder>={min_plddt_binder}, '
      f'RMSD<={max_bb_rmsd}, PAE_interaction<={max_pae_interaction}', flush=True)


def read_segment_records(pdb_path):
    """(first, last, count) per receptor segment, as node 01 recorded them."""
    segments = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('REMARK 999 SEGMENT '):
                fields = line.split()
                segments.append((int(fields[5]), int(fields[6]), int(fields[7])))
            elif not line.startswith('REMARK'):
                break
    return segments


def read_hotspot_numbers(pdb_path):
    """Residue numbers of the hotspots node 01 resolved, in its output numbering."""
    numbers = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith('REMARK 999 HOTSPOTS'):
                for hotspot in line[len('REMARK 999 HOTSPOTS'):].strip().split(','):
                    hotspot = hotspot.strip().lstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZ')
                    if hotspot:
                        numbers.append(int(hotspot))
            elif not line.startswith('REMARK'):
                break
    return numbers


def choose_receptor_segment(segments, hotspots):
    """The segment carrying the epitope, or the largest when none is named."""
    def rank(segment):
        first, last, count = segment
        return (sum(1 for h in hotspots if first <= h <= last), count)

    return max(segments, key=rank)


def get_receptor_sequence(pdb_path, bounds=None):
    """Receptor sequence from validated_target.pdb chain A using PPBuilder.

    PPBuilder splits on the C-N distance and the fragments are concatenated, so a
    receptor spanning a chain break would be folded as though the gap were a
    peptide bond. bounds restricts the sequence to one contiguous segment.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('receptor', pdb_path)
    ppb = PPBuilder()
    seqs = []
    for pp in ppb.build_peptides(structure[0]['A']):
        sequence = str(pp.get_sequence())
        if bounds is None:
            seqs.append(sequence)
            continue
        first, last = bounds
        seqs.append(''.join(letter for res, letter in zip(pp, sequence)
                            if first <= res.id[1] <= last))
    return ''.join(seqs)


def get_ca_atoms_by_position(chain):
    """Return list of Cα atoms from ATOM residues, ordered by sequence position."""
    ca_atoms = []
    for res in chain.get_residues():
        if res.id[0] != ' ':
            continue
        if 'CA' in res:
            ca_atoms.append(res['CA'])
    return ca_atoms


def compute_rmsd(folded_pdb_path, backbone_pdb_path):
    """
    Compute Cα RMSD between folded binder (chain A) and backbone (chain A).
    Alignment is by sequence position, not residue number.
    Returns (rmsd, warning_message) where warning_message is None on success.
    """
    parser = PDBParser(QUIET=True)

    try:
        folded = parser.get_structure('folded', folded_pdb_path)
        backbone = parser.get_structure('backbone', backbone_pdb_path)
    except Exception as e:
        return None, f'Could not parse PDB: {e}'

    try:
        folded_ca = get_ca_atoms_by_position(folded[0]['A'])
        backbone_ca = get_ca_atoms_by_position(backbone[0]['A'])
    except KeyError as e:
        return None, f'Chain A not found: {e}'

    len_diff = abs(len(folded_ca) - len(backbone_ca))
    if len_diff > 5:
        return None, (f'Cα length mismatch too large: folded={len(folded_ca)}, '
                      f'backbone={len(backbone_ca)}, diff={len_diff} > 5')

    min_len = min(len(folded_ca), len(backbone_ca))
    if min_len == 0:
        return None, 'No Cα atoms found in one or both structures'

    folded_ca = folded_ca[:min_len]
    backbone_ca = backbone_ca[:min_len]

    sup = Superimposer()
    sup.set_atoms(backbone_ca, folded_ca)
    return float(sup.rms), None


# Extract receptor sequence once. A receptor in several segments is trimmed to
# the one carrying the epitope: ColabFold has no way to express a gap inside a
# chain, so the alternative is folding two lobes as though they were bonded.
receptor_bounds = None
segment_records = read_segment_records('./inputs/validated_target.pdb')
if len(segment_records) > 1:
    first, last, count = choose_receptor_segment(
        segment_records, read_hotspot_numbers('./inputs/validated_target.pdb'))
    receptor_bounds = (first, last)
    print(f'Receptor spans {len(segment_records)} segments; folding A{first}-{last} '
          f'({count} residues) — the segment carrying the epitope', flush=True)

receptor_seq = get_receptor_sequence('./inputs/validated_target.pdb', receptor_bounds)
if not receptor_seq:
    print('ERROR: Could not extract receptor sequence from validated_target.pdb', flush=True)
    sys.exit(1)
print(f'Receptor sequence length: {len(receptor_seq)}', flush=True)

# Gather FASTA files sorted by design ID
fasta_files = sorted(
    glob.glob('./inputs/sequences/design_*.fasta'),
    key=lambda p: os.path.basename(p)
)
if not fasta_files:
    print('ERROR: No FASTA files found in ./inputs/sequences/', flush=True)
    sys.exit(1)

print(f'Found {len(fasta_files)} FASTA file(s)', flush=True)

filter_report = []

for fasta_path in fasta_files:
    backbone_id = os.path.splitext(os.path.basename(fasta_path))[0]  # e.g. design_001
    # Map padded design ID back to original RFdiffusion unpadded index for RMSD
    backbone_index = int(backbone_id.split('_')[1])
    backbone_pdb = f'./inputs/backbones/design_{backbone_index}.pdb'

    # Parse all sequence records from the FASTA
    records = []
    with open(fasta_path) as f:
        header, seq_lines = None, []
        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq_lines)))
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        if header is not None:
            records.append((header, ''.join(seq_lines)))

    print(f'{backbone_id}: {len(records)} sequence record(s)', flush=True)

    for seq_idx, (header, full_seq) in enumerate(records):
        if seq_idx == 0:
            continue  # skip ProteinMPNN native/input sequence (poly-Gly for RFdiffusion backbones)
        design_seq_id = f'{backbone_id}_seq_{seq_idx}'

        # ProteinMPNN multi-chain sequences are /-separated; take the first part (chain A = binder)
        binder_seq = full_seq.split('/')[0].replace('-', '')
        if not binder_seq:
            print(f'  WARNING: Empty binder sequence for {design_seq_id} — skipping', flush=True)
            continue

        binder_length = len(binder_seq)
        merged_seq = f'{binder_seq}:{receptor_seq}'

        if ':' not in merged_seq or not merged_seq.split(':')[0] or not merged_seq.split(':')[1]:
            print(f'  ERROR: Malformed merged FASTA for {design_seq_id}: "{merged_seq[:60]}..."', flush=True)
            sys.exit(1)

        with tempfile.TemporaryDirectory() as tmpdir:
            cf_input = os.path.join(tmpdir, 'cf_input')
            cf_output = os.path.join(tmpdir, 'cf_output')
            os.makedirs(cf_input)
            os.makedirs(cf_output)

            fasta_out = os.path.join(cf_input, f'{design_seq_id}.fasta')
            with open(fasta_out, 'w') as f:
                f.write(f'>{design_seq_id}\n{merged_seq}\n')

            cmd = [
                'colabfold_batch',
                cf_input,
                cf_output,
                '--num-recycle', str(num_recycle),
                '--msa-mode', msa_mode,
                '--model-type', 'alphafold2_multimer_v3',
            ]
            print(f'  Folding {design_seq_id}...', flush=True)
            subprocess.run(cmd, check=True)

            # Find rank 001 score JSON — all scores including PAE are in this single file
            score_files = glob.glob(os.path.join(cf_output, '*_scores_rank_001_*.json'))
            if not score_files:
                print(f'  WARNING: No rank-001 score JSON found for {design_seq_id}', flush=True)
                filter_report.append({
                    'design_id': design_seq_id,
                    'iptm': None, 'plddt_binder': None, 'plddt_receptor': None,
                    'bb_rmsd': None, 'pae_interaction': None,
                    'pass': False, 'fail_reasons': ['no_score_file'],
                })
                continue

            with open(score_files[0]) as sf:
                scores = json.load(sf)

            iptm = float(scores.get('iptm', 0.0))

            # pLDDT binder: first binder_length elements of flat plddt array
            plddt_all = scores.get('plddt', [])
            plddt_binder = float(np.mean(plddt_all[:binder_length])) if plddt_all else 0.0
            # The receptor is the remainder of the same array. Diagnostic only, never
            # filtered on: without it a failing funnel cannot distinguish a bad design
            # from a receptor that was never modelled.
            plddt_receptor = (float(np.mean(plddt_all[binder_length:]))
                              if len(plddt_all) > binder_length else None)

            # PAE interaction: mean of the off-diagonal cross-chain block
            # Block rows 0..binder_length-1, cols binder_length..end (and its transpose)
            pae_matrix = scores.get('pae', [])
            if pae_matrix:
                pae_np = np.array(pae_matrix)
                n = pae_np.shape[0]
                b = binder_length
                if b < n:
                    block_br = pae_np[:b, b:]   # binder rows, receptor cols
                    block_rb = pae_np[b:, :b]   # receptor rows, binder cols
                    pae_interaction = float(np.mean(np.concatenate([block_br.flatten(), block_rb.flatten()])))
                else:
                    pae_interaction = float(np.mean(pae_np))
            else:
                pae_interaction = float('inf')

            # Find rank 001 unrelaxed PDB
            pdb_files = glob.glob(os.path.join(cf_output, '*_unrelaxed_rank_001*.pdb'))
            if not pdb_files:
                pdb_files = glob.glob(os.path.join(cf_output, '*rank_001*.pdb'))

            folded_pdb_path = pdb_files[0] if pdb_files else None

            # Compute Cα RMSD vs original backbone
            bb_rmsd = None
            rmsd_warning = None
            if folded_pdb_path and os.path.exists(backbone_pdb):
                bb_rmsd, rmsd_warning = compute_rmsd(folded_pdb_path, backbone_pdb)
                if rmsd_warning:
                    print(f'  RMSD warning for {design_seq_id}: {rmsd_warning}', flush=True)
            elif not os.path.exists(backbone_pdb):
                print(f'  WARNING: Backbone PDB not found: {backbone_pdb}', flush=True)

            # Apply filters — collect all failing reasons
            fail_reasons = []
            if iptm < min_iptm:
                fail_reasons.append('iptm')
            if plddt_binder < min_plddt_binder:
                fail_reasons.append('plddt_binder')
            if bb_rmsd is None:
                fail_reasons.append('bb_rmsd_unavailable')
            elif bb_rmsd > max_bb_rmsd:
                fail_reasons.append('bb_rmsd')
            if pae_interaction > max_pae_interaction:
                fail_reasons.append('pae_interaction')

            passed = len(fail_reasons) == 0

            record = {
                'design_id': design_seq_id,
                'iptm': round(iptm, 4),
                'plddt_binder': round(plddt_binder, 3),
                'plddt_receptor': round(plddt_receptor, 3) if plddt_receptor is not None else None,
                'bb_rmsd': round(bb_rmsd, 3) if bb_rmsd is not None else None,
                'pae_interaction': round(pae_interaction, 3) if pae_interaction != float('inf') else None,
                'pass': passed,
                'fail_reasons': fail_reasons,
            }
            filter_report.append(record)

            status = 'PASS' if passed else f'FAIL ({",".join(fail_reasons)})'
            rmsd_str = f'{bb_rmsd:.2f}' if bb_rmsd is not None else 'N/A'
            receptor_str = f'{plddt_receptor:.1f}' if plddt_receptor is not None else 'N/A'
            print(f'  {design_seq_id}: iPTM={iptm:.3f} pLDDT={plddt_binder:.1f} '
                  f'pLDDT_receptor={receptor_str} RMSD={rmsd_str} '
                  f'PAE={pae_interaction:.1f} → {status}', flush=True)

            if passed and folded_pdb_path:
                dest = f'./outputs/folded/{design_seq_id}_rank001.pdb'
                shutil.copy(folded_pdb_path, dest)
                print(f'  Copied to {dest}', flush=True)

with open('./outputs/folded/filter_report.json', 'w') as f:
    json.dump(filter_report, f, indent=2)

n_pass = sum(1 for r in filter_report if r['pass'])
print(f'\nFilter summary: {n_pass}/{len(filter_report)} designs passed', flush=True)
print('filter_report.json written', flush=True)
