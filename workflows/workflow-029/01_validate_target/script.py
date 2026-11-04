#!/usr/bin/env python3
import os
import sys

from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, Select

# An ideal trans peptide bond spans 3.80 Å Cα–Cα; 4.5 Å leaves headroom for
# strained geometry without admitting a real break.
CA_BOND_MAX = 4.5


def detect_segments(residues, verbose=True, validate=True):
    """Contiguous runs of the chain, and the breaks between them.

    RFdiffusion indexes CA-bearing residues only and ignores insertion codes, so
    a residue it cannot address must not fall inside a contig range. validate is
    off for the pass taken before renumbering, which is about to replace the very
    numbering the checks reject.
    """
    usable = []
    skipped_icode = 0
    skipped_no_ca = 0
    for res in residues:
        if res.id[2] != ' ':
            skipped_icode += 1
        elif 'CA' not in res:
            skipped_no_ca += 1
        else:
            usable.append(res)

    if verbose and skipped_icode:
        print(f'WARNING: {skipped_icode} residue(s) carry an insertion code and are left out of '
              'the contig — RFdiffusion collapses them onto the residue they follow', flush=True)
    if verbose and skipped_no_ca:
        print(f'WARNING: {skipped_no_ca} residue(s) have no CA atom and are left out of the '
              'contig', flush=True)

    if not usable:
        print('ERROR: No CA-bearing residues remain in the receptor chain', flush=True)
        sys.exit(1)

    numbers = [r.id[1] for r in usable]
    if validate and numbers[0] < 1:
        print(f'ERROR: Residue numbering starts at {numbers[0]}. A contig cannot express a '
              'non-positive residue number — set renumber_residues=true for this target.',
              flush=True)
        sys.exit(1)
    if any(b <= a for a, b in zip(numbers, numbers[1:])):
        if verbose:
            print('WARNING: Residue numbering is not ascending in file order — sorting', flush=True)
        usable.sort(key=lambda r: r.id[1])
        numbers = [r.id[1] for r in usable]
        if validate and any(b == a for a, b in zip(numbers, numbers[1:])):
            print('ERROR: Duplicate residue numbers in the receptor chain. A contig cannot name '
                  'them unambiguously — set renumber_residues=true for this target.', flush=True)
            sys.exit(1)

    segments = [[usable[0]]]
    breaks = []
    for prev, nxt in zip(usable, usable[1:]):
        missing = nxt.id[1] - prev.id[1] - 1
        distance = nxt['CA'] - prev['CA']
        if missing or distance > CA_BOND_MAX:
            if missing and distance > CA_BOND_MAX:
                kind = 'BOTH'
            elif missing:
                kind = 'NUMBERING'
            else:
                kind = 'DISTANCE'
            breaks.append((prev.id[1], nxt.id[1], missing, distance, kind))
        if missing:
            segments.append([nxt])
        else:
            # A break with no gap in numbering is inexpressible: the index jump a
            # contig encodes is next - prev - 1, which is zero here. Reported only.
            segments[-1].append(nxt)

    return [(s[0].id[1], s[-1].id[1], len(s)) for s in segments], breaks


def main():
    os.makedirs('./outputs', exist_ok=True)

    receptor_chain = os.environ.get('PARAM_RECEPTOR_CHAIN', 'A')
    hotspot_residues = os.environ.get('PARAM_HOTSPOT_RESIDUES', '')
    strip_heteroatoms = os.environ.get('PARAM_STRIP_HETEROATOMS', 'true').lower() == 'true'
    renumber_residues = os.environ.get('PARAM_RENUMBER_RESIDUES', 'false').lower() == 'true'

    print(f'Input receptor chain: {receptor_chain}', flush=True)
    print(f'Strip heteroatoms: {strip_heteroatoms}', flush=True)
    print(f'Renumber residues: {renumber_residues}', flush=True)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('target', './inputs/target.pdb')

    model = structure[0]
    chain_ids = [c.id for c in model.get_chains()]
    print(f'Chains found in input: {chain_ids}', flush=True)

    if receptor_chain not in chain_ids:
        print(f'ERROR: Chain {receptor_chain} not found in input PDB. Available chains: {chain_ids}', flush=True)
        sys.exit(1)

    chain = model[receptor_chain]

    if strip_heteroatoms:
        hetero_residues = [r for r in chain.get_residues() if r.id[0] != ' ']
        for res in hetero_residues:
            chain.detach_child(res.id)
        print(f'Stripped {len(hetero_residues)} heteroatom residue(s)', flush=True)

    atom_residues = [r for r in chain.get_residues() if r.id[0] == ' ']
    residue_count = len(atom_residues)

    # Resolve hotspots against the original numbering the user supplied. With
    # renumbering off this is an identity mapping, which is what RFdiffusion wants —
    # it matches hotspots against the input numbering. Node 02 reads the result from
    # the REMARK records below either way.
    hotspots_requested = bool(hotspot_residues.strip())
    hotspot_mapping = []

    if hotspots_requested:
        hotspot_list = [h.strip() for h in hotspot_residues.split(',') if h.strip()]
        # Position in this list is the renumbered ID (1-based), since renumbering is sequential
        original_ids = [res.id[1] for res in atom_residues]
        for hotspot in hotspot_list:
            chain_prefix = hotspot[:len(hotspot) - len(hotspot.lstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))]
            try:
                res_num = int(hotspot[len(chain_prefix):])
            except ValueError:
                print(f'WARNING: Could not parse hotspot residue ID "{hotspot}" — skipping', flush=True)
                continue
            # A bare number means the receptor chain; a letter naming another chain is
            # a request this node cannot honour, since only receptor_chain is kept
            if chain_prefix and chain_prefix != receptor_chain:
                print(f'WARNING: Hotspot residue {hotspot} names chain {chain_prefix}, but the '
                      f'receptor chain is {receptor_chain} — dropping it from the hotspot set', flush=True)
                continue
            if res_num not in original_ids:
                print(f'WARNING: Hotspot residue {hotspot} (number {res_num}) is not present in the '
                      f'chain after processing — dropping it from the hotspot set', flush=True)
                continue
            mapped_num = original_ids.index(res_num) + 1 if renumber_residues else res_num
            hotspot_mapping.append((hotspot, f'A{mapped_num}'))
            print(f'Hotspot residue {hotspot} → A{mapped_num}', flush=True)

        if not hotspot_mapping:
            print('WARNING: None of the requested hotspot residues could be mapped — downstream '
                  'design will run without hotspots', flush=True)

    if renumber_residues:
        _, fused_breaks = detect_segments(atom_residues, verbose=False, validate=False)
        if fused_breaks:
            fused = ', '.join(f'{p}|{n} ({m} missing, {d:.1f} Å)'
                              for p, n, m, d, _ in fused_breaks)
            print(f'WARNING: renumbering fuses {len(fused_breaks)} chain break(s) — RFdiffusion will '
                  f'be given these as bonded residues: {fused}', flush=True)
        for res in atom_residues:
            chain.detach_child(res.id)
        for i, res in enumerate(atom_residues, start=1):
            res.id = (' ', i, ' ')
            chain.add(res)
        print(f'Renumbered {residue_count} residues starting from 1', flush=True)

    # The remarks describe the file as written, so this runs after any renumbering
    segments, breaks = detect_segments(atom_residues)

    # Always write output as chain A regardless of input chain ID
    chain.id = 'A'

    # Build a minimal structure containing only the cleaned chain
    new_structure = PDB.Structure.Structure('validated')
    new_model = PDB.Model.Model(0)
    new_model.add(chain)
    new_structure.add(new_model)

    io = PDBIO()
    io.set_structure(new_structure)
    io.save('./outputs/validated_target.pdb')

    # PDBIO writes no remarks, so they are prepended to the saved file. Each summary
    # line is written even when it carries nothing new, so node 02 can tell "node 01
    # handled this" from "older input".
    remarks = []
    if hotspots_requested:
        remarks += [f'REMARK 999 HOTSPOT {original} {mapped}\n'
                    for original, mapped in hotspot_mapping]
        remarks.append('REMARK 999 HOTSPOTS %s\n' % ','.join(m for _, m in hotspot_mapping))

    remarks.append('REMARK 999 NUMBERING %s\n' % ('RENUMBERED' if renumber_residues else 'ORIGINAL'))
    for index, (first, last, count) in enumerate(segments, start=1):
        remarks.append(f'REMARK 999 SEGMENT {index} A {first} {last} {count}\n')
    remarks.append('REMARK 999 SEGMENTS %s\n' % '/'.join(f'A{f}-{l}' for f, l, _ in segments))
    for prev, nxt, missing, distance, kind in breaks:
        remarks.append(f'REMARK 999 BREAK {prev} {nxt} {missing} {distance:.2f} {kind}\n')

    with open('./outputs/validated_target.pdb') as f:
        body = f.read()
    with open('./outputs/validated_target.pdb', 'w') as f:
        f.write(''.join(remarks) + body)

    segment_residues = sum(count for _, _, count in segments)
    summary = f'Receptor segments: {len(segments)} ({segment_residues} residues'
    if breaks:
        summary += ', %d chain break(s): %s' % (
            len(breaks), ', '.join(f'{p}|{n} {d:.1f} Å' for p, n, _, d, _ in breaks))
    print(summary + ')', flush=True)
    print(f'Chain normalization: {receptor_chain} → A', flush=True)
    print(f'Residue count: {residue_count}', flush=True)
    print(f'Output written to ./outputs/validated_target.pdb', flush=True)


if __name__ == '__main__':
    main()
