#!/usr/bin/env python3
import os
import sys

from Bio import PDB
from Bio.PDB import PDBParser, PDBIO, Select

os.makedirs('./outputs', exist_ok=True)

receptor_chain = os.environ.get('PARAM_RECEPTOR_CHAIN', 'A')
hotspot_residues = os.environ.get('PARAM_HOTSPOT_RESIDUES', '')
strip_heteroatoms = os.environ.get('PARAM_STRIP_HETEROATOMS', 'true').lower() == 'true'
renumber_residues = os.environ.get('PARAM_RENUMBER_RESIDUES', 'true').lower() == 'true'

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

# Resolve hotspots against the original numbering the user supplied, then carry
# the result through this node's renumbering: node 02 reads the mapped IDs from
# the REMARK records below, because by then the original numbers are gone.
hotspots_requested = bool(hotspot_residues.strip())
hotspot_mapping = []

if hotspots_requested:
    hotspot_list = [h.strip() for h in hotspot_residues.split(',') if h.strip()]
    # Position in this list is the renumbered ID (1-based), since renumbering is sequential
    original_ids = [res.id[1] for res in atom_residues]
    for hotspot in hotspot_list:
        try:
            res_num = int(hotspot.lstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))
        except ValueError:
            print(f'WARNING: Could not parse hotspot residue ID "{hotspot}" — skipping', flush=True)
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
    for res in atom_residues:
        chain.detach_child(res.id)
    for i, res in enumerate(atom_residues, start=1):
        res.id = (' ', i, ' ')
        chain.add(res)
    print(f'Renumbered {residue_count} residues starting from 1', flush=True)

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

# PDBIO writes no remarks, so the mapping is prepended to the saved file. The
# summary line is written whenever hotspots were requested, even when nothing
# mapped, so node 02 can tell "node 01 handled this" from "older input".
if hotspots_requested:
    remarks = [f'REMARK 999 HOTSPOT {original} {mapped}\n'
               for original, mapped in hotspot_mapping]
    remarks.append('REMARK 999 HOTSPOTS %s\n' % ','.join(m for _, m in hotspot_mapping))
    with open('./outputs/validated_target.pdb') as f:
        body = f.read()
    with open('./outputs/validated_target.pdb', 'w') as f:
        f.write(''.join(remarks) + body)

print(f'Chain normalization: {receptor_chain} → A', flush=True)
print(f'Residue count: {residue_count}', flush=True)
print(f'Output written to ./outputs/validated_target.pdb', flush=True)
