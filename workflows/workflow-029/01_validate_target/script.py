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

# Validate hotspot residues BEFORE renumbering (user provides original PDB numbers)
if hotspot_residues.strip():
    hotspot_list = [h.strip() for h in hotspot_residues.split(',') if h.strip()]
    residue_ids = {res.id[1] for res in atom_residues}
    for hotspot in hotspot_list:
        try:
            res_num = int(hotspot.lstrip('ABCDEFGHIJKLMNOPQRSTUVWXYZ'))
        except ValueError:
            print(f'WARNING: Could not parse hotspot residue ID "{hotspot}" — skipping', flush=True)
            continue
        if res_num not in residue_ids:
            print(f'WARNING: Hotspot residue {hotspot} (number {res_num}) not found in chain after processing', flush=True)
        else:
            print(f'Hotspot residue {hotspot} validated OK', flush=True)

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

print(f'Chain normalization: {receptor_chain} → A', flush=True)
print(f'Residue count: {residue_count}', flush=True)
print(f'Output written to ./outputs/validated_target.pdb', flush=True)
