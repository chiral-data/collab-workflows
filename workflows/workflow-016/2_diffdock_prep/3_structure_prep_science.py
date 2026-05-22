#!/usr/bin/env python3
import argparse
import json
import os
import sys
from typing import Dict, List, Tuple

from Bio import PDB
from Bio.PDB.Polypeptide import is_aa


def classify_residue(residue) -> str:
    het_flag = residue.get_id()[0]
    resname = residue.get_resname().strip().upper()

    if het_flag == "W" or resname in {"HOH", "WAT", "H2O"}:
        return "water"
    if het_flag != " ":
        return "hetero"
    if is_aa(residue, standard=True):
        return "standard_aa"
    return "nonstandard_aa"


def process_structure(input_pdb: str, output_pdb: str, role: str) -> Dict[str, object]:
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure(role, input_pdb)
    model = next(structure.get_models())

    chain_summaries: List[Dict[str, int]] = []
    total_input_residues = 0
    total_kept_residues = 0
    total_input_atoms = 0
    total_kept_atoms = 0
    removed_water = 0
    removed_hetero = 0
    removed_nonstandard_aa = 0

    for chain in model.get_chains():
        residues = list(chain.get_residues())
        input_residue_count = len(residues)
        input_atom_count = sum(1 for _ in chain.get_atoms())

        total_input_residues += input_residue_count
        total_input_atoms += input_atom_count

        keep_residues = []
        chain_removed_water = 0
        chain_removed_hetero = 0
        chain_removed_nonstandard = 0

        for residue in residues:
            kind = classify_residue(residue)
            if kind == "standard_aa":
                keep_residues.append(residue)
            elif kind == "water":
                chain_removed_water += 1
            elif kind == "hetero":
                chain_removed_hetero += 1
            else:
                chain_removed_nonstandard += 1

        for residue in residues:
            if residue not in keep_residues:
                chain.detach_child(residue.get_id())

        new_index = 1
        for residue in keep_residues:
            residue.id = (" ", new_index, " ")
            new_index += 1

        kept_atom_count = sum(1 for _ in chain.get_atoms())
        kept_residue_count = len(keep_residues)

        total_kept_residues += kept_residue_count
        total_kept_atoms += kept_atom_count
        removed_water += chain_removed_water
        removed_hetero += chain_removed_hetero
        removed_nonstandard_aa += chain_removed_nonstandard

        chain_summaries.append(
            {
                "chain_id": chain.get_id(),
                "input_residues": input_residue_count,
                "kept_residues": kept_residue_count,
                "removed_water_residues": chain_removed_water,
                "removed_hetero_residues": chain_removed_hetero,
                "removed_nonstandard_aa_residues": chain_removed_nonstandard,
                "input_atoms": input_atom_count,
                "kept_atoms": kept_atom_count,
            }
        )

    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

    return {
        "role": role,
        "input": os.path.abspath(input_pdb),
        "processed": os.path.abspath(output_pdb),
        "chain_count": len(chain_summaries),
        "totals": {
            "input_residues": total_input_residues,
            "kept_residues": total_kept_residues,
            "removed_residues": total_input_residues - total_kept_residues,
            "input_atoms": total_input_atoms,
            "kept_atoms": total_kept_atoms,
            "removed_water_residues": removed_water,
            "removed_hetero_residues": removed_hetero,
            "removed_nonstandard_aa_residues": removed_nonstandard_aa,
        },
        "chain_summaries": chain_summaries,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Clean and renumber antibody/antigen PDBs for DiffDock-PP input preparation."
    )
    parser.add_argument("--antibody_pdb", required=True, help="Path to antibody.pdb from Node 1")
    parser.add_argument("--antigen_pdb", required=True, help="Path to antigen.pdb from Node 1")
    parser.add_argument("--output_dir", default="2_structure_prep_outputs", help="Output directory")
    args = parser.parse_args()

    if not os.path.isfile(args.antibody_pdb):
        raise FileNotFoundError(f"Antibody PDB not found: {args.antibody_pdb}")
    if not os.path.isfile(args.antigen_pdb):
        raise FileNotFoundError(f"Antigen PDB not found: {args.antigen_pdb}")

    os.makedirs(args.output_dir, exist_ok=True)

    processed_antibody = os.path.join(args.output_dir, "processed_antibody.pdb")
    processed_antigen = os.path.join(args.output_dir, "processed_antigen.pdb")

    antibody = process_structure(args.antibody_pdb, processed_antibody, role="receptor")
    antigen = process_structure(args.antigen_pdb, processed_antigen, role="ligand")

    data = {
        "status": "success",
        "node": "structure_preparation",
        "description": "Water/hetero removed, standard amino acids kept, residues renumbered per chain from 1.",
        "geometric_modification": False,
        "antibody": antibody,
        "antigen": antigen,
        "outputs": {
            "processed_antibody": os.path.abspath(processed_antibody),
            "processed_antigen": os.path.abspath(processed_antigen),
            "data_json": os.path.abspath(os.path.join(args.output_dir, "data.json")),
        },
    }

    with open(os.path.join(args.output_dir, "data.json"), "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

    print(json.dumps(data, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())