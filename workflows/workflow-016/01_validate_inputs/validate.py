#!/usr/bin/env python3
"""Validate receptor and ligand PDB files for LightDock docking."""

import json
import os
import shutil
import sys


def parse_pdb(path):
    """Return (stats_dict, error_str). error_str is None on success."""
    if not os.path.exists(path):
        return None, f"File not found: {path}"

    atoms = []
    residues = set()
    chains = set()
    has_ca = False

    with open(path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                try:
                    atom_name = line[12:16].strip()
                    chain = line[21]
                    resseq = int(line[22:26].strip())
                    float(line[30:38])  # x — confirm parseable
                    atoms.append(atom_name)
                    residues.add((chain, resseq))
                    chains.add(chain)
                    if atom_name == "CA":
                        has_ca = True
                except (ValueError, IndexError):
                    pass

    if not atoms:
        return None, f"No ATOM/HETATM records found in {path}"
    if not has_ca:
        return None, f"No CA atoms found in {path} — may not be a protein structure"

    return {
        "num_atoms": len(atoms),
        "num_residues": len(residues),
        "num_chains": len(chains),
        "chains": sorted(chains),
    }, None


def main():
    receptor_path = os.environ.get("PARAM_RECEPTOR_FILE", "inputs/protein2_barnase.pdb")
    ligand_path = os.environ.get("PARAM_LIGAND_FILE", "inputs/protein1_barstar.pdb")

    print(f"Receptor: {receptor_path}", flush=True)
    print(f"Ligand:   {ligand_path}", flush=True)

    rec_info, rec_err = parse_pdb(receptor_path)
    if rec_err:
        print(f"ERROR (receptor): {rec_err}", flush=True)
        sys.exit(1)

    lig_info, lig_err = parse_pdb(ligand_path)
    if lig_err:
        print(f"ERROR (ligand): {lig_err}", flush=True)
        sys.exit(1)

    print(
        f"Receptor: {rec_info['num_atoms']} atoms, "
        f"{rec_info['num_residues']} residues, "
        f"chains={rec_info['chains']}",
        flush=True,
    )
    print(
        f"Ligand:   {lig_info['num_atoms']} atoms, "
        f"{lig_info['num_residues']} residues, "
        f"chains={lig_info['chains']}",
        flush=True,
    )

    shutil.copy(receptor_path, "receptor.pdb")
    shutil.copy(ligand_path, "ligand.pdb")
    print("Copied receptor.pdb and ligand.pdb.", flush=True)

    report = {
        "receptor": {"path": receptor_path, **rec_info},
        "ligand": {"path": ligand_path, **lig_info},
        "status": "ok",
    }
    with open("validation_report.json", "w") as f:
        json.dump(report, f, indent=2)

    print("Validation complete.", flush=True)


if __name__ == "__main__":
    main()
