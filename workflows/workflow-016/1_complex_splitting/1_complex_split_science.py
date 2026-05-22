#!/usr/bin/env python3
import argparse
import json
import os
import shutil
import sys
from typing import Dict, List, Tuple

from Bio import PDB
from Bio.PDB import PDBIO, Select
from Bio.PDB.Polypeptide import is_aa


class ChainSelector(Select):
    def __init__(self, chain_ids: List[str]):
        self.chain_ids = set(chain_ids)

    def accept_chain(self, chain):
        return chain.get_id() in self.chain_ids


def parse_chain_list(chain_csv: str) -> List[str]:
    if not chain_csv:
        return []
    parts = [c.strip() for c in chain_csv.replace(";", ",").split(",")]
    return [c for c in parts if c]


def get_chain_ids(structure) -> List[str]:
    # Use the first model only for deterministic chain order.
    model = next(structure.get_models())
    return [chain.get_id() for chain in model.get_chains()]


def summarize_chain(chain) -> Dict[str, int]:
    residues = [res for res in chain.get_residues() if is_aa(res, standard=True)]
    residue_count = len(residues)
    atom_count = sum(1 for _ in chain.get_atoms())
    hetero_residue_count = sum(1 for res in chain.get_residues() if res.get_id()[0] != " ")
    return {
        "residues": residue_count,
        "atoms": atom_count,
        "hetero_residues": hetero_residue_count,
    }


def collect_chain_statistics(structure) -> Dict[str, Dict[str, int]]:
    model = next(structure.get_models())
    stats: Dict[str, Dict[str, int]] = {}
    for chain in model.get_chains():
        stats[chain.get_id()] = summarize_chain(chain)
    return stats


def detect_antibody_chains(chain_ids: List[str]) -> Tuple[List[str], str]:
    if "H" in chain_ids and "L" in chain_ids:
        return ["H", "L"], "heuristic_H_L"
    if "VH" in chain_ids and "VL" in chain_ids:
        return ["VH", "VL"], "heuristic_VH_VL"
    if len(chain_ids) >= 2:
        return chain_ids[:2], "heuristic_first_two_chains"
    raise ValueError("Cannot auto-detect antibody chains: fewer than two chains in structure.")


def build_chain_info(
    chain_ids: List[str],
    chain_stats: Dict[str, Dict[str, int]],
    antibody_set: set,
) -> Dict[str, Dict[str, int]]:
    chain_info: Dict[str, Dict[str, int]] = {}
    for cid in chain_ids:
        info = dict(chain_stats[cid])
        info["type"] = "antibody" if cid in antibody_set else "antigen"
        chain_info[cid] = info
    return chain_info


def list_group_chain_info(chain_info: Dict[str, Dict[str, int]], chain_ids: List[str]) -> List[Dict[str, int]]:
    return [{"id": cid, **chain_info[cid]} for cid in chain_ids]


def sum_metric(chain_info: Dict[str, Dict[str, int]], chain_ids: List[str], metric: str) -> int:
    return sum(chain_info[cid][metric] for cid in chain_ids)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Split antibody-antigen complex PDB into antibody and antigen structures."
    )
    parser.add_argument("--input_pdb", required=True, help="Input co-crystal complex PDB path")
    parser.add_argument(
        "--antibody_chains",
        default="",
        help="Comma-separated antibody chain IDs (example: H,L). If omitted, auto-detection is used.",
    )
    parser.add_argument("--output_dir", default="1_complex_split_outputs", help="Output directory")
    args = parser.parse_args()

    if not os.path.isfile(args.input_pdb):
        raise FileNotFoundError(f"Input PDB not found: {args.input_pdb}")

    os.makedirs(args.output_dir, exist_ok=True)

    pdb_parser = PDB.PDBParser(QUIET=True)
    structure = pdb_parser.get_structure("complex", args.input_pdb)

    all_chain_ids = get_chain_ids(structure)
    if len(all_chain_ids) < 2:
        raise ValueError("Complex must contain at least two chains for antibody/antigen splitting.")

    user_antibody = parse_chain_list(args.antibody_chains)
    if user_antibody:
        missing = [cid for cid in user_antibody if cid not in all_chain_ids]
        if missing:
            raise ValueError(
                f"User-specified antibody chain(s) not found in structure: {', '.join(missing)}. "
                f"Available chains: {', '.join(all_chain_ids)}"
            )
        antibody_chains = user_antibody
        detection_method = "user_specified"
    else:
        antibody_chains, detection_method = detect_antibody_chains(all_chain_ids)

    antigen_chains = [cid for cid in all_chain_ids if cid not in set(antibody_chains)]
    if not antigen_chains:
        raise ValueError("Antigen chain set is empty after antibody assignment.")

    io = PDBIO()
    io.set_structure(structure)

    antibody_pdb = os.path.join(args.output_dir, "antibody.pdb")
    antigen_pdb = os.path.join(args.output_dir, "antigen.pdb")
    original_complex_pdb = os.path.join(args.output_dir, "original_complex.pdb")

    io.save(antibody_pdb, ChainSelector(antibody_chains))
    io.save(antigen_pdb, ChainSelector(antigen_chains))
    shutil.copy2(args.input_pdb, original_complex_pdb)

    chain_stats = collect_chain_statistics(structure)
    chain_info = build_chain_info(all_chain_ids, chain_stats, set(antibody_chains))

    chain_info_path = os.path.join(args.output_dir, "chain_info.json")
    chain_info_payload = {
        "input_file": os.path.abspath(args.input_pdb),
        "detection_method": detection_method,
        "all_chain_ids": all_chain_ids,
        "chain_info": chain_info,
    }
    with open(chain_info_path, "w", encoding="utf-8") as f:
        json.dump(chain_info_payload, f, indent=2)

    data = {
        "status": "success",
        "node": "complex_splitting",
        "input_file": os.path.abspath(args.input_pdb),
        "detection_method": detection_method,
        "antibody": {
            "chains": antibody_chains,
            "chain_count": len(antibody_chains),
            "residue_count": sum_metric(chain_info, antibody_chains, "residues"),
            "atom_count": sum_metric(chain_info, antibody_chains, "atoms"),
            "chain_info": list_group_chain_info(chain_info, antibody_chains),
            "pdb": os.path.abspath(antibody_pdb),
        },
        "antigen": {
            "chains": antigen_chains,
            "chain_count": len(antigen_chains),
            "residue_count": sum_metric(chain_info, antigen_chains, "residues"),
            "atom_count": sum_metric(chain_info, antigen_chains, "atoms"),
            "chain_info": list_group_chain_info(chain_info, antigen_chains),
            "pdb": os.path.abspath(antigen_pdb),
        },
        "outputs": {
            "antibody_pdb": os.path.abspath(antibody_pdb),
            "antigen_pdb": os.path.abspath(antigen_pdb),
            "original_complex": os.path.abspath(original_complex_pdb),
            "chain_info_json": os.path.abspath(chain_info_path),
            "data_json": os.path.abspath(os.path.join(args.output_dir, "data.json")),
        },
        "viewer_defaults": {
            "engine": "icn3d",
            "default_file": "antibody.pdb",
            "available_files": ["antibody.pdb", "antigen.pdb", "original_complex.pdb"],
            "width": "100%",
            "height": "560px",
            "style": "cartoon",
            "color_scheme": "chain",
        },
    }

    data_json_path = os.path.join(args.output_dir, "data.json")
    with open(data_json_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

    print(json.dumps(data, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())