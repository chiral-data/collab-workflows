#!/usr/bin/env python3
import argparse
import json
import os
import stat
import sys
from typing import Dict, List, Tuple

from Bio import PDB
from Bio import SeqIO
from Bio.Data.IUPACData import protein_letters_3to1_extended
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


AA_3TO1 = {k.upper(): v.upper() for k, v in protein_letters_3to1_extended.items()}
AA_3TO1.update({"MSE": "M"})


def residue_to_one_letter(residue) -> str:
    resname = residue.get_resname().strip().upper()
    return AA_3TO1.get(resname, "X")


def extract_sequences(pdb_file: str) -> Tuple[Dict[str, str], List[Dict[str, object]]]:
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)
    model = next(structure.get_models())

    sequences: Dict[str, str] = {}
    chain_info: List[Dict[str, object]] = []

    for chain in model.get_chains():
        chain_id = chain.get_id()
        residues = [res for res in chain.get_residues() if res.get_id()[0] == " "]
        one_letter = [residue_to_one_letter(res) for res in residues]
        seq = "".join(one_letter)

        if not seq:
            continue

        unknown_count = sum(1 for aa in one_letter if aa == "X")
        sequences[chain_id] = seq
        chain_info.append(
            {
                "chain_id": chain_id,
                "length": len(seq),
                "unknown_residue_count": unknown_count,
                "known_residue_fraction": round((len(seq) - unknown_count) / len(seq), 4),
                "sequence": seq,
            }
        )

    chain_info.sort(key=lambda x: x["chain_id"])
    return sequences, chain_info


def write_fasta(entity_name: str, sequences: Dict[str, str], output_fasta: str) -> None:
    records = [
        SeqRecord(Seq(seq), id=f"{entity_name}_{chain_id}", description=f"chain={chain_id}")
        for chain_id, seq in sorted(sequences.items())
    ]
    SeqIO.write(records, output_fasta, "fasta")


def write_placeholder_pt(path: str, entity_name: str, chain_info: List[Dict[str, object]]) -> None:
    payload = {
        "placeholder": True,
        "entity": entity_name,
        "format": "pt_placeholder_text",
        "note": "Replace with real ESM-2 embeddings using generate_embeddings.py or esm-extract.",
        "chain_lengths": {row["chain_id"]: row["length"] for row in chain_info},
    }
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


def write_embedding_helper_script(output_dir: str) -> str:
    script_path = os.path.join(output_dir, "generate_embeddings.py")
    content = """#!/usr/bin/env python3
import argparse
import os
import shutil
import subprocess
import sys


def run_esm_extract(model, fasta, out_path, repr_layers):
    cmd = [
        "esm-extract",
        model,
        fasta,
        out_path,
        "--repr_layers",
        str(repr_layers),
        "--include",
        "mean",
        "per_tok",
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(description="Generate real ESM-2 embeddings from FASTA files.")
    parser.add_argument("--model", default="esm2_t36_3B_UR50D")
    parser.add_argument("--repr_layers", type=int, default=36)
    parser.add_argument("--antibody_fasta", default="antibody.fasta")
    parser.add_argument("--antigen_fasta", default="antigen.fasta")
    parser.add_argument("--antibody_out", default="antibody_features.pt")
    parser.add_argument("--antigen_out", default="antigen_features.pt")
    args = parser.parse_args()

    if shutil.which("esm-extract") is None:
        print("ERROR: esm-extract is not available in PATH.")
        print("Install with: pip install fair-esm")
        sys.exit(1)

    for fasta in (args.antibody_fasta, args.antigen_fasta):
        if not os.path.isfile(fasta):
            print(f"ERROR: FASTA file not found: {fasta}")
            sys.exit(1)

    run_esm_extract(args.model, args.antibody_fasta, args.antibody_out, args.repr_layers)
    run_esm_extract(args.model, args.antigen_fasta, args.antigen_out, args.repr_layers)
    print("Embedding generation complete.")


if __name__ == "__main__":
    main()
"""
    with open(script_path, "w", encoding="utf-8") as f:
        f.write(content)
    current_mode = os.stat(script_path).st_mode
    os.chmod(script_path, current_mode | stat.S_IXUSR)
    return script_path


def summarize_entity(name: str, input_pdb: str, fasta_path: str, feature_path: str, chain_info: List[Dict[str, object]]) -> Dict[str, object]:
    total_length = sum(row["length"] for row in chain_info)
    unknown_total = sum(row["unknown_residue_count"] for row in chain_info)
    return {
        "input_pdb": os.path.abspath(input_pdb),
        "role": "receptor" if name == "antibody" else "ligand",
        "chains": [row["chain_id"] for row in chain_info],
        "chain_count": len(chain_info),
        "total_length": total_length,
        "unknown_residue_count": unknown_total,
        "fasta": os.path.abspath(fasta_path),
        "features_placeholder": os.path.abspath(feature_path),
        "chain_info": chain_info,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Extract sequences and prepare ESM-2 embedding inputs.")
    parser.add_argument("--antibody_pdb", required=True)
    parser.add_argument("--antigen_pdb", required=True)
    parser.add_argument("--output_dir", default="3_feature_extract_outputs")
    args = parser.parse_args()

    if not os.path.isfile(args.antibody_pdb):
        raise FileNotFoundError(f"Antibody PDB not found: {args.antibody_pdb}")
    if not os.path.isfile(args.antigen_pdb):
        raise FileNotFoundError(f"Antigen PDB not found: {args.antigen_pdb}")

    os.makedirs(args.output_dir, exist_ok=True)

    ab_seq, ab_chain_info = extract_sequences(args.antibody_pdb)
    ag_seq, ag_chain_info = extract_sequences(args.antigen_pdb)

    if not ab_seq:
        raise ValueError("No valid amino-acid residues found in processed antibody PDB.")
    if not ag_seq:
        raise ValueError("No valid amino-acid residues found in processed antigen PDB.")

    ab_fasta = os.path.join(args.output_dir, "antibody.fasta")
    ag_fasta = os.path.join(args.output_dir, "antigen.fasta")
    ab_feat = os.path.join(args.output_dir, "antibody_features.pt")
    ag_feat = os.path.join(args.output_dir, "antigen_features.pt")

    write_fasta("antibody", ab_seq, ab_fasta)
    write_fasta("antigen", ag_seq, ag_fasta)
    write_placeholder_pt(ab_feat, "antibody", ab_chain_info)
    write_placeholder_pt(ag_feat, "antigen", ag_chain_info)
    embed_script = write_embedding_helper_script(args.output_dir)

    antibody_summary = summarize_entity("antibody", args.antibody_pdb, ab_fasta, ab_feat, ab_chain_info)
    antigen_summary = summarize_entity("antigen", args.antigen_pdb, ag_fasta, ag_feat, ag_chain_info)

    sequence_info = {
        "antibody": {
            "chains": antibody_summary["chains"],
            "chain_lengths": {row["chain_id"]: row["length"] for row in ab_chain_info},
            "total_length": antibody_summary["total_length"],
        },
        "antigen": {
            "chains": antigen_summary["chains"],
            "chain_lengths": {row["chain_id"]: row["length"] for row in ag_chain_info},
            "total_length": antigen_summary["total_length"],
        },
    }

    sequence_info_path = os.path.join(args.output_dir, "sequence_info.json")
    with open(sequence_info_path, "w", encoding="utf-8") as f:
        json.dump(sequence_info, f, indent=2)

    data = {
        "status": "success",
        "node": "feature_extraction",
        "description": "Sequence extraction completed. ESM-2 feature files are placeholders.",
        "esm2": {
            "embeddings_generated": False,
            "placeholder_files": [os.path.abspath(ab_feat), os.path.abspath(ag_feat)],
            "generation_script": os.path.abspath(embed_script),
            "recommended_model": "esm2_t36_3B_UR50D",
            "recommended_repr_layers": [36],
        },
        "antibody": antibody_summary,
        "antigen": antigen_summary,
        "outputs": {
            "antibody_fasta": os.path.abspath(ab_fasta),
            "antigen_fasta": os.path.abspath(ag_fasta),
            "sequence_info_json": os.path.abspath(sequence_info_path),
            "data_json": os.path.abspath(os.path.join(args.output_dir, "data.json")),
        },
    }

    data_path = os.path.join(args.output_dir, "data.json")
    with open(data_path, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

    print(json.dumps(data, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())