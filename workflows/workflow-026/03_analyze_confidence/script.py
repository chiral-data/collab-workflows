#!/usr/bin/env python3
import glob
import json
import os
import re
import sys

import numpy as np
from Bio import SeqIO


def get_chain_lengths(fasta_path):
    if not os.path.exists(fasta_path):
        return []
    records = list(SeqIO.parse(fasta_path, "fasta"))
    return [len(r.seq) for r in records]


def interface_pae_mean(pae_matrix, chain_lengths):
    """Mean PAE of off-diagonal (inter-chain) blocks."""
    pae = np.array(pae_matrix)
    n = pae.shape[0]
    intra_mask = np.zeros((n, n), dtype=bool)
    offset = 0
    for length in chain_lengths:
        intra_mask[offset:offset + length, offset:offset + length] = True
        offset += length
    off_diag = pae[~intra_mask]
    return float(np.mean(off_diag)) if off_diag.size > 0 else None


def main():
    os.makedirs("./outputs", exist_ok=True)

    top_n = int(os.environ.get("PARAM_TOP_N_MODELS", "5"))
    colabfold_version = os.environ.get("PARAM_COLABFOLD_VERSION", "1.6.1")
    job_id = os.environ.get("PARAM_JOB_ID", "")

    # Find score JSONs
    score_files = glob.glob("./inputs/predictions/*_scores_rank_*.json")
    if not score_files:
        print("ERROR: No score JSON files found in ./inputs/", flush=True)
        sys.exit(1)

    print(f"Found {len(score_files)} score file(s)", flush=True)

    # Derive job_id from first score filename if not set
    if not job_id:
        basename = os.path.basename(score_files[0])
        m = re.match(r"(.+?)_scores_rank_", basename)
        job_id = m.group(1) if m else "unknown"

    # Get chain lengths for interface PAE (multimer only)
    chain_lengths = get_chain_lengths("./inputs/validated_sequences.fasta")

    # Parse each score file and extract rank
    parsed = []
    for path in score_files:
        m = re.search(r"rank_(\d+)", os.path.basename(path))
        if not m:
            continue
        rank = int(m.group(1))
        with open(path) as f:
            data = json.load(f)
        parsed.append((rank, path, data))

    parsed.sort(key=lambda x: x[0])
    parsed = parsed[:top_n]

    # Detect mode
    is_multimer = any("iptm" in d for _, _, d in parsed)
    mode = "multimer" if is_multimer else "monomer"
    print(f"Mode: {mode}", flush=True)

    models = []
    for rank, path, data in parsed:
        model_name = os.path.basename(path).replace("_scores_", "_unrelaxed_")
        model_name = re.sub(r"\.json$", ".pdb", model_name)

        plddt_arr = data["plddt"]
        pae_mat = data["pae"]
        plddt_mean = float(np.mean(plddt_arr))
        pae_mean = float(np.mean(pae_mat))
        max_pae = float(data["max_pae"])
        ptm = float(data["ptm"])
        iptm = float(data["iptm"]) if "iptm" in data else None

        iface_pae = None
        if is_multimer and len(chain_lengths) > 1:
            iface_pae = interface_pae_mean(pae_mat, chain_lengths)

        models.append({
            "rank": rank,
            "model_name": model_name,
            "plddt_mean": plddt_mean,
            "plddt_per_residue": [float(v) for v in plddt_arr],
            "pae_mean": pae_mean,
            "max_pae": max_pae,
            "pae": [[float(v) for v in row] for row in pae_mat],
            "interface_pae_mean": iface_pae,
            "ptm": ptm,
            "iptm": iptm,
        })

        print(
            f"  rank {rank}: pLDDT_mean={plddt_mean:.1f}, PAE_mean={pae_mean:.2f}, "
            f"PTM={ptm:.3f}" + (f", iPTM={iptm:.3f}" if iptm is not None else ""),
            flush=True,
        )

    summary = {
        "colabfold_version": colabfold_version,
        "job_id": job_id,
        "mode": mode,
        "models": models,
    }

    out_path = "./outputs/confidence_summary.json"
    with open(out_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"Written: {out_path}", flush=True)


if __name__ == "__main__":
    main()
