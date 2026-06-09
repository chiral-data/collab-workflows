#!/usr/bin/env python3
import glob
import os
import shutil
import subprocess
import sys

from Bio import SeqIO


def check_gpu():
    result = subprocess.run(["nvidia-smi"], capture_output=True)
    if result.returncode == 0:
        return
    # Fallback: try torch
    try:
        import torch
        if torch.cuda.is_available():
            return
        print("ERROR: No GPU available (torch.cuda.is_available() returned False)", flush=True)
        sys.exit(1)
    except ImportError:
        pass
    print("ERROR: No GPU available (nvidia-smi failed, torch not importable)", flush=True)
    sys.exit(1)


def detect_mode(fasta_path):
    records = list(SeqIO.parse(fasta_path, "fasta"))
    is_multimer = len(records) > 1 or any(
        ":" in r.id or ":" in r.description for r in records
    )
    return "multimer" if is_multimer else "monomer", records


def main():
    dry_run = "--dry-run" in sys.argv

    os.makedirs("./outputs", exist_ok=True)

    fasta_path = "./inputs/validated_sequences.fasta"
    if not os.path.exists(fasta_path):
        print(f"ERROR: Input FASTA not found: {fasta_path}", flush=True)
        sys.exit(1)

    mode, _ = detect_mode(fasta_path)
    print(f"Detected mode: {mode}", flush=True)

    if not dry_run:
        check_gpu()

    # Read parameters
    num_models = int(os.environ.get("PARAM_NUM_MODELS", "5"))
    num_recycle = int(os.environ.get("PARAM_NUM_RECYCLE", "3"))
    msa_mode = os.environ.get("PARAM_MSA_MODE", "mmseqs2_uniref_env")
    pair_mode = os.environ.get("PARAM_PAIR_MODE", "auto")
    use_templates = os.environ.get("PARAM_USE_TEMPLATES", "false").lower() == "true"
    host_url = os.environ.get("PARAM_HOST_URL", "https://api.colabfold.com")

    # Auto-resolve pair_mode
    if pair_mode == "auto":
        pair_mode = "unpaired_paired" if mode == "multimer" else "unpaired"

    # colabfold_batch requires a directory as input
    os.makedirs("./cf_input", exist_ok=True)
    shutil.copy(fasta_path, "./cf_input/sequences.fasta")

    cmd = [
        "colabfold_batch",
        "./cf_input",
        "./cf_output",
        "--num-models", str(num_models),
        "--num-recycle", str(num_recycle),
        "--msa-mode", msa_mode,
        "--pair-mode", pair_mode,
        "--host-url", host_url,
    ]
    if use_templates:
        cmd.append("--templates")

    print(f"Resolved command: {' '.join(cmd)}", flush=True)

    if dry_run:
        print("DRY RUN: exiting without running prediction", flush=True)
        sys.exit(0)

    subprocess.run(cmd, check=True)

    # Copy prediction outputs flat into ./outputs/
    copied = 0
    for src in glob.glob("./cf_output/**/*", recursive=True):
        if os.path.isfile(src):
            dest = f"./outputs/{os.path.basename(src)}"
            shutil.copy2(src, dest)
            print(f"  Collected: {os.path.basename(src)}", flush=True)
            copied += 1

    print(f"Collected {copied} output file(s)", flush=True)

    if copied == 0:
        print("ERROR: No output files found from colabfold_batch", flush=True)
        print("Keeping ./cf_output/ for debugging", flush=True)
        shutil.rmtree("./cf_input", ignore_errors=True)
        sys.exit(1)

    # Pass validated_sequences.fasta through so node 03 can read chain lengths
    shutil.copy(fasta_path, "./outputs/validated_sequences.fasta")
    print("Pass-through: validated_sequences.fasta", flush=True)

    # Clean up intermediate dirs
    shutil.rmtree("./cf_input", ignore_errors=True)
    shutil.rmtree("./cf_output", ignore_errors=True)


if __name__ == "__main__":
    main()
