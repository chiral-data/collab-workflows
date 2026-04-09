#!/usr/bin/env python3
"""Download mDeepFRI model weights (GCN + CNN) and pass through validated FASTA.

Node 02 prepares the structural alignment environment: the GCN models perform
graph-based inference on contact maps derived from FoldComp structural alignments,
while CNN models handle sequence-only prediction for proteins without structural hits.
"""

import argparse
import json
import os
import shutil
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="Download mDeepFRI model weights")
    parser.add_argument("--model-version", default="1.1", choices=["1.0", "1.1"])
    args = parser.parse_args()

    os.makedirs("./outputs", exist_ok=True)

    # Download model weights to outputs/
    print(f"Downloading mDeepFRI model weights v{args.model_version} ...", flush=True)
    result = subprocess.run(
        ["mDeepFRI", "get-models", "-o", "./outputs", "-v", args.model_version],
        text=True,
    )
    if result.returncode != 0:
        print(f"ERROR: Model download failed (exit code {result.returncode})", flush=True)
        sys.exit(result.returncode)

    # List downloaded model files
    model_files = [
        f for f in os.listdir("./outputs")
        if f.endswith(".onnx") or f.endswith(".json")
    ]
    print(f"Downloaded {len(model_files)} model file(s):", flush=True)
    for fname in sorted(model_files):
        size_mb = os.path.getsize(os.path.join("./outputs", fname)) / 1_048_576
        print(f"  {fname} ({size_mb:.1f} MB)", flush=True)

    # Pass validated.fasta through to the next node
    src = "./inputs/validated.fasta"
    dst = "./outputs/validated.fasta"
    shutil.copy(src, dst)
    with open(src) as f:
        n_seqs = sum(1 for line in f if line.startswith(">"))
    print(f"Passed through {n_seqs} validated sequence(s) -> {dst}", flush=True)

    print("Node 02 complete.", flush=True)


if __name__ == "__main__":
    main()
