#!/usr/bin/env python3
"""Download mDeepFRI model weights (GCN + CNN) and pass through validated FASTA.

Node 02 prepares the structural alignment environment: the GCN models perform
graph-based inference on contact maps derived from FoldComp structural alignments,
while CNN models handle sequence-only prediction for proteins without structural hits.
"""

import argparse
import os
import shutil
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description="Download mDeepFRI model weights")
    parser.add_argument("--model-version", default="1.1", choices=["1.0", "1.1"])
    args = parser.parse_args()

    # mDeepFRI get-models requires a non-existing or empty target directory.
    # Download to a subdir, then flatten model files to CWD for silva output collection.
    weights_subdir = "./models"
    if os.path.exists(weights_subdir):
        shutil.rmtree(weights_subdir)

    print(f"Downloading mDeepFRI model weights v{args.model_version} ...", flush=True)
    result = subprocess.run(
        ["mDeepFRI", "get-models", "-o", weights_subdir, "-v", args.model_version],
        text=True,
    )
    if result.returncode != 0:
        print(f"ERROR: Model download failed (exit code {result.returncode})", flush=True)
        sys.exit(result.returncode)

    # Flatten model files to CWD so silva's output glob (*.onnx, *.json) finds them
    model_files = [f for f in os.listdir(weights_subdir)
                   if f.endswith(".onnx") or f.endswith(".json")]
    print(f"Copying {len(model_files)} model file(s) to output:", flush=True)
    for fname in sorted(model_files):
        src = os.path.join(weights_subdir, fname)
        shutil.copy(src, f"./{fname}")
        size_mb = os.path.getsize(f"./{fname}") / 1_048_576
        print(f"  {fname} ({size_mb:.1f} MB)", flush=True)

    shutil.rmtree(weights_subdir, ignore_errors=True)

    # Pass validated.fasta through to the next node
    src = "./inputs/validated.fasta"
    shutil.copy(src, "./validated.fasta")
    with open(src) as f:
        n_seqs = sum(1 for line in f if line.startswith(">"))
    print(f"Passed through {n_seqs} validated sequence(s)", flush=True)

    print("Node 02 complete.", flush=True)


if __name__ == "__main__":
    main()
