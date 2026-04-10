#!/usr/bin/env python3
"""Run mDeepFRI predict-function on validated protein sequences.

mDeepFRI pipeline:
1. MMseqs2 searches the FoldComp structural database(s) for similar proteins.
2. Structural hits → contact maps → GCN prediction (structure-aware).
3. Sequence-only proteins → CNN prediction (sequence-based fallback).
4. Results written to results.tsv and alignment_summary.tsv.
"""

import os
import shutil
import subprocess
import sys
import tempfile


def main():
    # Parameters from silva environment variables
    modes_str = os.environ.get("PARAM_PREDICTION_MODES", "mf bp cc")
    skip_pdb = os.environ.get("PARAM_SKIP_PDB", "true").lower() in ("true", "1", "yes")
    threads = os.environ.get("PARAM_THREADS", "1")
    sensitivity = os.environ.get("PARAM_MMSEQS_SENSITIVITY", "5.7")

    fasta = "./inputs/validated.fasta"
    if not os.path.exists(fasta):
        print("ERROR: validated.fasta not found in ./inputs/", flush=True)
        sys.exit(1)

    with open(fasta) as f:
        n_seqs = sum(1 for line in f if line.startswith(">"))
    print(f"Input: {n_seqs} protein sequence(s)", flush=True)

    # Verify model weight files are present
    onnx_files = [f for f in os.listdir("./inputs") if f.endswith(".onnx")]
    if not onnx_files:
        print("ERROR: No .onnx model files found in ./inputs/", flush=True)
        sys.exit(1)
    print(f"Model weights: {len(onnx_files)} .onnx file(s) found", flush=True)

    # Use a temp subdir for mDeepFRI output to avoid conflicts with silva's output collection
    outdir = tempfile.mkdtemp(prefix="mdf_out_", dir=".")

    # Build predict-function command
    cmd = [
        "mDeepFRI", "predict-function",
        "-i", fasta,
        "-w", "./inputs",
        "-o", outdir,
        "-t", threads,
        "-s", sensitivity,
        "--remove-intermediate",
        "--skip-matrix",
    ]

    # Add prediction modes (mf, bp, cc, ec)
    for mode in modes_str.split():
        mode = mode.strip().lower()
        if mode:
            cmd.extend(["-p", mode])

    if skip_pdb:
        cmd.append("--skip-pdb")

    print(f"\nRunning: {' '.join(cmd)}\n", flush=True)
    result = subprocess.run(cmd, text=True)

    if result.returncode != 0:
        print(f"ERROR: mDeepFRI failed (exit code {result.returncode})", flush=True)
        shutil.rmtree(outdir, ignore_errors=True)
        sys.exit(result.returncode)

    # Copy output files to CWD (silva collects from CWD, not subdirs)
    copied = []
    for fname in ("results.tsv", "alignment_summary.tsv"):
        src = os.path.join(outdir, fname)
        if os.path.exists(src):
            shutil.copy(src, f"./{fname}")
            with open(f"./{fname}") as f:
                n_data = max(0, sum(1 for _ in f) - 1)
            print(f"  {fname}: {n_data} data row(s)", flush=True)
            copied.append(fname)
        else:
            print(f"  WARNING: {fname} not found in mDeepFRI output", flush=True)

    shutil.rmtree(outdir, ignore_errors=True)

    if not copied:
        print("ERROR: No output files produced by mDeepFRI", flush=True)
        sys.exit(1)

    print("Prediction complete.", flush=True)


if __name__ == "__main__":
    main()
