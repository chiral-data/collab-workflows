#!/usr/bin/env python3
"""Validate FASTA inputs for mDeepFRI — check format, amino acid content, and sequence length."""

import argparse
import json
import os
import sys

# Standard amino acids + IUPAC ambiguity codes + gap/stop
VALID_AA = set("ACDEFGHIKLMNPQRSTVWYBXZJUO*-")
NUCLEOTIDE_ONLY = set("ACGTU")


def parse_fasta(path):
    """Yield (header, sequence) from a FASTA file."""
    header, parts = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header = line[1:].strip()
                parts = []
            elif line:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def is_protein(seq):
    """Return True if the sequence is likely a protein (not nucleotide)."""
    seq_upper = seq.upper().replace("-", "").replace("*", "")
    if not seq_upper:
        return False
    unique = set(seq_upper)
    # Purely ACGTU characters → likely nucleotide
    if unique.issubset(NUCLEOTIDE_ONLY):
        return False
    # If >90% of unique chars are nucleotide chars, still flag as nucleotide
    nt_overlap = unique & NUCLEOTIDE_ONLY
    if len(nt_overlap) / max(len(unique), 1) > 0.9 and len(unique) <= 5:
        return False
    return True


def main():
    parser = argparse.ArgumentParser(description="Validate FASTA protein sequences")
    parser.add_argument("--min-length", type=int, default=10, help="Minimum sequence length")
    parser.add_argument("--max-length", type=int, default=5000, help="Maximum sequence length")
    args = parser.parse_args()

    if not os.path.isdir("./inputs"):
        print("ERROR: ./inputs directory not found", flush=True)
        sys.exit(1)

    input_files = []
    for fname in sorted(os.listdir("./inputs")):
        if fname.lower().endswith((".fasta", ".fa", ".faa")):
            input_files.append(os.path.join("./inputs", fname))

    if not input_files:
        print("ERROR: No FASTA files found in ./inputs/", flush=True)
        sys.exit(1)

    print(f"Found {len(input_files)} input file(s): {', '.join(os.path.basename(f) for f in input_files)}", flush=True)
    print(f"Length filter: {args.min_length}–{args.max_length} aa", flush=True)

    total = valid = too_short = too_long = not_protein = invalid_chars = 0
    seen_ids = set()
    duplicates = 0

    with open("./validated.fasta", "w") as out:
        for fpath in input_files:
            for header, seq in parse_fasta(fpath):
                total += 1
                seq_clean = seq.upper().replace(" ", "").replace("\t", "")

                if not is_protein(seq_clean):
                    not_protein += 1
                    print(f"  SKIP (not protein): {header[:60]}", flush=True)
                    continue

                bad_chars = set(seq_clean) - VALID_AA
                if bad_chars:
                    invalid_chars += 1
                    print(f"  SKIP (invalid chars {bad_chars}): {header[:60]}", flush=True)
                    continue

                # Strip gap/stop for length check
                seq_residues = seq_clean.replace("-", "").replace("*", "")

                if len(seq_residues) < args.min_length:
                    too_short += 1
                    continue
                if len(seq_residues) > args.max_length:
                    too_long += 1
                    continue

                # Dedup by header ID (first word)
                seq_id = header.split()[0]
                if seq_id in seen_ids:
                    duplicates += 1
                    continue
                seen_ids.add(seq_id)

                # Normalize Swiss-Prot/TrEMBL headers (sp|ACC|NAME or tr|ACC|NAME)
                # MMseqs2 extracts only ACC from these, so we normalize to match.
                parts = seq_id.split("|")
                out_id = parts[1] if len(parts) == 3 and parts[0] in ("sp", "tr") else seq_id

                # Write in 80-char lines
                out.write(f">{out_id}\n")
                for i in range(0, len(seq_residues), 80):
                    out.write(seq_residues[i:i + 80] + "\n")
                valid += 1

    report = {
        "total_input": total,
        "valid": valid,
        "filtered_not_protein": not_protein,
        "filtered_invalid_chars": invalid_chars,
        "filtered_too_short": too_short,
        "filtered_too_long": too_long,
        "filtered_duplicates": duplicates,
        "min_length": args.min_length,
        "max_length": args.max_length,
        "input_files": [os.path.basename(f) for f in input_files],
    }
    with open("./validation_report.json", "w") as f:
        json.dump(report, f, indent=2)

    print(f"\nValidation summary: {valid}/{total} sequences passed", flush=True)
    for k, v in report.items():
        if k.startswith("filtered_") and v > 0:
            print(f"  {k}: {v}", flush=True)

    if valid == 0:
        print("ERROR: No valid sequences remain after filtering", flush=True)
        sys.exit(1)

    print("Validation complete.", flush=True)


if __name__ == "__main__":
    main()
