#!/usr/bin/env python3
import os
import sys

from Bio import SeqIO

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYX")


def main():
    os.makedirs("./outputs", exist_ok=True)

    min_len = int(os.environ.get("PARAM_MIN_LENGTH", "10"))
    max_len = int(os.environ.get("PARAM_MAX_LENGTH", "2500"))

    input_path = "./inputs/sequences.fasta"
    if not os.path.exists(input_path):
        print(f"ERROR: Input file not found: {input_path}", flush=True)
        sys.exit(1)

    records = list(SeqIO.parse(input_path, "fasta"))
    if not records:
        print("ERROR: No sequences found in FASTA file", flush=True)
        sys.exit(1)

    # Detect monomer vs. multimer
    is_multimer = len(records) > 1 or any(":" in r.id or ":" in r.description for r in records)
    mode = "multimer" if is_multimer else "monomer"

    errors = []
    for rec in records:
        seq = str(rec.seq).upper()
        bad_chars = set(seq) - VALID_AA
        if bad_chars:
            errors.append(
                f"Sequence '{rec.id}': invalid characters {sorted(bad_chars)}"
            )
        if len(seq) < min_len:
            errors.append(
                f"Sequence '{rec.id}': length {len(seq)} is below minimum {min_len}"
            )
        if len(seq) > max_len:
            errors.append(
                f"Sequence '{rec.id}': length {len(seq)} exceeds maximum {max_len}"
            )

    if errors:
        for msg in errors:
            print(f"ERROR: {msg}", flush=True)
        sys.exit(1)

    print(f"Sequences: {len(records)}", flush=True)
    print(f"Mode: {mode}", flush=True)
    for rec in records:
        print(f"  {rec.id}: {len(rec.seq)} aa", flush=True)

    SeqIO.write(records, "./outputs/validated_sequences.fasta", "fasta")
    print("Validation passed. Output: ./outputs/validated_sequences.fasta", flush=True)


if __name__ == "__main__":
    main()
