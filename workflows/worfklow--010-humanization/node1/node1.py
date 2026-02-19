"""
Node 1 – Input Chain Preparation
=================================
Reads heavy and light chain FASTA files, validates sequences,
pairs them by order, converts each pair to an ABB3 ab_input dict
using abodybuilder3's string_to_input utility, and saves each
pair as a .pt file for downstream nodes.

Outputs
-------
results/node1/<pair_id>.pt
    Each file contains a dict:
        {
            "id":       str,
            "heavy":    str,
            "light":    str,
            "ab_input": Dict[str, Tensor],   # from string_to_input()
        }
"""

from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import torch


# ---------------------------------------------------------------------------
# FASTA helpers
# ---------------------------------------------------------------------------

def read_fasta(path: Path) -> List[Tuple[str, str]]:
    """Return list of (header, sequence) from a FASTA file."""
    records: List[Tuple[str, str]] = []
    cur_id: Optional[str] = None
    cur_seq: List[str] = []

    with path.open("r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_id is not None:
                    records.append((cur_id, "".join(cur_seq)))
                cur_id = line[1:].split()[0]   # use first token as ID
                cur_seq = []
            else:
                cur_seq.append(line.upper())

    if cur_id is not None:
        records.append((cur_id, "".join(cur_seq)))

    return records


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def validate_seq(seq: str, label: str) -> str:
    """Validate sequence; return uppercased, stripped sequence on success."""
    seq = seq.upper().strip()
    if not seq:
        raise ValueError(f"{label}: sequence is empty.")
    bad = set(c for c in seq if c not in VALID_AA)
    if bad:
        raise ValueError(
            f"{label}: sequence contains invalid characters: {', '.join(sorted(bad))}. "
            f"Only standard 20 amino-acid single-letter codes are accepted."
        )
    return seq


# ---------------------------------------------------------------------------
# Pairing
# ---------------------------------------------------------------------------

def pair_by_order(
    heavy: List[Tuple[str, str]],
    light: List[Tuple[str, str]],
) -> List[Tuple[str, str, str]]:
    """
    Pair heavy and light chains by position.
    Returns list of (pair_id, heavy_seq, light_seq).
    pair_id is constructed as  <heavy_id>-<light_id>.
    """
    if len(heavy) != len(light):
        raise ValueError(
            f"Heavy FASTA has {len(heavy)} sequences but light FASTA has "
            f"{len(light)}.  They must be the same length for order-based pairing."
        )
    pairs = []
    for (h_id, h_seq), (l_id, l_seq) in zip(heavy, light):
        pair_id = f"{h_id}-{l_id}"
        pairs.append((pair_id, h_seq, l_seq))
    return pairs


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Node 1: Read, validate, and pair heavy/light FASTA files."
    )
    parser.add_argument("--heavy_fasta", required=True, help="Path to heavy-chain FASTA file")
    parser.add_argument("--light_fasta", required=True, help="Path to light-chain FASTA file")
    parser.add_argument("--out_dir", required=True, help="Directory to write .pt output files")
    args = parser.parse_args()

    heavy_path = Path(args.heavy_fasta)
    light_path = Path(args.light_fasta)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1. Read
    print(f"[Node 1] Reading heavy chains from: {heavy_path}")
    heavy_records = read_fasta(heavy_path)
    print(f"[Node 1] Reading light chains from: {light_path}")
    light_records = read_fasta(light_path)

    print(f"[Node 1] Found {len(heavy_records)} heavy / {len(light_records)} light sequences.")

    # 2. Validate
    heavy_validated = []
    for h_id, h_seq in heavy_records:
        seq = validate_seq(h_seq, label=f"Heavy [{h_id}]")
        heavy_validated.append((h_id, seq))

    light_validated = []
    for l_id, l_seq in light_records:
        seq = validate_seq(l_seq, label=f"Light [{l_id}]")
        light_validated.append((l_id, seq))

    # 3. Pair
    pairs = pair_by_order(heavy_validated, light_validated)
    print(f"[Node 1] Paired {len(pairs)} antibody chain(s).")

    # 4. Convert to ab_input and save
    # Import here so the script fails informatively if abodybuilder3 is missing.
    try:
        from abodybuilder3.utils import string_to_input
    except ImportError as e:
        raise ImportError(
            "abodybuilder3 is not installed or not on the Python path. "
            "Please install it before running this node."
        ) from e

    for pair_id, h_seq, l_seq in pairs:
        print(f"[Node 1]   Processing pair: {pair_id}")

        ab_input = string_to_input(heavy=h_seq, light=l_seq)

        out_path = out_dir / f"{pair_id}.pt"
        torch.save(
            {
                "id": pair_id,
                "heavy": h_seq,
                "light": l_seq,
                "ab_input": ab_input,
            },
            out_path,
        )
        print(f"[Node 1]   Saved → {out_path}")

    print(f"[Node 1] Done. {len(pairs)} file(s) written to {out_dir}")


if __name__ == "__main__":
    main()