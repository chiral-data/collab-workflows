#!/usr/bin/env python3
"""Measure a design CIF independently of every node that reads one.

Node 01 is the thing being checked, so the suite cannot use node 01's parser to decide
what node 01 should have found. This file hand-parses the _atom_site loop and shares no
code with any node -- the same separation 037's probe_entry.py keeps, and for the same
reason: a check that measures its subject with its subject's own instrument passes for
the wrong reason (§8.4).

What it reports, per file, per chain, in FILE ORDER:

    chain          label_asym_id
    n_residues     distinct label_seq_id values
    first_id       the lowest, which is not always 1
    last_id
    gaps           [(after, before, missing), ...] from the numbering
    sequence       one-letter, in label order, X where there is no code
    n_unknown      how many of those X there are
    n_ca           residues carrying a CA atom
    ca_md5         a digest of the CA coordinates, in order

ca_md5 is there so a check can ask whether node 01 recorded the SAME coordinates in the
SAME order, without this file and node 01 having to agree on a float format.

Run:  python3 -I tests/probe_design.py <dir-or-file> [--json]

-I because the files may be downloaded data and must not be able to put a module on the
import path.
"""
import glob
import hashlib
import json
import os
import sys

# Deliberately the same table node 01 uses, because the MAPPING is not what is being
# checked here -- the chains, the order, the numbering and the coordinates are. If this
# table and node 01's ever disagree the suite would report a sequence mismatch, which is
# the right answer to that disagreement.
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}


def atom_site(path):
    """({column_name: index}, [row_fields]) for the _atom_site loop, parsed by hand.

    Hand-parsed rather than through Biopython so this tool can run outside the container
    and cannot be fooled by the same parser the nodes use. Handles the one-loop,
    space-separated form that both RFD4 and tests/make_designs.py emit; it is not a
    general mmCIF reader and does not pretend to be.
    """
    columns, rows, state = [], [], "search"
    with open(path, errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if state == "search":
                if line.startswith("_atom_site."):
                    columns = [line]
                    state = "head"
                continue
            if state == "head":
                if line.startswith("_atom_site."):
                    columns.append(line)
                    continue
                state = "rows"
            if state == "rows":
                if not line or line.startswith(("#", "loop_", "_", "data_")):
                    break
                fields = line.split()
                if len(fields) == len(columns):
                    rows.append(fields)
    index = {name.split(".", 1)[1]: i for i, name in enumerate(columns)}
    return index, rows


def probe(path):
    index, rows = atom_site(path)
    if not rows:
        return {"file": os.path.basename(path), "error": "no _atom_site rows"}

    def field(row, name, default="."):
        at = index.get(name)
        return row[at] if at is not None else default

    order, chains = [], {}
    for row in rows:
        seq_raw = field(row, "label_seq_id")
        if seq_raw in (".", "?", ""):
            continue
        try:
            seq_id = int(seq_raw)
        except ValueError:
            continue
        chain = field(row, "label_asym_id")
        if chain not in chains:
            chains[chain] = {}
            order.append(chain)
        entry = chains[chain].setdefault(
            seq_id, {"comp": field(row, "label_comp_id"), "ca": None}
        )
        if field(row, "label_atom_id") == "CA":
            try:
                entry["ca"] = (
                    float(field(row, "Cartn_x", "0")),
                    float(field(row, "Cartn_y", "0")),
                    float(field(row, "Cartn_z", "0")),
                )
            except ValueError:
                entry["ca"] = None

    out = []
    for chain in order:
        ids = sorted(chains[chain])
        residues = [chains[chain][i] for i in ids]
        gaps = [
            {"after": a, "before": b, "missing": b - a - 1}
            for a, b in zip(ids, ids[1:]) if b - a > 1
        ]
        sequence = "".join(THREE_TO_ONE.get(r["comp"], "X") for r in residues)
        digest = hashlib.md5()
        for r in residues:
            if r["ca"] is None:
                digest.update(b"none|")
            else:
                digest.update(("%.3f,%.3f,%.3f|" % r["ca"]).encode())
        out.append({
            "chain": chain,
            "n_residues": len(ids),
            "first_id": ids[0],
            "last_id": ids[-1],
            "residue_ids": ids,
            "gaps": gaps,
            "sequence": sequence,
            "n_unknown": sequence.count("X"),
            "n_ca": sum(1 for r in residues if r["ca"] is not None),
            "ca_md5": digest.hexdigest(),
        })
    return {"file": os.path.basename(path), "chains": out}


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    target = sys.argv[1]
    as_json = "--json" in sys.argv[2:]
    paths = (
        sorted(glob.glob(os.path.join(target, "*.cif")))
        if os.path.isdir(target) else [target]
    )
    if not paths:
        print(f"probe_design: no .cif files in {target}", file=sys.stderr)
        return 1

    results = [probe(path) for path in paths]
    if as_json:
        print(json.dumps(results, indent=2))
        return 0
    for result in results:
        print(result["file"])
        if "error" in result:
            print(f"  ERROR: {result['error']}")
            continue
        for chain in result["chains"]:
            gaps = ", ".join(
                "%d->%d (%d missing)" % (g["after"], g["before"], g["missing"])
                for g in chain["gaps"]
            ) or "none"
            print(
                "  chain %-3s n=%-4d ids %d..%d  gaps: %s" %
                (chain["chain"], chain["n_residues"], chain["first_id"],
                 chain["last_id"], gaps)
            )
            print(
                "      ca=%d/%d  X=%d  seq=%s" %
                (chain["n_ca"], chain["n_residues"], chain["n_unknown"],
                 chain["sequence"][:60] + ("..." if len(chain["sequence"]) > 60 else ""))
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
