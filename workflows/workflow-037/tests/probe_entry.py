#!/usr/bin/env python3
"""Measure the input shapes that break this pipeline, from an mmCIF's _atom_site loop.

Every number this prints is read from the file. Nothing here is a claim from the
literature or from memory -- that is the point. The fixture manifest records what this
tool reported, so a reviewer can re-run it and check that the fixtures still carry the
shapes they were chosen for (RCSB does re-release entries).

Run:  python3 -I tests/probe_entry.py <dir-of-cif-files>
      python3 -I tests/probe_entry.py <one.cif>

-I because the .cif files are downloaded data and must not be able to put a module on
the import path.
"""
import glob
import os
import sys
from collections import defaultdict

STANDARD = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

KEYS = ["models", "chains", "ins_codes", "dup_auth", "modified", "het_in_poly",
        "nonpolymer", "zero_occ", "nan_xyz", "altlocs"]


def atom_site(path):
    """Return ({column: index}, [row fields]) for the _atom_site loop, parsed by hand.

    Hand-parsed rather than via Biopython so the harness can run outside the container
    and cannot be fooled by the same parser the nodes use.
    """
    cols, rows, state = [], [], "search"
    with open(path, errors="replace") as handle:
        for line in handle:
            s = line.strip()
            if state == "search":
                if s.startswith("_atom_site."):
                    cols = [s]
                    state = "head"
                continue
            if state == "head":
                if s.startswith("_atom_site."):
                    cols.append(s)
                    continue
                state = "body"
            if state == "body":
                if s.startswith("#") or s.startswith("loop_") or s.startswith("_"):
                    break
                fields = s.split()
                if len(fields) == len(cols):
                    rows.append(fields)
    return {name: i for i, name in enumerate(cols)}, rows


def probe(path):
    """Measured shape of one mmCIF. Raises nothing; returns {'parse': 'FAILED'} instead."""
    idx, rows = atom_site(path)
    if not rows:
        return {"parse": "FAILED"}

    def col(name):
        i = idx.get("_atom_site." + name)
        return [r[i] for r in rows] if i is not None else [""] * len(rows)

    group, label_seq = col("group_PDB"), col("label_seq_id")
    auth_seq, ins = col("auth_seq_id"), col("pdbx_PDB_ins_code")
    model, occ = col("pdbx_PDB_model_num"), col("occupancy")
    comp, asym = col("label_comp_id"), col("label_asym_id")
    altloc, xs = col("label_alt_id"), col("Cartn_x")

    blank = (".", "?", "")
    poly = [i for i in range(len(rows)) if label_seq[i] not in blank]
    modified = sorted({comp[i] for i in poly if comp[i] not in STANDARD})
    nonpoly = sorted({comp[i] for i in range(len(rows))
                      if label_seq[i] in blank and comp[i] not in ("HOH", "")})
    ins_codes = sorted({ins[i] for i in range(len(rows)) if ins[i] not in blank})

    # author numbers covering more than one residue of the same chain = insertion codes
    seen = defaultdict(set)
    for i in poly:
        seen[(asym[i], auth_seq[i])].add(label_seq[i])
    dup_auth = sum(1 for v in seen.values() if len(v) > 1)

    def is_zero(value):
        try:
            return float(value) == 0.0
        except ValueError:
            return False

    models = sorted({model[i] for i in range(len(rows)) if model[i] not in blank})
    alts = sorted({altloc[i] for i in range(len(rows)) if altloc[i] not in blank})

    return {
        "models": len(models) or 1,
        "chains": len({asym[i] for i in poly}),
        "n_atom_rows": len(rows),
        "ins_codes_list": ins_codes,
        "nonpolymer_list": nonpoly,
        "modified_list": modified,
        "chain_ids": sorted({asym[i] for i in poly}),
        "ins_codes": f"{len(ins_codes)} ({','.join(ins_codes[:4])})" if ins_codes else "-",
        "dup_auth": dup_auth or "-",
        "modified": f"{len(modified)} ({','.join(modified[:4])})" if modified else "-",
        "het_in_poly": sum(1 for i in poly if group[i] == "HETATM") or "-",
        "nonpolymer": f"{len(nonpoly)} ({','.join(nonpoly[:4])})" if nonpoly else "-",
        "zero_occ": sum(1 for i in range(len(rows)) if is_zero(occ[i])) or "-",
        "nan_xyz": sum(1 for i in range(len(rows))
                       if xs[i].lower() in ("nan", "-nan", "inf", "-inf")) or "-",
        "altlocs": ",".join(alts) if alts else "-",
    }


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    target = sys.argv[1]
    paths = ([target] if target.lower().endswith(".cif")
             else sorted(glob.glob(os.path.join(target, "*.cif"))))
    if not paths:
        print(f"No .cif files under {target}")
        return 1
    print(f"{'entry':8s} " + " ".join(f"{k:>14s}" for k in KEYS))
    for path in paths:
        name = os.path.basename(path).replace(".cif", "")
        result = probe(path)
        if result.get("parse") == "FAILED":
            print(f"{name:8s} PARSE FAILED")
            continue
        print(f"{name:8s} " + " ".join(f"{str(result[k]):>14s}" for k in KEYS))
    return 0


if __name__ == "__main__":
    sys.exit(main())
