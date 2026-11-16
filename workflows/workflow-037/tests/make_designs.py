#!/usr/bin/env python3
"""Synthesise minimal 'generated' structures so node 04 can be tested without a GPU.

Node 03 needs a GPU and minutes per design, which would put node 04 out of reach of a
suite anyone actually runs. But node 04 only reads chain ids, residue names, sequence
numbers and CA coordinates -- so a hand-built two-chain mmCIF exercises it exactly as a
real design does, and does it deterministically, which a sampled design never could.

Chain B is the design, chain A the target, matching the chain order RFD4 produces from
this workflow's contigs (designed segment first).

Run:  python3 -I tests/make_designs.py <out-dir> <n-designs>
"""
import os
import sys

COLUMNS = [
    "_atom_site.group_PDB",
    "_atom_site.id",
    "_atom_site.type_symbol",
    "_atom_site.label_atom_id",
    "_atom_site.label_alt_id",
    "_atom_site.label_comp_id",
    "_atom_site.label_asym_id",
    "_atom_site.label_entity_id",
    "_atom_site.label_seq_id",
    "_atom_site.pdbx_PDB_ins_code",
    "_atom_site.Cartn_x",
    "_atom_site.Cartn_y",
    "_atom_site.Cartn_z",
    "_atom_site.occupancy",
    "_atom_site.B_iso_or_equiv",
    "_atom_site.auth_seq_id",
    "_atom_site.auth_asym_id",
    "_atom_site.pdbx_PDB_model_num",
]

DESIGNED_SEQ = ["ALA", "LEU", "GLU", "LYS", "ILE", "ALA", "GLU", "LEU", "LYS", "ALA"]
TARGET_SEQ = ["GLY", "SER", "THR", "VAL", "PRO", "ASN", "ASP", "TYR", "PHE", "TRP",
              "GLY", "SER", "THR", "VAL", "PRO", "ASN", "ASP", "TYR", "PHE", "TRP"]


def chain_rows(chain_id, entity, residues, x0, start_id):
    """One CA per residue, 3.8 A apart along z: a plausible extended backbone."""
    rows, atom_id = [], start_id
    for i, comp in enumerate(residues, start=1):
        rows.append([
            "ATOM", str(atom_id), "C", "CA", ".", comp, chain_id, entity, str(i), "?",
            f"{x0:.3f}", "0.000", f"{i * 3.8:.3f}", "1.00", "20.00", str(i), chain_id, "1",
        ])
        atom_id += 1
    return rows, atom_id


def write_design(path, swap_chains=False):
    # 6 A apart: close enough to be an interface, far enough not to read as a clash.
    designed, next_id = chain_rows("B", "1", DESIGNED_SEQ, 0.0, 1)
    target, _ = chain_rows("A", "0", TARGET_SEQ, 6.0, next_id)
    blocks = [target, designed] if swap_chains else [designed, target]
    with open(path, "w") as handle:
        handle.write(f"data_{os.path.basename(path)}\n#\nloop_\n")
        for column in COLUMNS:
            handle.write(column + "\n")
        for block in blocks:
            for row in block:
                handle.write(" ".join(row) + "\n")
        handle.write("#\n")


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    out_dir, count = sys.argv[1], int(sys.argv[2])
    swap = len(sys.argv) > 3 and sys.argv[3] == "swap"
    os.makedirs(out_dir, exist_ok=True)
    for i in range(count):
        write_design(os.path.join(out_dir, f"design_{i}_generated.cif"), swap_chains=swap)
    print(f"wrote {count} design(s) to {out_dir}"
          + (" with chains in reversed order" if swap else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
