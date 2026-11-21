#!/usr/bin/env python3
"""Synthesise 'generated' structures so nodes 04-06 can be tested without a GPU.

Node 03 needs a GPU and minutes per design, which would put everything after it out of
reach of a suite anyone actually runs. But nodes 04, 05 and 06 read only chain ids,
residue names, sequence numbers and CA coordinates, so a hand-built mmCIF exercises them
exactly as a real design does -- and deterministically, which a sampled design cannot be.

Chain B is the design, chain A the target, matching the chain order RFD4 produces from
this workflow's contigs (designed segment first).

Shapes, chosen so the geometry checks in node 04 have something real to fail on:

  compact   a serpentine walk through a cubic lattice at 3.8 A. Every consecutive pair
            is one lattice step apart, so the CA-CA distances are exactly right and the
            whole thing is globular. This is what a good design looks like.
  extended  a straight line. CA-CA distances are still correct -- it is the shape that
            is wrong, which is precisely the failure a broken checkpoint produces.
  helix     a textbook alpha-helix: CA radius 2.3 A, rise 1.5 A, 100 deg per residue.
            Neither compact nor broken -- a real secondary structure, and the shape a
            single-helix binder arrives in. The lattice and straight-line shapes between
            them cannot produce it, so without this shape nothing in the suite asks what
            the compactness check does to a fold that is genuinely elongated and
            genuinely folded. Cross-checked against 2ZTA chain A, a 32-residue leucine
            zipper helix: measured Rg 13.67 A, this shape at 32 gives 14.04 A.
  broken    compact, with one residue displaced to open a gap no peptide bond spans.
  clash     a compact design interpenetrating the target rather than docked against it.
  swap      compact, but with the chains written target-first.

Run:  python3 -I tests/make_designs.py <out-dir> <n-designs> [shape] [length]
"""
import math
import os
import sys

COLUMNS = [
    "_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol",
    "_atom_site.label_atom_id", "_atom_site.label_alt_id", "_atom_site.label_comp_id",
    "_atom_site.label_asym_id", "_atom_site.label_entity_id", "_atom_site.label_seq_id",
    "_atom_site.pdbx_PDB_ins_code", "_atom_site.Cartn_x", "_atom_site.Cartn_y",
    "_atom_site.Cartn_z", "_atom_site.occupancy", "_atom_site.B_iso_or_equiv",
    "_atom_site.auth_seq_id", "_atom_site.auth_asym_id", "_atom_site.pdbx_PDB_model_num",
]

CA_STEP_A = 3.8
# Repeating, so a design of any length has a defined sequence and node 06 can be checked
# against it exactly.
DESIGN_CYCLE = ["ALA", "LEU", "GLU", "LYS", "ILE", "ALA", "GLU", "LEU", "LYS", "ALA"]
TARGET_CYCLE = ["GLY", "SER", "THR", "VAL", "PRO", "ASN", "ASP", "TYR", "PHE", "TRP"]
ONE_LETTER = {
    "ALA": "A", "LEU": "L", "GLU": "E", "LYS": "K", "ILE": "I",
    "GLY": "G", "SER": "S", "THR": "T", "VAL": "V", "PRO": "P",
    "ASN": "N", "ASP": "D", "TYR": "Y", "PHE": "F", "TRP": "W",
}


def lattice_walk(n, step=CA_STEP_A, origin=(0.0, 0.0, 0.0)):
    """A self-avoiding serpentine path through a cube: every step is one lattice edge.

    x reverses at the end of each row, y at the end of each layer, so consecutive points
    always differ in exactly one coordinate by exactly one step.
    """
    side = 1
    while side ** 3 < n:
        side += 1
    points, x, y, dx, dy = [], 0, 0, 1, 1
    for z in range(side):
        for row in range(side):
            for column in range(side):
                points.append((x, y, z))
                if column < side - 1:
                    x += dx
            dx = -dx
            if row < side - 1:
                y += dy
        dy = -dy
    ox, oy, oz = origin
    return [(ox + px * step, oy + py * step, oz + pz * step) for px, py, pz in points[:n]]


def lattice_extent(n, step=CA_STEP_A):
    """How far a lattice walk of n points reaches along one axis."""
    side = 1
    while side ** 3 < n:
        side += 1
    return (side - 1) * step


def straight_walk(n, step=CA_STEP_A, origin=(0.0, 0.0, 0.0)):
    ox, oy, oz = origin
    return [(ox, oy, oz + i * step) for i in range(n)]


# Textbook alpha-helix CA geometry. The three together put consecutive alpha carbons
# 3.83 A apart, inside the 2.8-4.8 A window node 04 accepts, so this shape exercises the
# compactness check with the chain-break check deliberately silent.
HELIX_RADIUS_A = 2.3
HELIX_RISE_A = 1.5
HELIX_TWIST_DEG = 100.0


def helix_walk(n, origin=(0.0, 0.0, 0.0)):
    ox, oy, oz = origin
    points = []
    for i in range(n):
        theta = math.radians(HELIX_TWIST_DEG * i)
        points.append((ox + HELIX_RADIUS_A * math.cos(theta),
                       oy + HELIX_RADIUS_A * math.sin(theta),
                       oz + HELIX_RISE_A * i))
    return points


def rows_for(chain_id, entity, cycle, points, start_id):
    rows, atom_id = [], start_id
    for i, point in enumerate(points, start=1):
        comp = cycle[(i - 1) % len(cycle)]
        rows.append([
            "ATOM", str(atom_id), "C", "CA", ".", comp, chain_id, entity, str(i), "?",
            f"{point[0]:.3f}", f"{point[1]:.3f}", f"{point[2]:.3f}",
            "1.00", "20.00", str(i), chain_id, "1",
        ])
        atom_id += 1
    return rows, atom_id


def design_sequence(length):
    return "".join(ONE_LETTER[DESIGN_CYCLE[i % len(DESIGN_CYCLE)]] for i in range(length))


def write_design(path, shape, length, target_length):
    if shape == "extended":
        designed = straight_walk(length)
    elif shape == "helix":
        designed = helix_walk(length)
    else:
        designed = lattice_walk(length)
    if shape == "broken" and len(designed) > 2:
        middle = len(designed) // 2
        x, y, z = designed[middle]
        # 12 A from its neighbour: no peptide bond reaches that far.
        designed[middle] = (x + 12.0, y, z)

    # The target sits beside the design, or inside it for the clash case. The gap has to
    # be measured from the design's own extent: a fixed offset left only 2.6 A at
    # length 30, so every design came out flagged as clashing and the ranking case was
    # testing nothing. 6 A clears the 3.8 A clash threshold and still falls inside the
    # 8 A contact cutoff, so the designs have a real interface to be scored on.
    # A helix reaches along z, not across the lattice box, so its own extent along x is
    # one radius. Using the lattice figure here would push the target 80 A clear of a
    # 50-mer helix and leave the design with no interface to be scored on at all.
    if shape == "helix":
        separation = HELIX_RADIUS_A + 6.0
    else:
        separation = lattice_extent(length) + 6.0
    offset = (1.0, 0.0, 0.0) if shape == "clash" else (separation, 0.0, 0.0)
    target = lattice_walk(target_length, origin=offset)

    design_rows, next_id = rows_for("B", "1", DESIGN_CYCLE, designed, 1)
    target_rows, _ = rows_for("A", "0", TARGET_CYCLE, target, next_id)
    blocks = [target_rows, design_rows] if shape == "swap" else [design_rows, target_rows]

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
    shape = sys.argv[3] if len(sys.argv) > 3 else "compact"
    length = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    target_length = max(length * 2, 20)
    os.makedirs(out_dir, exist_ok=True)
    for i in range(count):
        # Only the first design carries the shape; the rest are ordinary, so a case can
        # ask whether one bad design is ranked and reported correctly among good ones.
        this = shape if i == 0 else "compact"
        write_design(os.path.join(out_dir, f"design_{i}_generated.cif"),
                     this, length, target_length)
    print(f"wrote {count} design(s) to {out_dir}: shape={shape} length={length} "
          f"sequence={design_sequence(length)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
