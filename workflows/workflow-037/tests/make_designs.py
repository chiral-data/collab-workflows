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

Added for workflow-040, which reads the same structures (the five shapes above describe
the DESIGN; these five describe the TARGET or the residue names, which is what 040's
ingest node has to get right):

  monomer   the design chain alone, no target. What a monomer or motif_scaffold run comes
            back as, and the only input for which 'designed_chain_source = single-chain'
            is the correct answer rather than a fallback. 040 must accept it on the
            self-consistency path and refuse it on the complex path -- there is no target
            to score an interface against.
  multichain
            one design chain against TWO target chains, A and C. A real target is
            routinely more than one chain: 037's own 1BRS fixture is three copies of the
            barnase-barstar complex, and five of its six chains vanished silently before
            53d1d460. The chain ids are deliberately not consecutive, so nothing can pass
            by assuming the targets are the alphabetical run after the design.
  gapped    a target chain with an internal numbering gap, and the second segment
            displaced so the junction is not a peptide bond -- which is what a
            disordered loop absent from the coordinates actually looks like. Ordinary in
            anything from RCSB. ColabFold cannot express a gap inside a chain, and the
            measured answer is to fold the chain whole and disclose the fusion: trimming
            to the segment carrying the epitope gave ~42 pLDDT where the whole chain
            reaches 81 (8e4772a7).
  renumbered
            a target chain whose label numbering starts at 3, not 1. Three renumberings
            stack in 040 -- the deposited entry's, RFD4's renumber-from-1, and
            ColabFold's own renumber-from-1-per-chain -- so any correspondence
            established by matching residue NUMBERS is wrong somewhere. Residue ids are
            recorded in chain order and the mapping is by position (ef74646e).
  nonstandard
            a design chain carrying modified residues with no one-letter code, taken from
            037's 5O45 fixture (9KK, CCS, MEA, NH2, SAR). They score X, and X has to be
            COUNTED: ColabFold folds X as a sequence gap, which depresses confidence for
            reasons that have nothing to do with the design, so a low pLDDT here must not
            be reported as a design defect.
  allnonstandard
            every residue of the design modified, so the sequence is entirely X. There is
            nothing for a predictor to fold, and the design must be skipped and marked
            rather than sent to the GPU (54bf21b6).

--manifest also writes generation_manifest.json beside the structures, carrying the
target chain ids this generator actually used -- which is the field node 03 records and
node 04 resolves the designed chain from. 037's own run_fixtures.sh writes that file
inline with target_chains = ["A"], which is correct for the five original shapes and
wrong for 'multichain' and 'monomer'; the flag is the general form, and it means the
harness does not restate a fact this file decides.

Run:  python3 -I tests/make_designs.py <out-dir> <n-designs> [shape] [length] [--manifest]
"""
import json
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
# Modified polymer components from 037's 5O45 fixture: real PDB chemical component ids
# with no one-letter code, so a parser must score them X rather than recognise them.
NONSTANDARD_CYCLE = ["9KK", "CCS", "MEA", "NH2", "SAR"]
# A few of them among ordinary residues: enough that the unknown count is non-zero while
# the sequence is still foldable, which is the case that must be counted rather than
# skipped.
MIXED_CYCLE = ["ALA", "9KK", "LEU", "GLU", "CCS", "LYS", "ILE", "SAR", "GLU", "ALA"]
ONE_LETTER = {
    "ALA": "A", "LEU": "L", "GLU": "E", "LYS": "K", "ILE": "I",
    "GLY": "G", "SER": "S", "THR": "T", "VAL": "V", "PRO": "P",
    "ASN": "N", "ASP": "D", "TYR": "Y", "PHE": "F", "TRP": "W",
}

# Where a target chain's numbering starts, for the 'renumbered' shape. 3, because that is
# where barnase starts in 1BRS -- the entry 037 pinned the renumbering problem against.
RENUMBER_FIRST_ID = 3
# The 'gapped' shape drops this many residues from the middle of the target's numbering.
GAP_MISSING_RESIDUES = 4
# And moves the second segment this far, so the junction is not a peptide bond. The same
# 12 A the 'broken' shape uses, and for the same reason: no peptide bond reaches that far.
GAP_DISPLACEMENT_A = 12.0


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


def rows_for(chain_id, entity, cycle, points, start_id, seq_ids=None):
    """mmCIF _atom_site rows for one chain, one CA per residue.

    seq_ids overrides label_seq_id, which is otherwise 1..n. The composition cycle is
    indexed by POSITION in the chain, not by label_seq_id, so renumbering a chain does
    not also change its sequence -- those are two different fixtures and conflating them
    would make either one unreadable.
    """
    rows, atom_id = [], start_id
    for i, point in enumerate(points, start=1):
        comp = cycle[(i - 1) % len(cycle)]
        label = str(seq_ids[i - 1]) if seq_ids else str(i)
        rows.append([
            "ATOM", str(atom_id), "C", "CA", ".", comp, chain_id, entity, label, "?",
            f"{point[0]:.3f}", f"{point[1]:.3f}", f"{point[2]:.3f}",
            "1.00", "20.00", label, chain_id, "1",
        ])
        atom_id += 1
    return rows, atom_id


def design_sequence(length, cycle=DESIGN_CYCLE):
    """The one-letter sequence of a design of this length. X where there is no code."""
    return "".join(ONE_LETTER.get(cycle[i % len(cycle)], "X") for i in range(length))


def design_cycle_for(shape):
    if shape == "nonstandard":
        return MIXED_CYCLE
    if shape == "allnonstandard":
        return NONSTANDARD_CYCLE
    return DESIGN_CYCLE


def write_design(path, shape, length, target_length):
    """Write one structure. Returns the target chain ids it used, in file order."""
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

    design_rows, next_id = rows_for(
        "B", "1", design_cycle_for(shape), designed, 1
    )

    target_blocks, target_chains = [], []
    if shape == "monomer":
        pass
    elif shape == "multichain":
        # Two target chains, A and C. The ids skip B -- which the design holds -- so a
        # consumer cannot pass by assuming the targets are the letters following it.
        half = max(target_length // 2, 1)
        first = lattice_walk(half, origin=offset)
        second_offset = (offset[0] + lattice_extent(half) + 6.0, offset[1], offset[2])
        second = lattice_walk(target_length - half, origin=second_offset)
        rows_a, next_id = rows_for("A", "0", TARGET_CYCLE, first, next_id)
        rows_c, next_id = rows_for("C", "2", TARGET_CYCLE, second, next_id)
        target_blocks, target_chains = [rows_a, rows_c], ["A", "C"]
    elif shape == "gapped":
        # One chain, two segments, with GAP_MISSING_RESIDUES absent from the numbering and
        # the second segment moved so the junction is not a bond.
        split = max(target_length // 2, 1)
        head = lattice_walk(split, origin=offset)
        tail_offset = (offset[0], offset[1] + GAP_DISPLACEMENT_A, offset[2])
        tail = lattice_walk(target_length - split, origin=tail_offset)
        ids = list(range(1, split + 1)) + [
            split + GAP_MISSING_RESIDUES + 1 + i for i in range(target_length - split)
        ]
        rows_a, next_id = rows_for(
            "A", "0", TARGET_CYCLE, head + tail, next_id, seq_ids=ids
        )
        target_blocks, target_chains = [rows_a], ["A"]
    elif shape == "renumbered":
        points = lattice_walk(target_length, origin=offset)
        ids = [RENUMBER_FIRST_ID + i for i in range(target_length)]
        rows_a, next_id = rows_for(
            "A", "0", TARGET_CYCLE, points, next_id, seq_ids=ids
        )
        target_blocks, target_chains = [rows_a], ["A"]
    else:
        points = lattice_walk(target_length, origin=offset)
        rows_a, next_id = rows_for("A", "0", TARGET_CYCLE, points, next_id)
        target_blocks, target_chains = [rows_a], ["A"]

    if shape == "swap":
        blocks = target_blocks + [design_rows]
    else:
        blocks = [design_rows] + target_blocks

    with open(path, "w") as handle:
        handle.write(f"data_{os.path.basename(path)}\n#\nloop_\n")
        for column in COLUMNS:
            handle.write(column + "\n")
        for block in blocks:
            for row in block:
                handle.write(" ".join(row) + "\n")
        handle.write("#\n")
    return target_chains


def main():
    argv = [a for a in sys.argv[1:] if a != "--manifest"]
    write_manifest = "--manifest" in sys.argv[1:]
    if len(argv) < 2:
        print(__doc__)
        return 2
    out_dir, count = argv[0], int(argv[1])
    shape = argv[2] if len(argv) > 2 else "compact"
    length = int(argv[3]) if len(argv) > 3 else 10
    target_length = max(length * 2, 20)
    os.makedirs(out_dir, exist_ok=True)
    target_chains = []
    for i in range(count):
        # Only the first design carries the shape; the rest are ordinary, so a case can
        # ask whether one bad design is ranked and reported correctly among good ones.
        this = shape if i == 0 else "compact"
        used = write_design(os.path.join(out_dir, f"design_{i}_generated.cif"),
                            this, length, target_length)
        if i == 0:
            target_chains = used
    if write_manifest:
        # The ids this file actually used, not a restatement of them. A 'monomer' run
        # correctly records none, which is what makes 'single-chain' the right answer
        # downstream rather than a fallback.
        with open(os.path.join(out_dir, "generation_manifest.json"), "w") as handle:
            json.dump({
                "target_chains": target_chains,
                "note": "written by tests/make_designs.py; synthetic, not a real run",
            }, handle, indent=2)
    print(f"wrote {count} design(s) to {out_dir}: shape={shape} length={length} "
          f"target_chains={target_chains or 'none'} "
          f"sequence={design_sequence(length, design_cycle_for(shape))}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
