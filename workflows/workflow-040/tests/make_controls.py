#!/usr/bin/env python3
"""Turn a deposited structure into something 040 can ingest, for the Tier 2 controls.

040 consumes workflow-037's output: `*_generated.cif` plus a generation manifest. A
native complex from the PDB is not in that shape, so the controls need a converter, and
this is it. What comes out is a CA-only mmCIF in exactly the column order
`workflow-037/tests/make_designs.py` writes, which is the format both workflows' nodes
already read.

CA-only on purpose, not as a shortcut: node 040/04 superposes CA, node 037/04 counts
contacts between CA, and nothing in either workflow reads another backbone atom. Carrying
side chains would mean this file had to decide what to do with altlocs, partial occupancy
and missing atoms on a path where nobody would ever look at the answer.

The four Phase 1 controls, all from the same two entries:

    C1  native sequence on its own backbone            must PASS
    C2  a second native complex, numbered from 3       must PASS
    C3  C1's backbone, sequence shuffled               must FAIL
    C4  C1's backbone, sequence all alanine            must FAIL

C3 and C4 keep the backbone and change the sequence, which is the only way a
self-consistency score can be wrong about them: the backbone is a real folded protein, so
a scorer that rewards plausible geometry rather than sequence-structure agreement passes
all four and the control set catches it.

Chains are chosen by auth_asym_id
---------------------------------
Because "chain A" means the auth id to everyone who reads a PDB entry, while mmCIF's
label_asym_id splits an entry into one id per entity instance -- polymers, ligands and
each water group alike. 037 lost five of 1BRS's six chains to exactly that confusion
(53d1d460). Both ids are recorded for every chain that comes through, and --list prints
them so a control is chosen by measurement rather than from memory, which is the rule in
tests/README.md.

Numbering is kept, not reset
----------------------------
label_seq_id comes through from the source. Barnase in 1BRS starts at 3, and that is the
single most useful property in the control set: a correspondence built by matching residue
numbers instead of positions is wrong there and nowhere else (ef74646e). Resetting it
would throw away the one control that can catch it.

Run
---
  python3 -I tests/make_controls.py --list <entry.cif>
  python3 -I tests/make_controls.py <entry.cif> <out-dir> \\
      --name design_0 --binder <auth-chain> --target <auth-chain>[,<auth-chain>...] \\
      [--sequence native|scrambled|poly-alanine] [--seed 0]

-I because a file downloaded from the PDB is untrusted input and must not be able to put
a module on the import path. Nothing here imports anything but the standard library, so
-I costs nothing.
"""
import hashlib
import json
import os
import random
import sys

# Byte-for-byte the column order make_designs.py writes, so a node cannot tell a converted
# entry from a generated design by its header.
COLUMNS = [
    "_atom_site.group_PDB", "_atom_site.id", "_atom_site.type_symbol",
    "_atom_site.label_atom_id", "_atom_site.label_alt_id", "_atom_site.label_comp_id",
    "_atom_site.label_asym_id", "_atom_site.label_entity_id", "_atom_site.label_seq_id",
    "_atom_site.pdbx_PDB_ins_code", "_atom_site.Cartn_x", "_atom_site.Cartn_y",
    "_atom_site.Cartn_z", "_atom_site.occupancy", "_atom_site.B_iso_or_equiv",
    "_atom_site.auth_seq_id", "_atom_site.auth_asym_id", "_atom_site.pdbx_PDB_model_num",
]

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}

# Where the binder and the targets land in the output. Binder B and target A, matching
# make_designs.py's compact shape: the manifest then records the targets, the design chain
# is the one absent from that list, and node 01 resolves it down the good path rather than
# the first-chain fallback -- which is the behaviour C1 is partly there to confirm.
BINDER_LABEL = "B"
TARGET_LABELS = ["A", "C", "D", "E", "F", "G", "H"]


def die(message):
    print(f"make_controls: ERROR: {message}", file=sys.stderr)
    sys.exit(1)


def read_atom_site(path):
    """({column: index}, [row fields]) for the _atom_site loop. Hand-parsed.

    Handles the one-row-per-line form RCSB serves. Quoted fields containing spaces are
    rejected rather than mis-split: no column this file reads is ever quoted in an RCSB
    _atom_site loop, so a row that needs quote handling is a row this file has not been
    shown and should not guess at.
    """
    columns, rows, state = [], [], "search"
    with open(path, errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            stripped = line.strip()
            if state == "search":
                if stripped.startswith("_atom_site."):
                    columns = [stripped]
                    state = "head"
                continue
            if state == "head":
                if stripped.startswith("_atom_site."):
                    columns.append(stripped)
                    continue
                state = "rows"
            if state == "rows":
                if not stripped or stripped.startswith(("#", "loop_", "_", "data_")):
                    break
                if '"' in stripped or "'" in stripped:
                    continue
                fields = stripped.split()
                if len(fields) == len(columns):
                    rows.append(fields)
    if not rows:
        die(f"no usable _atom_site rows in {path}")
    return {name.split(".", 1)[1]: i for i, name in enumerate(columns)}, rows


def chains_in(path):
    """Per auth chain, in file order: the CA residues this file would take from it.

    One CA per residue. HETATM rows are kept when they carry a numeric label_seq_id,
    because that is what a modified polymer residue looks like -- 037's 5O45 fixture has
    106 of them -- and dropped otherwise, which is what removes waters and ligands.
    Altloc is taken as '.' or the first letter, and only model 1.
    """
    index, rows = read_atom_site(path)

    def field(row, name, default="."):
        at = index.get(name)
        return row[at] if at is not None else default

    order, chains = [], {}
    skipped = {"no_seq_id": 0, "not_ca": 0, "altloc": 0, "model": 0}
    for row in rows:
        if field(row, "pdbx_PDB_model_num", "1") not in ("1", "."):
            skipped["model"] += 1
            continue
        if field(row, "label_atom_id") != "CA":
            skipped["not_ca"] += 1
            continue
        if field(row, "label_alt_id") not in (".", "?", "A"):
            skipped["altloc"] += 1
            continue
        raw_seq = field(row, "label_seq_id")
        try:
            seq_id = int(raw_seq)
        except ValueError:
            skipped["no_seq_id"] += 1
            continue
        auth = field(row, "auth_asym_id", field(row, "label_asym_id"))
        if auth not in chains:
            chains[auth] = {"auth": auth, "label": set(), "residues": [], "seen": set()}
            order.append(auth)
        entry = chains[auth]
        if seq_id in entry["seen"]:
            continue
        entry["seen"].add(seq_id)
        entry["label"].add(field(row, "label_asym_id"))
        entry["residues"].append({
            "seq_id": seq_id,
            "auth_seq_id": field(row, "auth_seq_id", raw_seq),
            "comp": field(row, "label_comp_id"),
            "group": field(row, "group_PDB", "ATOM"),
            "xyz": (
                float(field(row, "Cartn_x", "0")),
                float(field(row, "Cartn_y", "0")),
                float(field(row, "Cartn_z", "0")),
            ),
        })
    out = []
    for auth in order:
        entry = chains[auth]
        entry["residues"].sort(key=lambda r: r["seq_id"])
        ids = [r["seq_id"] for r in entry["residues"]]
        entry["label"] = sorted(entry["label"])
        entry.pop("seen")
        entry["n_residues"] = len(ids)
        entry["first_id"] = ids[0] if ids else None
        entry["last_id"] = ids[-1] if ids else None
        entry["gaps"] = [
            {"after": a, "before": b, "missing": b - a - 1}
            for a, b in zip(ids, ids[1:]) if b - a > 1
        ]
        entry["sequence"] = "".join(
            THREE_TO_ONE.get(r["comp"], "X") for r in entry["residues"]
        )
        entry["n_unknown"] = entry["sequence"].count("X")
        entry["n_hetatm"] = len([r for r in entry["residues"] if r["group"] != "ATOM"])
        out.append(entry)
    return out, skipped


def transform(residues, how, seed):
    """Return (residues, note). The BACKBONE never moves; only the residue identities do."""
    if how == "native":
        return residues, "the deposited sequence, unchanged"
    if how == "poly-alanine":
        out = [dict(r, comp="ALA") for r in residues]
        return out, f"every one of {len(residues)} residues set to ALA"
    if how == "scrambled":
        # A permutation of the chain's own residues, so composition is preserved exactly
        # and only the ORDER changes. That is the point of having C3 as well as C4: a
        # scorer could be fooled by a sequence with the right composition and not by one
        # with the wrong composition, and the two controls separate those.
        comps = [r["comp"] for r in residues]
        shuffled = list(comps)
        rng = random.Random(seed)
        for _ in range(1000):
            rng.shuffle(shuffled)
            if shuffled != comps:
                break
        else:
            die("could not find a permutation that differs from the original")
        moved = sum(1 for a, b in zip(comps, shuffled) if a != b)
        out = [dict(r, comp=c) for r, c in zip(residues, shuffled)]
        return out, (
            f"the same {len(comps)} residues permuted with random.Random({seed}); "
            f"{moved} of {len(comps)} positions changed identity, composition identical"
        )
    die(f"unknown --sequence '{how}'")


def write_cif(path, blocks):
    """blocks: [(label_asym_id, auth_asym_id, entity, residues)]. One CA per residue."""
    atom_id = 1
    with open(path, "w") as handle:
        handle.write(f"data_{os.path.basename(path)}\n#\nloop_\n")
        for column in COLUMNS:
            handle.write(column + "\n")
        for label, auth, entity, residues in blocks:
            for residue in residues:
                x, y, z = residue["xyz"]
                handle.write(" ".join([
                    "ATOM", str(atom_id), "C", "CA", ".", residue["comp"],
                    label, str(entity), str(residue["seq_id"]), "?",
                    f"{x:.3f}", f"{y:.3f}", f"{z:.3f}", "1.00", "20.00",
                    str(residue["auth_seq_id"]), label, "1",
                ]) + "\n")
                atom_id += 1
        handle.write("#\n")
    return atom_id - 1


def parse_args(argv):
    options = {
        "name": "design_0", "binder": None, "target": None,
        "sequence": "native", "seed": "0",
    }
    positional, i = [], 0
    while i < len(argv):
        token = argv[i]
        if token.startswith("--"):
            key = token[2:]
            if key not in options:
                die(f"unknown option '{token}'")
            if i + 1 >= len(argv):
                die(f"'{token}' needs a value")
            options[key] = argv[i + 1]
            i += 2
        else:
            positional.append(token)
            i += 1
    return positional, options


def main():
    argv = sys.argv[1:]
    if not argv:
        print(__doc__)
        return 2

    if argv[0] == "--list":
        if len(argv) < 2:
            die("--list needs a path")
        chains, skipped = chains_in(argv[1])
        print(f"{argv[1]}  sha256={sha256(argv[1])[:16]}")
        print(f"  rows not taken: {skipped}")
        for chain in chains:
            gaps = ", ".join(
                "%d->%d (%d missing)" % (g["after"], g["before"], g["missing"])
                for g in chain["gaps"]
            ) or "none"
            print(f"  auth {chain['auth']:<3} label {','.join(chain['label']):<12} "
                  f"n={chain['n_residues']:<5} ids {chain['first_id']}..{chain['last_id']}"
                  f"  X={chain['n_unknown']}  HETATM={chain['n_hetatm']}  gaps: {gaps}")
            print(f"       {chain['sequence'][:72]}"
                  + ("..." if len(chain["sequence"]) > 72 else ""))
        return 0

    positional, options = parse_args(argv)
    if len(positional) != 2:
        die("needs <entry.cif> <out-dir>, or --list <entry.cif>")
    source, out_dir = positional
    if not options["binder"] or not options["target"]:
        die("--binder and --target are required: name them from --list output, not from "
            "memory. 'chain A' is not a sufficient description of an entry (53d1d460)")

    chains, skipped = chains_in(source)
    by_auth = {c["auth"]: c for c in chains}
    binder_id = options["binder"]
    target_ids = [t.strip() for t in options["target"].split(",") if t.strip()]
    for wanted in [binder_id] + target_ids:
        if wanted not in by_auth:
            die(f"auth chain '{wanted}' is not in {source}. It has "
                f"{sorted(by_auth)} -- run --list and choose from those")
    if binder_id in target_ids:
        die(f"'{binder_id}' is named as both the binder and a target")
    if len(target_ids) > len(TARGET_LABELS):
        die(f"{len(target_ids)} targets is more than this file relabels")

    binder = by_auth[binder_id]
    if binder["n_residues"] < 4:
        die(f"the binder chain has {binder['n_residues']} residues; a superposition over "
            "fewer than 4 is determined by its own inputs")

    residues, note = transform(binder["residues"], options["sequence"],
                               int(options["seed"]))
    before = binder["sequence"]
    after = "".join(THREE_TO_ONE.get(r["comp"], "X") for r in residues)
    if options["sequence"] != "native" and after == before:
        die(f"--sequence {options['sequence']} left the sequence unchanged, so this "
            "control could not fail for the reason it exists to test")
    # The backbone must be untouched: that is the entire design of C3 and C4, and it is
    # checked here rather than assumed, because a converter that moved an atom would make
    # a negative control fail for the wrong reason.
    for original, moved in zip(binder["residues"], residues):
        if original["xyz"] != moved["xyz"] or original["seq_id"] != moved["seq_id"]:
            die("the transform moved a backbone atom or changed a residue number")

    os.makedirs(out_dir, exist_ok=True)
    blocks = [(BINDER_LABEL, binder_id, 1, residues)]
    target_map = {}
    for i, auth in enumerate(target_ids):
        label = TARGET_LABELS[i]
        target_map[auth] = label
        blocks.append((label, auth, i + 2, by_auth[auth]["residues"]))

    cif_path = os.path.join(out_dir, f"{options['name']}_generated.cif")
    n_atoms = write_cif(cif_path, blocks)

    # The ids this run actually used, which is what 037's node 04 reads to decide which
    # chain is the design. Not a restatement of the request: the labels, as written.
    manifest_path = os.path.join(out_dir, "generation_manifest.json")
    existing = {}
    if os.path.exists(manifest_path):
        with open(manifest_path) as handle:
            existing = json.load(handle)
    recorded = sorted(set((existing.get("target_chains") or [])
                          + [target_map[a] for a in target_ids]))
    with open(manifest_path, "w") as handle:
        json.dump({
            "target_chains": recorded,
            "note": ("written by tests/make_controls.py from deposited structures; the "
                     "designed chain is the one absent from this list"),
        }, handle, indent=2)

    provenance = {
        "name": options["name"],
        "source_file": os.path.basename(source),
        "source_sha256": sha256(source),
        "rows_not_taken": skipped,
        "sequence_transform": options["sequence"],
        "transform_seed": int(options["seed"]),
        "transform_note": note,
        "binder": {
            "auth_asym_id": binder_id,
            "label_asym_id_in_source": binder["label"],
            "label_asym_id_written": BINDER_LABEL,
            "n_residues": binder["n_residues"],
            "first_id": binder["first_id"],
            "last_id": binder["last_id"],
            "gaps": binder["gaps"],
            "n_unknown": binder["n_unknown"],
            "n_hetatm": binder["n_hetatm"],
            "sequence_deposited": before,
            "sequence_written": after,
        },
        "targets": [{
            "auth_asym_id": auth,
            "label_asym_id_in_source": by_auth[auth]["label"],
            "label_asym_id_written": target_map[auth],
            "n_residues": by_auth[auth]["n_residues"],
            "first_id": by_auth[auth]["first_id"],
            "last_id": by_auth[auth]["last_id"],
            "gaps": by_auth[auth]["gaps"],
            "n_unknown": by_auth[auth]["n_unknown"],
            "sequence": by_auth[auth]["sequence"],
        } for auth in target_ids],
        "wrote": {"cif": os.path.basename(cif_path), "n_atoms": n_atoms},
    }
    prov_path = os.path.join(out_dir, f"{options['name']}_provenance.json")
    with open(prov_path, "w") as handle:
        json.dump(provenance, handle, indent=2)

    print(f"make_controls: {options['name']} from {os.path.basename(source)}")
    print(f"  binder  auth {binder_id} -> {BINDER_LABEL}  {binder['n_residues']} residues"
          f"  ids {binder['first_id']}..{binder['last_id']}"
          f"  gaps {len(binder['gaps'])}")
    print(f"  targets " + ", ".join(
        f"auth {a} -> {target_map[a]} ({by_auth[a]['n_residues']} residues, ids "
        f"{by_auth[a]['first_id']}..{by_auth[a]['last_id']})" for a in target_ids
    ))
    print(f"  sequence: {note}")
    if options["sequence"] != "native":
        print(f"    deposited {before[:60]}")
        print(f"    written   {after[:60]}")
    print(f"  wrote {cif_path} ({n_atoms} CA atoms) and {prov_path}")
    print(f"  manifest target_chains={recorded}")
    return 0


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


if __name__ == "__main__":
    sys.exit(main())
