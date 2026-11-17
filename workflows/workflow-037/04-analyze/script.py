#!/usr/bin/env python3
"""Node 04 — Analyze Designs.

Parse each *_generated.cif into a per-design metrics row and write design_metrics.csv /
design_metrics.json.

Generation runs with `~metrics`, so RFD4 emits no scores of its own — every number here is
computed from the geometry and composition of the output structures. None of it is a
confidence estimate: there is no refolding oracle in this workflow, so nothing here says a
design folds or binds. Self-consistency scoring is a follow-up that chains into
workflow-026.

Interface convention: contacts here are CA-CA pairs within 8.0 A. That is deliberately NOT
node 02's epitope definition, which is any heavy atom within 4.5 A -- upstream's trained
hotspot criterion. The two measure different things and should not be aligned. Node 02 has
to reproduce the convention the model was conditioned with, on a known complex, because a
wrong epitope aims the design at the wrong patch. Node 04 is describing a structure the
model has just produced, where a coarse count of how much of the target the design sits
against is the useful summary, and alpha carbons give that cheaply without pretending to
resolve individual atom contacts in a backbone that was generated, not measured.

So the number in the report's "Contacts" column is not the number of hotspot atoms the
design was aimed with, and it is not meant to be.

Chain convention: the designed chain is the one that was NOT in the target. Node 03
records the target's chain ids in generation_manifest.json, and the chain of a generated
file whose id is not among them is the design; the rest are target.

This node used to take the first chain in the file. That is true of what node 02's contig
produces ("50-120, /0, A1-129"), but it is a fact about node 02 being asserted two nodes
away from where it is decided, and when it did not hold -- chains written target-first --
every metric was attributed to the target without a word. Reading the answer costs one
field in the manifest and removes the assumption entirely.

Where the manifest is missing, or the ids do not separate the chains cleanly, the first
chain is used, the log says so, and the row records designed_chain_source, so a fallback
is visible rather than indistinguishable from a resolved answer. All chains are reported
individually in chain_lengths either way.
"""
import csv
import glob
import json
import os
import re
import sys

CONTACT_CUTOFF_A = 8.0  # CA-CA. See the note below on why this is not node 02's 4.5.

# Below this, two chains are interpenetrating rather than touching. The anchor is the
# CA-CA virtual bond: consecutive, covalently bonded residues in a trans peptide sit
# 3.80-3.81 A apart, and that is as close as two alpha carbons get in real polypeptide.
# Two CAs on *different* chains closer than a covalent CA-CA bond is not a tight
# interface, it is backbone passing through backbone. Generation is free to produce it,
# nothing upstream rejects it, and it scores *better* on a raw contact count than a good
# design does, so it has to be caught here or it reaches the top of the report.
CLASH_CA_DISTANCE_A = 3.8

HYDROPHOBIC = set("AVILMFWYC")
POSITIVE = set("KR")
NEGATIVE = set("DE")

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}

COLUMNS = [
    "design",
    "n_chains",
    "designed_chain",
    "designed_chain_source",
    "designed_length",
    "total_length",
    "designed_sequence",
    "hydrophobic_fraction",
    "net_charge",
    "radius_of_gyration_A",
    "expected_radius_of_gyration_A",
    "compactness_ratio",
    "overextended",
    "chain_breaks",
    "n_interface_contacts",
    "n_target_residues_contacted",
    "min_interchain_ca_distance_A",
    "steric_clash",
    "chain_lengths",
]


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def read_chains(path):
    """Return [(chain_id, [(one_letter, ca_xyz_or_None), ...]), ...] in file order."""
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    try:
        cif = MMCIF2Dict(path)
    except Exception as exc:
        fail(f"Could not parse '{path}': {exc}")

    if "_atom_site.label_asym_id" not in cif:
        fail(f"'{path}' has no _atom_site records; not a structure file.")

    n_atoms = len(cif["_atom_site.label_asym_id"])

    def column(key, default="."):
        return cif.get(key, [default] * n_atoms)

    asym = column("_atom_site.label_asym_id")
    comp = column("_atom_site.label_comp_id")
    seq_id = column("_atom_site.label_seq_id")
    atom_id = column("_atom_site.label_atom_id")
    xs = column("_atom_site.Cartn_x", "0")
    ys = column("_atom_site.Cartn_y", "0")
    zs = column("_atom_site.Cartn_z", "0")

    order = []
    residues = {}
    for i in range(n_atoms):
        if seq_id[i] in (".", "?", ""):
            continue
        chain = asym[i]
        if chain not in residues:
            residues[chain] = {}
            order.append(chain)
        key = int(seq_id[i])
        entry = residues[chain].setdefault(
            key, {"one_letter": THREE_TO_ONE.get(comp[i], "X"), "ca": None}
        )
        if atom_id[i] == "CA":
            try:
                entry["ca"] = (float(xs[i]), float(ys[i]), float(zs[i]))
            except ValueError:
                entry["ca"] = None

    chains = []
    for chain in order:
        items = [residues[chain][k] for k in sorted(residues[chain])]
        chains.append((chain, items))
    if not chains:
        fail(f"'{path}' contains no polymer residues.")
    return chains


# A trans peptide bond puts consecutive alpha carbons 3.80-3.81 A apart; cis prolines
# and ordinary distortion stay well inside this window. Outside it there is no bond.
CA_BOND_MIN_A = 2.8
CA_BOND_MAX_A = 4.8

# A compact folded protein has Rg ~ 2.2 * N^0.38 A; a chemically denatured one follows
# Rg ~ 2.54 * N^0.522 A (Kohn et al. 2004). The second over the first is 2.06 at N=60 and
# 2.35 at N=150, so an unfolded chain sits at roughly twice the compact value whatever
# its length, and the ratio separates a fold from a coil without knowing the fold.
#
# The limit is measured, not chosen. Nineteen designs from six archived runs of this
# workflow ranged 0.94 to 1.49, the loose end being monomer and motif_scaffold runs,
# which produce more elongated folds than binder runs do. 1.8 clears the worst real
# design by 0.31 and sits 0.26 below the tightest coil. A straight chain, for reference,
# comes out at 6.3.
#
# This is a broken-model check, not a quality check: it catches output that is not a
# folded protein at all. It will not notice a design that folds and does not bind.
COMPACT_COEFFICIENT = 2.2
COMPACT_EXPONENT = 0.38
COMPACTNESS_LIMIT = 1.8
MIN_LENGTH_FOR_COMPACTNESS = 30


def count_chain_breaks(coords):
    """Consecutive alpha carbons that no peptide bond could join.

    Computed over the CA positions in sequence order. A design arrives as one continuous
    chain, so any break means the model emitted something that is not a polypeptide.
    """
    import numpy as np

    if len(coords) < 2:
        return 0
    points = np.array(coords, dtype=float)
    distances = np.sqrt(((points[1:] - points[:-1]) ** 2).sum(axis=1))
    return int(((distances < CA_BOND_MIN_A) | (distances > CA_BOND_MAX_A)).sum())


def expected_radius_of_gyration(n_residues):
    if n_residues < 1:
        return None
    return COMPACT_COEFFICIENT * (n_residues ** COMPACT_EXPONENT)


def radius_of_gyration(coords):
    import numpy as np

    if len(coords) < 2:
        return None
    points = np.array(coords, dtype=float)
    centre = points.mean(axis=0)
    return float(np.sqrt(((points - centre) ** 2).sum(axis=1).mean()))


def interface(designed_coords, target_coords):
    """Contact count, distinct target residues contacted, and minimum CA-CA distance.

    With no second chain there is no interface at all, so the three columns are None
    rather than 0: a monomer run reporting "0 contacts" reads as "this design failed to
    engage the target", which is a judgement about a target that does not exist.
    """
    import numpy as np

    if not designed_coords or not target_coords:
        return None, None, None
    a = np.array(designed_coords, dtype=float)
    b = np.array(target_coords, dtype=float)
    distances = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    close = distances < CONTACT_CUTOFF_A
    return (
        int(close.sum()),
        int(close.any(axis=0).sum()),
        round(float(distances.min()), 2),
    )


def resolve_designed_chain(chains, target_chains, path):
    """Pick the designed chain, and say how it was picked.

    Returns (index, source). With the target's chain ids in hand this is a set
    difference, not a guess. Two cases decline to answer: no ids recorded, and more than
    one chain left over after removing the target's -- a multi-segment design, or chains
    relabelled by upstream. Choosing silently in either case would be worse than the
    first-chain fallback plus a warning, because the answer would look resolved.
    """
    if target_chains:
        known = set(target_chains)
        candidates = [i for i, (cid, _) in enumerate(chains) if cid not in known]
        if len(candidates) == 1:
            return candidates[0], "manifest-target-chains"
        if len(chains) == 1:
            # monomer and motif_scaffold runs come back as the design alone. Even if its
            # id collides with a target chain id, there is only one chain to report.
            return 0, "single-chain"
        print(
            f"WARNING: '{os.path.basename(path)}' has chains "
            f"{[cid for cid, _ in chains]} and this run's target chains were "
            f"{sorted(known)}, which leaves {len(candidates)} candidate(s) for the "
            "design. Falling back to the first chain in the file; check chain_lengths.",
            flush=True,
        )
    elif len(chains) > 1:
        print(
            "WARNING: generation_manifest.json records no target chains, so the design "
            f"in '{os.path.basename(path)}' is being taken as the first chain in the "
            "file. That holds for the contigs this workflow builds, and is not checked "
            "here.",
            flush=True,
        )
    return 0, "assumed-first-chain"


def analyze(path, target_chains=None):
    chains = read_chains(path)
    index, chain_source = resolve_designed_chain(chains, target_chains, path)
    designed_id, designed = chains[index]
    targets = [pair for i, pair in enumerate(chains) if i != index]

    sequence = "".join(residue["one_letter"] for residue in designed)
    designed_coords = [r["ca"] for r in designed if r["ca"] is not None]
    target_coords = [
        r["ca"] for _, residues in targets for r in residues if r["ca"] is not None
    ]

    contacts, contacted, min_distance = interface(designed_coords, target_coords)
    gyration = radius_of_gyration(designed_coords)
    breaks = count_chain_breaks(designed_coords)
    expected_gyration = expected_radius_of_gyration(len(designed))
    ratio = (
        round(gyration / expected_gyration, 2)
        if gyration is not None and expected_gyration
        else None
    )
    # Only a verdict where the relation means something; the number is still reported.
    overextended = (
        ratio is not None
        and len(designed) >= MIN_LENGTH_FOR_COMPACTNESS
        and ratio > COMPACTNESS_LIMIT
    )
    if len(chains) > 2:
        print(
            f"WARNING: '{os.path.basename(path)}' has {len(chains)} chains. This node "
            f"treats {designed_id} as the design and the rest as target, so with several "
            "designed segments the extra ones are counted as target and the interface "
            "numbers mix designed-to-designed contacts in. Check chain_lengths.",
            flush=True,
        )
    # X is an unrecognised residue; excluded from composition so it cannot dilute the ratio.
    known = [c for c in sequence if c != "X"]
    hydrophobic = sum(1 for c in known if c in HYDROPHOBIC)

    return {
        "design": os.path.basename(path),
        "n_chains": len(chains),
        "designed_chain": designed_id,
        "designed_chain_source": chain_source,
        "designed_length": len(designed),
        "total_length": sum(len(residues) for _, residues in chains),
        "designed_sequence": sequence,
        "hydrophobic_fraction": round(hydrophobic / len(known), 3) if known else None,
        "net_charge": sum(1 for c in sequence if c in POSITIVE)
        - sum(1 for c in sequence if c in NEGATIVE),
        "radius_of_gyration_A": round(gyration, 2) if gyration is not None else None,
        "expected_radius_of_gyration_A": (
            round(expected_gyration, 2) if expected_gyration is not None else None
        ),
        "compactness_ratio": ratio,
        # The verdict, not just the number. Node 05 renders this rather than applying
        # COMPACTNESS_LIMIT itself: the threshold and the minimum length it is only
        # meaningful above are this node's to decide, and a report that re-derived them
        # could disagree with the warning printed here about the same design.
        "overextended": overextended,
        "chain_breaks": breaks,
        "n_interface_contacts": contacts,
        "n_target_residues_contacted": contacted,
        "min_interchain_ca_distance_A": min_distance,
        "steric_clash": (
            min_distance is not None and min_distance < CLASH_CA_DISTANCE_A
        ),
        "chain_lengths": ";".join(f"{cid}:{len(res)}" for cid, res in chains),
    }


def chain_convention(rows):
    """One phrase for the report's run info, naming how the designed chain was chosen.

    Spelled out per source rather than summarised, because a batch where some designs
    resolved and others fell back would otherwise be reported as whichever reads better.
    """
    wording = {
        "manifest-target-chains": (
            "the chain absent from the run's target, read from the generation manifest"
        ),
        "single-chain": "the only chain in the file",
        "assumed-first-chain": (
            "first chain in file, assumed: no target chain ids to check against"
        ),
    }
    sources = sorted({row.get("designed_chain_source") for row in rows})
    return "; ".join(wording.get(source, str(source)) for source in sources)


def natural_key(name):
    """Sort key that reads runs of digits as numbers.

    RFD4 names designs with a trailing index, so a plain lexicographic sort puts
    design_10 between design_1 and design_2. Every node downstream keys off the
    design name, so the order only ever affected what a reader sees -- but what a
    reader sees is the whole point of the next two nodes.
    """
    return [int(p) if p.isdigit() else p for p in re.split(r"(\d+)", name)]


def main():
    os.makedirs("./outputs", exist_ok=True)

    designs = sorted(glob.glob("./inputs/*_generated.cif"),
                     key=lambda p: natural_key(os.path.basename(p)))
    if not designs:
        fail("No generated structures found in ./inputs/ (*_generated.cif)")

    manifest = {}
    manifest_path = "./inputs/generation_manifest.json"
    if os.path.exists(manifest_path):
        try:
            with open(manifest_path) as handle:
                manifest = json.load(handle)
        except Exception as exc:
            print(f"WARNING: could not read generation manifest: {exc}", flush=True)

    print(f"Found {len(designs)} generated design(s)", flush=True)

    target_chains = manifest.get("target_chains")
    if target_chains:
        print(f"Target chains recorded by node 03: {target_chains}", flush=True)

    rows = [analyze(path, target_chains) for path in designs]

    with open("./outputs/design_metrics.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    payload = {
        "n_designs": len(rows),
        "contact_cutoff_angstrom": CONTACT_CUTOFF_A,
        "designed_chain_convention": chain_convention(rows),
        "metrics_are_geometric_only": True,
        "seed": manifest.get("seed"),
        "hf_ckpt_revision": manifest.get("hf_ckpt_revision"),
        "rfd4_proteina_version": manifest.get("rfd4_proteina_version"),
        "securedna_screened": False,
        "designs": rows,
    }
    with open("./outputs/design_metrics.json", "w") as handle:
        json.dump(payload, handle, indent=2)

    broken = [r for r in rows if r["chain_breaks"]]
    if broken:
        named = ", ".join(f"{r['design']} ({r['chain_breaks']})" for r in broken)
        print(
            f"WARNING: {len(broken)} design(s) contain a gap between consecutive alpha "
            f"carbons that no peptide bond spans: {named}. A design arrives as one "
            "continuous chain, so this is not a borderline geometry -- it means the "
            "model emitted something that is not a polypeptide. Check the checkpoint "
            "and the generation log before reading any other number here.",
            flush=True,
        )
    loose = [r for r in rows if r["overextended"]]
    if loose:
        named = ", ".join(
            f"{r['design']} ({r['radius_of_gyration_A']} A vs "
            f"{r['expected_radius_of_gyration_A']} A expected)" for r in loose
        )
        print(
            f"WARNING: {len(loose)} design(s) are far more extended than a folded protein "
            f"of their length: {named}. Rg should be near 2.2*N^0.38 A for a compact "
            "fold; several times that is a chain that has not folded. This is what a "
            "wrong checkpoint or a failed generation looks like, and no other metric "
            "here would show it -- contacts, charge and composition are all still "
            "computable on an unfolded chain.",
            flush=True,
        )
    for row in rows:
        print(
            f"  {row['design']}: designed chain {row['designed_chain']} "
            f"len={row['designed_length']}, contacts={row['n_interface_contacts']}, "
            f"Rg={row['radius_of_gyration_A']} A"
            + (" [STERIC CLASH]" if row["steric_clash"] else ""),
            flush=True,
        )
    clashing = [row["design"] for row in rows if row["steric_clash"]]
    if clashing:
        print(
            f"WARNING: {len(clashing)} design(s) have an inter-chain CA-CA distance under "
            f"{CLASH_CA_DISTANCE_A} A, which is shorter than the CA-CA virtual bond between "
            "two bonded residues: the chains interpenetrate rather than dock. "
            f"({', '.join(clashing)}) A clash inflates the contact count rather than "
            "lowering it, so these would otherwise read as the best designs in the batch.",
            flush=True,
        )
    print("Wrote ./outputs/design_metrics.csv and ./outputs/design_metrics.json", flush=True)


if __name__ == "__main__":
    main()
