#!/usr/bin/env python3
"""Node 04 — Analyze Designs.

Parse each *_generated.cif into a per-design metrics row and write design_metrics.csv /
design_metrics.json.

Generation runs with `~metrics`, so RFD4 emits no scores of its own — every number here is
computed from the geometry and composition of the output structures. None of it is a
confidence estimate: there is no refolding oracle in this workflow, so nothing here says a
design folds or binds. Self-consistency scoring is a follow-up that chains into
workflow-026.

Chain convention: the designed chain is taken to be the FIRST chain in the file, because
the contig built by node 02 puts the designed segment ahead of the chain break
("50-120, /0, A1-129"). Remaining chains are treated as target. All chains are reported
individually too, so a layout change is visible rather than silently mis-attributed.
"""
import csv
import glob
import json
import os
import sys

CONTACT_CUTOFF_A = 8.0  # CA-CA, the usual coarse interface definition

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
    "designed_length",
    "total_length",
    "designed_sequence",
    "hydrophobic_fraction",
    "net_charge",
    "radius_of_gyration_A",
    "n_interface_contacts",
    "n_target_residues_contacted",
    "min_interchain_distance_A",
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


def radius_of_gyration(coords):
    import numpy as np

    if len(coords) < 2:
        return None
    points = np.array(coords, dtype=float)
    centre = points.mean(axis=0)
    return float(np.sqrt(((points - centre) ** 2).sum(axis=1).mean()))


def interface(designed_coords, target_coords):
    """Contact count, distinct target residues contacted, and minimum distance."""
    import numpy as np

    if not designed_coords or not target_coords:
        return 0, 0, None
    a = np.array(designed_coords, dtype=float)
    b = np.array(target_coords, dtype=float)
    distances = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2)
    close = distances < CONTACT_CUTOFF_A
    return (
        int(close.sum()),
        int(close.any(axis=0).sum()),
        round(float(distances.min()), 2),
    )


def analyze(path):
    chains = read_chains(path)
    designed_id, designed = chains[0]
    targets = chains[1:]

    sequence = "".join(residue["one_letter"] for residue in designed)
    designed_coords = [r["ca"] for r in designed if r["ca"] is not None]
    target_coords = [
        r["ca"] for _, residues in targets for r in residues if r["ca"] is not None
    ]

    contacts, contacted, min_distance = interface(designed_coords, target_coords)
    gyration = radius_of_gyration(designed_coords)
    # X is an unrecognised residue; excluded from composition so it cannot dilute the ratio.
    known = [c for c in sequence if c != "X"]
    hydrophobic = sum(1 for c in known if c in HYDROPHOBIC)

    return {
        "design": os.path.basename(path),
        "n_chains": len(chains),
        "designed_chain": designed_id,
        "designed_length": len(designed),
        "total_length": sum(len(residues) for _, residues in chains),
        "designed_sequence": sequence,
        "hydrophobic_fraction": round(hydrophobic / len(known), 3) if known else None,
        "net_charge": sum(1 for c in sequence if c in POSITIVE)
        - sum(1 for c in sequence if c in NEGATIVE),
        "radius_of_gyration_A": round(gyration, 2) if gyration is not None else None,
        "n_interface_contacts": contacts,
        "n_target_residues_contacted": contacted,
        "min_interchain_distance_A": min_distance,
        "chain_lengths": ";".join(f"{cid}:{len(res)}" for cid, res in chains),
    }


def main():
    os.makedirs("./outputs", exist_ok=True)

    designs = sorted(glob.glob("./inputs/*_generated.cif"))
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

    rows = [analyze(path) for path in designs]

    with open("./outputs/design_metrics.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    payload = {
        "n_designs": len(rows),
        "contact_cutoff_angstrom": CONTACT_CUTOFF_A,
        "designed_chain_convention": "first chain in file",
        "metrics_are_geometric_only": True,
        "seed": manifest.get("seed"),
        "hf_ckpt_revision": manifest.get("hf_ckpt_revision"),
        "rfd4_proteina_version": manifest.get("rfd4_proteina_version"),
        "securedna_screened": False,
        "designs": rows,
    }
    with open("./outputs/design_metrics.json", "w") as handle:
        json.dump(payload, handle, indent=2)

    for row in rows:
        print(
            f"  {row['design']}: designed chain {row['designed_chain']} "
            f"len={row['designed_length']}, contacts={row['n_interface_contacts']}, "
            f"Rg={row['radius_of_gyration_A']} A",
            flush=True,
        )
    print("Wrote ./outputs/design_metrics.csv and ./outputs/design_metrics.json", flush=True)


if __name__ == "__main__":
    main()
