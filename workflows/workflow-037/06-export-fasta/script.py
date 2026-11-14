#!/usr/bin/env python3
"""Node 06 — Export FASTA.

Write the designed sequences as `sequences.fasta`, which is exactly the filename and format
workflow-026 (ColabFold) takes as its input, so refolding a batch of designs is a file copy
rather than a reformatting job.

Why this node exists: this workflow runs with `~metrics` and scores nothing, so a design's
only real test is independent refolding — fold the sequence, check the fold matches the
backbone RFD4 drew (self-consistency), and keep the designs that survive. That filter lives
in workflow-026; this node is the handoff.

Each record holds one design's generated chain alone. Refolding that answers "does this
sequence fold to the shape it was designed as?", which is the standard self-consistency
check and the only one this handoff can support.

There is deliberately no complex-export mode. Writing every chain colon-joined is the
only encoding that would make ColabFold predict the complex, because workflow-026 hands
the FASTA directory straight to colabfold_batch and each record becomes an independent
job — a chain break signalled in the header would be ignored. But workflow-026 validates
the sequence against ACDEFGHIKLMNPQRSTVWYX and rejects anything else, so a colon-joined
record fails its first node outright:

    ERROR: Sequence 'design_0_generated': invalid characters [':']

The two requirements are strictly incompatible from this side, and this workflow must not
change workflow-026 to suit itself. Re-enabling complex export means adding ":" to that
validator's VALID_AA first.
"""
import glob
import json
import os
import re
import sys

LINE_WIDTH = 60
OUTPUT_NAME = "sequences.fasta"

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def chain_sequences(path):
    """Return [(chain_id, sequence), ...] in the order the chains appear in the file."""
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    try:
        cif = MMCIF2Dict(path)
    except Exception as exc:
        fail(f"Could not parse '{path}': {exc}")

    if "_atom_site.label_asym_id" not in cif:
        fail(f"'{path}' has no _atom_site records; not a structure file.")

    asym = cif["_atom_site.label_asym_id"]
    seq_id = cif.get("_atom_site.label_seq_id", ["."] * len(asym))
    comp = cif.get("_atom_site.label_comp_id", ["UNK"] * len(asym))

    order = []
    residues = {}
    for i in range(len(asym)):
        if seq_id[i] in (".", "?", ""):
            continue
        chain = asym[i]
        if chain not in residues:
            residues[chain] = {}
            order.append(chain)
        residues[chain].setdefault(int(seq_id[i]), THREE_TO_ONE.get(comp[i], "X"))

    return [
        (chain, "".join(residues[chain][k] for k in sorted(residues[chain])))
        for chain in order
    ]


def safe_id(name):
    """ColabFold names output files after the record id, so keep it filesystem-safe."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", name)


def wrap(sequence):
    return "\n".join(
        sequence[i : i + LINE_WIDTH] for i in range(0, len(sequence), LINE_WIDTH)
    )


def main():
    os.makedirs("./outputs", exist_ok=True)

    # The parameter is gone, but fail loudly rather than silently ignoring a value that
    # someone set by hand expecting the old behaviour.
    mode = os.environ.get("PARAM_FASTA_MODE", "designed")
    if mode != "designed":
        fail(
            f"fasta_mode='{mode}' is no longer available; this node writes the designed "
            "chain alone. Colon-joined complex records are rejected by workflow-026's "
            "validator, which accepts only ACDEFGHIKLMNPQRSTVWYX, so the mode could not "
            "do what it promised."
        )

    designs = sorted(glob.glob("./inputs/*_generated.cif"))
    if not designs:
        fail("No generated structures found in ./inputs/ (*_generated.cif)")

    # Metrics are optional here: they only supply the ordering node 04 already decided on.
    order = {}
    metrics_path = "./inputs/design_metrics.json"
    if os.path.exists(metrics_path):
        try:
            with open(metrics_path) as handle:
                rows = json.load(handle).get("designs", [])
            order = {row["design"]: index for index, row in enumerate(rows)}
        except Exception as exc:
            print(f"WARNING: could not read design metrics: {exc}", flush=True)
    designs.sort(key=lambda p: (order.get(os.path.basename(p), 10**6), os.path.basename(p)))

    records = []
    skipped = []
    unknown_total = 0
    for path in designs:
        base = os.path.basename(path)[: -len(".cif")]
        chains = chain_sequences(path)
        if not chains:
            skipped.append((base, "no polymer chains"))
            continue
        sequence = chains[0][1]
        label = f"chain={chains[0][0]} len={len(sequence)}"
        residues = sequence
        if not residues:
            skipped.append((base, "empty sequence"))
            continue
        if set(residues) == {"X"}:
            skipped.append((base, "sequence is entirely unknown residues"))
            continue
        unknown_total += residues.count("X")
        records.append((safe_id(base), label, sequence))

    if not records:
        fail(
            "No design yielded a usable sequence: "
            + "; ".join(f"{name} ({why})" for name, why in skipped)
        )

    out_path = f"./outputs/{OUTPUT_NAME}"
    with open(out_path, "w") as handle:
        for record_id, label, sequence in records:
            handle.write(f">{record_id} {label}\n{wrap(sequence)}\n")

    for name, why in skipped:
        print(f"WARNING: skipped {name}: {why}", flush=True)
    if unknown_total:
        print(
            f"WARNING: {unknown_total} residue(s) across the batch are unknown (X). "
            "ColabFold accepts X but folds it as a gap in the sequence, so a design with "
            "many of them will refold poorly for that reason alone.",
            flush=True,
        )

    lengths = [len(seq) for _, _, seq in records]
    print(
        f"Wrote {len(records)} record(s) to {out_path}; "
        f"residues per record {min(lengths)}-{max(lengths)}",
        flush=True,
    )
    print(
        "Hand this file to workflow-026 as its sequences.fasta to refold the designs and "
        "filter on self-consistency; nothing in this workflow scores them.",
        flush=True,
    )


if __name__ == "__main__":
    main()
