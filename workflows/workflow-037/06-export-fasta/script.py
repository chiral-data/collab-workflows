#!/usr/bin/env python3
"""Node 06 — Export FASTA.

Write the designed sequences as `sequences.fasta`, which is exactly the filename and format
workflow-026 (ColabFold) takes as its input, so refolding a batch of designs is a file copy
rather than a reformatting job.

Why this node exists: this workflow runs with `~metrics` and scores nothing, so a design's
only real test is independent refolding — fold the sequence, check the fold matches the
backbone RFD4 drew (self-consistency), and keep the designs that survive. That filter lives
in workflow-026; this node is the handoff.

`fasta_mode` chooses what each record holds:

  designed  one record per design, the generated chain alone. Refolding this answers "does
            this sequence fold to the shape it was designed as?", which is the standard
            self-consistency check and the right default.
  complex   one record per design, every chain colon-joined (designed first, then target).
            ColabFold reads the colon as a chain break and folds the complex, which also
            gives interface confidence — slower, and only meaningful when there is a target.

Note on how workflow-026 reads the result: it labels any multi-record file a "multimer" and
therefore resolves `pair_mode=auto` to `unpaired_paired`. That is harmless for single-chain
records (there is nothing to pair) and correct for colon-joined ones, so both modes are safe
to hand over as-is.
"""
import glob
import json
import os
import re
import sys

MODES = ("designed", "complex")
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

    mode = os.environ.get("PARAM_FASTA_MODE", "designed")
    if mode not in MODES:
        fail(f"Unknown fasta_mode '{mode}'. Expected one of: {', '.join(MODES)}")

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
        if mode == "designed":
            sequence = chains[0][1]
            label = f"chain={chains[0][0]} len={len(sequence)}"
        else:
            sequence = ":".join(seq for _, seq in chains)
            label = (
                f"chains={','.join(c for c, _ in chains)} "
                f"len={','.join(str(len(s)) for _, s in chains)}"
            )
        residues = sequence.replace(":", "")
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

    lengths = [len(seq.replace(":", "")) for _, _, seq in records]
    print(f"Mode: {mode}", flush=True)
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
