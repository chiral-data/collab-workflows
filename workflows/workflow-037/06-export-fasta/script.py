#!/usr/bin/env python3
"""Node 06 — Export FASTA.

Write the designed sequences in the filename and format workflow-026 (ColabFold) reads.

Two artifacts, because they answer different questions:

  sequences.fasta        every design, one record each. The batch record, and what you
                         read or archive. NOT what you hand to workflow-026.
  refold_<design>.fasta  one design per file. This is the unit workflow-026 can process
                         correctly, and there is one run of workflow-026 per file.

Why not hand workflow-026 the batch directly, when its input is named sequences.fasta and
would otherwise be a file copy: it processes one job per run, and a multi-record FASTA
silently becomes one conflated job rather than N independent ones. The prediction step
itself is fine — colabfold_batch makes each record its own job — but everything after it
assumes a single job:

  * 026/01_validate_fasta and 026/02_predict both set is_multimer = len(records) > 1, so N
    independent monomers are labelled "multimer" and pair_mode auto-resolves to
    unpaired_paired. Harmless to a single-chain prediction, wrong in every log line.
  * 026/02_predict copies every colabfold output flat into ./outputs/, so N designs' files
    land in one directory with nothing separating them.
  * 026/03_analyze_confidence then globs *_scores_rank_*.json across all of them, takes
    job_id from whichever file sorts first, sorts by rank alone and truncates to top_n. At
    3 designs x 5 models the "top 5" are ranks 1,1,1,2,2 drawn from three different
    designs, all reported under one design's name. Its chain_lengths comes from the whole
    FASTA, so three designs' lengths are read as one job's three chains.
  * 026/04_pymol_render takes its structure from an unsorted glob and its pLDDT from a
    rank-sorted list — two different orderings over a mixed set.

None of that fails loudly. It produces a confidence report about a design that does not
exist. One design per file is what makes the handoff mean what it says.

What refolding does and does not tell you: workflow-026 returns pLDDT and PAE — a
confidence estimate. It does not compare the predicted structure against the backbone this
workflow generated, and it cannot, because only the sequence travels. That comparison is
self-consistency RMSD, it is the standard test for a generative design, and nothing in
either workflow performs it today. A sequence can refold confidently into a fold that is
not the one it was designed as, and this handoff would call that a success.

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
PER_DESIGN_PREFIX = "refold_"
MANIFEST_NAME = "generation_manifest.json"

# Carried on every record's description line. sequences.fasta is the only artifact built
# to leave this workflow, so it is the one that has to say the designs were not screened:
# the manifest and the report both record it, and neither travels with the file.
#
# This is a description and not a leading ";" comment banner, deliberately. Biopython
# accepts comment lines today with a DeprecationWarning that says it will become a
# ValueError, and workflow-026 calls SeqIO.parse without a try/except — so a banner would
# break workflow-026 the next time its image is rebuilt, and look like a bug in 026.
#
# It must also contain no colon. workflow-026 reads a colon in a record's id or
# description as a chain-break marker and flips the run to multimer pairing.
UNSCREENED_TAG = "UNSCREENED - computational design, requires biosecurity review"
SCREENED_TAG = "screening recorded upstream - see generation_manifest.json"

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


def screening_tag():
    """The biosafety tag every record carries.

    Taken from node 03's manifest when there is one, so that if screening is ever added
    the tag follows it instead of contradicting it. Absent or unreadable manifest means
    unscreened: that is both the truth for the RFD4 CLI path and the safe default.
    """
    screened = None
    path = f"./inputs/{MANIFEST_NAME}"
    if os.path.exists(path):
        try:
            with open(path) as handle:
                screened = json.load(handle).get("securedna_screened")
        except Exception as exc:
            print(f"WARNING: could not read {MANIFEST_NAME}: {exc}", flush=True)
    tag = SCREENED_TAG if screened is True else UNSCREENED_TAG
    if ":" in tag:
        fail(
            "The biosafety tag contains a colon, which workflow-026 reads as a chain "
            "break. Rewrite it without one."
        )
    return tag


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

    # Node 04 has already decided the two things this node would otherwise decide
    # again: what order the designs go in, and which chain of each file is the design.
    # Both are read from its metrics rather than recomputed. The order is cosmetic; the
    # chain is not. Exporting the wrong one writes the target's sequence into the FASTA
    # under the design's name, which is wrong, silent and entirely plausible-looking.
    order = {}
    designed_chain = {}
    metrics_path = "./inputs/design_metrics.json"
    if os.path.exists(metrics_path):
        try:
            with open(metrics_path) as handle:
                rows = json.load(handle).get("designs", [])
            order = {row["design"]: index for index, row in enumerate(rows)}
            designed_chain = {
                row["design"]: row.get("designed_chain")
                for row in rows
                if row.get("designed_chain")
            }
        except Exception as exc:
            print(f"WARNING: could not read design metrics: {exc}", flush=True)
    def natural_key(name):
        """Sort key that reads runs of digits as numbers, so design_2 precedes design_10."""
        return [int(p) if p.isdigit() else p for p in re.split(r"(\d+)", name)]

    designs.sort(key=lambda p: (order.get(os.path.basename(p), 10**6),
                                natural_key(os.path.basename(p))))

    records = []
    skipped = []
    unknown_total = 0
    for path in designs:
        base = os.path.basename(path)[: -len(".cif")]
        chains = chain_sequences(path)
        if not chains:
            skipped.append((base, "no polymer chains"))
            continue
        wanted = designed_chain.get(os.path.basename(path))
        picked = next((pair for pair in chains if pair[0] == wanted), None)
        if wanted and picked is None:
            print(
                f"WARNING: node 04 named chain {wanted} as the design in {base}, but "
                f"this file has {[cid for cid, _ in chains]}. Exporting the first chain.",
                flush=True,
            )
        if picked is None:
            if not wanted and len(chains) > 1:
                print(
                    f"WARNING: no designed chain recorded for {base}, so the first chain "
                    "in the file is being exported.",
                    flush=True,
                )
            picked = chains[0]
        chain_id, sequence = picked
        label = f"chain={chain_id} len={len(sequence)}"
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

    tag = screening_tag()
    out_path = f"./outputs/{OUTPUT_NAME}"
    with open(out_path, "w") as handle:
        for record_id, label, sequence in records:
            handle.write(f">{record_id} {label} {tag}\n{wrap(sequence)}\n")

    # One file per design: the unit workflow-026 can process without conflating designs.
    # See the module docstring for what goes wrong when it is handed the batch instead.
    for record_id, label, sequence in records:
        with open(f"./outputs/{PER_DESIGN_PREFIX}{record_id}.fasta", "w") as handle:
            handle.write(f">{record_id} {label} {tag}\n{wrap(sequence)}\n")

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
    print(f"Every record is tagged: {tag}", flush=True)
    print(
        f"Also wrote {len(records)} single-design file(s) {PER_DESIGN_PREFIX}*.fasta.",
        flush=True,
    )
    print(
        "To refold: run workflow-026 once per " + PER_DESIGN_PREFIX + "*.fasta file, "
        "renaming it to sequences.fasta. Do NOT hand it the batch file above — "
        "workflow-026 reports one job per run, so several records in one file are "
        "conflated into a single mislabelled result rather than refused.",
        flush=True,
    )
    print(
        "What comes back is pLDDT and PAE, a confidence estimate. It is not "
        "self-consistency: nothing compares the prediction against the backbone generated "
        "here, because only the sequence travels. A sequence can refold confidently into "
        "the wrong fold.",
        flush=True,
    )


if __name__ == "__main__":
    main()
