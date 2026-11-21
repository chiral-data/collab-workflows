#!/usr/bin/env python3
"""Read the designed structures and say, per design, which chain is the binder.

This node exists so that nothing downstream ever has to decide which chain is the
design. That decision was already made upstream, twice: workflow-037's node 03 recorded
the target's chain ids and its node 04 resolved the designed chain from them. Both
answers are in design_metrics.json. Read them.

Contract
--------
Reads   *_generated.cif          workflow-037 node 03, one structure per design
        design_metrics.json      workflow-037 node 04, optional
Writes  design_inventory.json    one entry per design
        design_metrics.json      passed through unchanged for nodes 04 and 05 to merge

design_inventory.json, per design:

    design                  the design's name, which is what everything downstream keys on
    source_file             the CIF this entry was read from
    binder_chain            chain id of the designed chain
    binder_sequence         one-letter, in chain order
    binder_residue_ids      label_seq_id values, IN CHAIN ORDER
    binder_ca               Ca coordinates, in that same order, null where absent
    target_chains           [{chain, sequence, residue_ids, gaps, n_residues}, ...]
    designed_chain_source   carried through VERBATIM from design_metrics.json
    n_unknown               residues with no one-letter code, scored X
    skipped / skip_reason   true when there is nothing foldable here

The Ca coordinates live here on purpose. Node 04 needs the designed backbone to compute
an RMSD against, and the alternative -- passing the CIFs through and letting node 04
parse them -- would make node 04 resolve the binder chain a second time. That is the
defect class this repo keeps fixing. One node reads the structures; nobody else parses
them.

binder_residue_ids, binder_sequence and binder_ca are index-aligned, always, including
where a residue has no CA atom: that position holds null rather than being dropped.
Dropping it would silently shift every later index, and node 04's RMSD maps the two
structures BY POSITION (§6, ef74646e). The count of such positions is reported as
n_binder_ca_missing, so node 04 can state how many residues its RMSD is over.

designed_chain_source
---------------------
Three values, and they are workflow-037's, not ours (04-analyze/script.py:248). Do not
invent a fourth and do not collapse them:

    manifest-target-chains   resolved by set difference against the recorded target
                             chains. The good path.
    single-chain             one chain in the file, so it is the design by construction:
                             a monomer or motif_scaffold run. Valid for self-consistency.
                             Node 02 must NOT write a complex input for it -- there is no
                             target to score an interface against.
    assumed-first-chain      037 declined to answer and fell back to chain order.
                             Proceed, say so loudly, and carry the value forward.

Absent design_metrics.json is the stated-fallback path, not an error: record
assumed-first-chain, warn, and continue.

Why the fallback matters more than it looks
-------------------------------------------
037's first-chain fallback was in force on EVERY real run for weeks. Its node 03 recorded
target_chains by listing every label_asym_id in target.cif, which on 4ZQK is the target
chain, the chain the contig dropped, a sodium and two solvent groups -- so the designed
chain was in that list, the set difference found no candidate, and the fallback fired
every time (ef74646e). And it was right every time, because the design does come out
first. Nothing failed. Nothing looked wrong.

So this node's own test is not that the fallback works. It is that the good path was
taken: tests/fixtures.json, ingest/resolved-chain-good-path.

Why this node validates parameters it never uses
------------------------------------------------
It is the only node guaranteed to run before the GPU node, and every threshold in this
workflow is read after the predictor has already spent its minutes. A nan in
max_sc_rmsd_a would otherwise be discovered by node 05, GPU-hours later, and §8.8 is the
record of what a nan threshold does when it is not refused: it passes "must be positive"
and passes every range check, because every comparison against nan is false. Nodes 04 and
05 parse their own thresholds with the same helper, so a node run standalone is still
safe; this node is the gate, not the only check.

Fails when
----------
  * no *_generated.cif inputs
  * design_metrics.json names a chain the structure does not contain
  * a design's binder chain resolves to zero residues
  * two designs resolve to the same name
  * any parameter is non-finite, out of range, or the wrong type
  * score_interface is true, which is Phase 2 and not built

Does NOT fail when design_metrics.json is absent, or when a target chain carries internal
gaps -- gaps are measured and reported, because node 03 has to disclose that it folds the
chain whole and fuses them. Trimming around a gap was measured to be much worse
(8e4772a7, §8.12).
"""
import glob
import json
import math
import os
import re
import shutil
import sys

# 037's node 04 table, plus the three non-standard codes it carries. A residue outside it
# becomes X, and X is counted rather than quietly folded: ColabFold treats X as a sequence
# gap, which depresses pLDDT for reasons that have nothing to do with the design
# (54bf21b6).
THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}

# ColabFold names every output file after the record id, so the design name has to survive
# a round trip through the filesystem unchanged (54bf21b6). Anything outside this set is
# replaced with '_', and a collision after replacement is an error rather than an
# overwrite.
SAFE_NAME_RE = re.compile(r"[^A-Za-z0-9._-]")

GENERATED_SUFFIX = "_generated.cif"


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def numeric_param(env_name, display_name, default, cast):
    """Read a numeric parameter, failing the same way every other bad parameter does.

    Casting inline at the point of use turns a typo into a raw ValueError traceback, which
    is the one error in a node that does not look like the others (975ab379).
    """
    raw = os.environ.get(env_name)
    if raw is None or raw.strip() == "":
        return default
    try:
        value = cast(raw)
    except ValueError:
        kind = "a whole number" if cast is int else "a number"
        fail(f"{display_name} must be {kind}, got '{raw}'.")
    # float() accepts 'nan', 'inf' and anything that overflows to inf, and none of them
    # are caught by a range check: every comparison against nan is false, so a nan cutoff
    # passes "must be positive" and "outside the expected range" alike and then makes
    # every distance test false (55edd8fc). inf does the opposite and takes everything.
    # The parsed value is named in the message because '1e400' does not look infinite and
    # a user who typed it has no other way to find out that it became inf.
    if not math.isfinite(value):
        fail(
            f"{display_name} must be a finite number, got '{raw}', "
            f"which parses as {value}."
        )
    return value


def positive_param(env_name, display_name, default, cast):
    value = numeric_param(env_name, display_name, default, cast)
    if value <= 0:
        fail(f"{display_name} must be positive, got {value}.")
    return value


def bool_param(env_name, display_name, default):
    raw = os.environ.get(env_name)
    if raw is None or raw.strip() == "":
        return default
    lowered = raw.strip().lower()
    if lowered in ("true", "1", "yes", "on"):
        return True
    if lowered in ("false", "0", "no", "off"):
        return False
    fail(f"{display_name} must be true or false, got '{raw}'.")


def check_parameters():
    """Refuse every bad parameter here, before the GPU node is reached.

    Returns nothing. The values are re-read by the nodes that use them; what this
    function is for is the failure, and the failure has to happen first.
    """
    positive_param("PARAM_MAX_SC_RMSD_A", "max_sc_rmsd_a", 2.0, float)
    plddt = numeric_param("PARAM_MIN_PLDDT_BINDER", "min_plddt_binder", 80.0, float)
    if not 0.0 <= plddt <= 100.0:
        fail(f"min_plddt_binder is a pLDDT, so it must be 0-100, got {plddt}.")
    positive_param("PARAM_NUM_MODELS", "num_models", 5, int)
    positive_param("PARAM_NUM_RECYCLE", "num_recycle", 3, int)
    positive_param("PARAM_PREDICT_TIMEOUT_S", "predict_timeout_s", 1800, int)
    positive_param("PARAM_TOP_N", "top_n", 10, int)

    msa_mode = (os.environ.get("PARAM_MSA_MODE") or "single_sequence").strip()
    if msa_mode != "single_sequence":
        fail(
            f"msa_mode '{msa_mode}' is not available. Phase 1 predicts the binder alone, "
            "where there is no alignment to find -- a de novo sequence has no homologs by "
            "construction (§4). The target-only alignment is Phase 2 and arrives with "
            "01b-align-target."
        )

    if bool_param("PARAM_SCORE_INTERFACE", "score_interface", False):
        fail(
            "score_interface = true is Phase 2 and is not built. Refusing rather than "
            "ignoring it: a run that silently dropped the flag would report "
            "self-consistency numbers under the heading the user asked for, and nothing "
            "in the output would say the interface was never scored (§8.7, §12)."
        )


def read_chains(path):
    """Return [(chain_id, [residue, ...]), ...] in file order, residues in label order.

    Each residue is {"seq_id": int, "one_letter": str, "ca": (x, y, z) or None}.
    """
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

    order, residues = [], {}
    for i in range(n_atoms):
        if seq_id[i] in (".", "?", ""):
            continue
        try:
            key = int(seq_id[i])
        except ValueError:
            continue
        chain = asym[i]
        if chain not in residues:
            residues[chain] = {}
            order.append(chain)
        entry = residues[chain].setdefault(
            key, {"seq_id": key, "one_letter": THREE_TO_ONE.get(comp[i], "X"), "ca": None}
        )
        if atom_id[i] == "CA":
            try:
                entry["ca"] = (float(xs[i]), float(ys[i]), float(zs[i]))
            except ValueError:
                entry["ca"] = None

    chains = []
    for chain in order:
        chains.append((chain, [residues[chain][k] for k in sorted(residues[chain])]))
    if not chains:
        fail(f"'{path}' contains no polymer residues.")
    return chains


def gaps_in(residues):
    """Internal breaks in label numbering, as [{after, before, missing}, ...].

    Numbering gaps, not distance gaps: a disordered loop absent from the coordinates is
    what this is for, and it is ordinary in anything from RCSB. Node 03 folds the chain
    whole and discloses that each break is fused, because trimming around one was
    measured at ~42 pLDDT against 81 for the whole chain (8e4772a7, §8.12).
    """
    found = []
    for before, after in zip(residues, residues[1:]):
        step = after["seq_id"] - before["seq_id"]
        if step > 1:
            found.append({
                "after": before["seq_id"],
                "before": after["seq_id"],
                "missing": step - 1,
            })
    return found


def resolve_binder(chains, row, path):
    """Return (index, source). Read 037's answer; do not re-derive it.

    The one thing this does beyond reading is check that the named chain exists. Reading
    an upstream answer does not mean trusting it blindly: a metrics file naming a chain
    the structure does not contain is a contract violation, and it has to fail here
    rather than resolve to something plausible two nodes later.
    """
    ids = [cid for cid, _ in chains]

    if row is not None and row.get("designed_chain"):
        named = row["designed_chain"]
        if named not in ids:
            fail(
                f"design_metrics.json names designed chain '{named}' for "
                f"'{os.path.basename(path)}', but that chain is not in the structure, "
                f"which has {ids}. Refusing to guess: the upstream answer and the "
                "structure disagree, and picking either one silently would make the "
                "mismatch invisible from here on."
            )
        source = row.get("designed_chain_source")
        if not source:
            # A row without the field is 037 before 252ee072. Taking its designed_chain is
            # still right; claiming it was resolved from the manifest is not.
            print(
                f"WARNING: the metrics row for '{os.path.basename(path)}' names a "
                "designed chain but carries no designed_chain_source. Recording "
                "'assumed-first-chain' would be wrong and 'manifest-target-chains' is a "
                "claim this node cannot check, so the value is passed through empty and "
                "the report shows it empty.",
                flush=True,
            )
        return ids.index(named), source

    if len(chains) == 1:
        # A monomer or motif_scaffold run comes back as the design alone. 037 deliberately
        # does not report this as a fallback, because doing so made a correct manifest
        # look like a missing one (04-analyze/script.py:265).
        return 0, "single-chain"

    print(
        f"WARNING: no metrics row for '{os.path.basename(path)}', so its designed chain "
        f"is being taken as the first chain in the file, '{ids[0]}'. That is true of what "
        "workflow-037's contigs produce, and it is not checked here. This is the "
        "assumed-first-chain path: it was in force on every real 037 run for weeks and "
        "was right every time, which is exactly why it has to stay visible (§8.2, "
        "ef74646e). The value reaches design_scores.json and the report.",
        flush=True,
    )
    return 0, "assumed-first-chain"


def ingest(path, row):
    chains = read_chains(path)
    index, source = resolve_binder(chains, row, path)
    binder_id, binder = chains[index]
    targets = [pair for i, pair in enumerate(chains) if i != index]

    if not binder:
        fail(
            f"The designed chain '{binder_id}' of '{os.path.basename(path)}' has no "
            "polymer residues. There is nothing to refold."
        )

    sequence = "".join(r["one_letter"] for r in binder)
    n_unknown = sequence.count("X")
    skip_reason = None
    if not sequence:
        skip_reason = "the designed chain has no sequence"
    elif n_unknown == len(sequence):
        # 037 node 06 skips rather than writing junk (54bf21b6). The design still appears
        # in the inventory, marked: a design that failed must appear as failed, not vanish
        # (§8.5).
        skip_reason = (
            f"all {len(sequence)} residues of the designed chain are non-standard, so the "
            "sequence is entirely X and there is nothing for a predictor to fold"
        )

    stem = os.path.basename(path)
    if stem.endswith(GENERATED_SUFFIX):
        stem = stem[: -len(GENERATED_SUFFIX)]
    name = SAFE_NAME_RE.sub("_", stem)
    if name != stem:
        print(
            f"NOTE: design name '{stem}' was rewritten to '{name}' for the filesystem; "
            "ColabFold names every output file after the record id.",
            flush=True,
        )

    binder_gaps = gaps_in(binder)
    if binder_gaps:
        print(
            f"WARNING: the designed chain '{binder_id}' of '{os.path.basename(path)}' has "
            f"{len(binder_gaps)} internal numbering gap(s). A design arrives from the "
            "generator as one continuous chain, so this means something renumbered or "
            "truncated it between there and here. The RMSD is computed over the residues "
            "present, by position.",
            flush=True,
        )

    return {
        "design": name,
        "source_file": os.path.basename(path),
        "binder_chain": binder_id,
        "binder_sequence": sequence,
        "binder_residue_ids": [r["seq_id"] for r in binder],
        "binder_ca": [list(r["ca"]) if r["ca"] else None for r in binder],
        "n_binder_residues": len(binder),
        "n_binder_ca_missing": sum(1 for r in binder if r["ca"] is None),
        "binder_gaps": binder_gaps,
        "target_chains": [
            {
                "chain": cid,
                "sequence": "".join(r["one_letter"] for r in residues),
                "residue_ids": [r["seq_id"] for r in residues],
                "gaps": gaps_in(residues),
                "n_residues": len(residues),
            }
            for cid, residues in targets
        ],
        "n_chains": len(chains),
        "designed_chain_source": source,
        "n_unknown": n_unknown,
        "skipped": skip_reason is not None,
        "skip_reason": skip_reason,
    }


def natural_key(name):
    """Sort runs of digits as numbers, so design_10 does not land between 1 and 2."""
    return [int(p) if p.isdigit() else p for p in re.split(r"(\d+)", name)]


def load_metrics(path):
    """Return (payload, {design_file_name: row}). An absent file is not an error."""
    if not os.path.exists(path):
        print(
            "WARNING: no design_metrics.json beside the structures. Workflow-037 node 04 "
            "writes one, and it carries the answer to which chain is the design. Without "
            "it every design falls back to the first chain in the file. That is a stated "
            "fallback and not an error, but it is recorded per design and it reaches the "
            "report (§8.2).",
            flush=True,
        )
        return None, {}
    try:
        with open(path) as handle:
            payload = json.load(handle)
    except Exception as exc:
        # Deliberately not fatal, and deliberately not silent. A corrupt metrics file
        # leaves us on the fallback path, which is where an absent one leaves us too, and
        # the distinction is in the log and in designed_chain_source.
        print(f"WARNING: could not read '{path}': {exc}. Falling back per design.",
              flush=True)
        return None, {}

    index = {}
    for row in payload.get("designs") or []:
        key = row.get("design")
        if key:
            index[os.path.basename(key)] = row
    print(f"Read design_metrics.json: {len(index)} row(s)", flush=True)
    return payload, index


def main():
    check_parameters()
    os.makedirs("./outputs", exist_ok=True)

    designs = sorted(
        glob.glob(f"./inputs/*{GENERATED_SUFFIX}"),
        key=lambda p: natural_key(os.path.basename(p)),
    )
    if not designs:
        fail(
            f"No generated structures found in ./inputs/ (*{GENERATED_SUFFIX}). This node "
            "reads workflow-037 node 03's output -- the structures, not the exported "
            "FASTA: a sequence alone cannot be compared against the backbone it was "
            "designed as, which is the whole measurement this workflow makes."
        )

    metrics_path = "./inputs/design_metrics.json"
    metrics_payload, rows = load_metrics(metrics_path)

    print(f"Found {len(designs)} generated design(s)", flush=True)

    entries, seen = [], {}
    for path in designs:
        entry = ingest(path, rows.get(os.path.basename(path)))
        if entry["design"] in seen:
            fail(
                f"'{entry['source_file']}' and '{seen[entry['design']]}' both resolve to "
                f"the design name '{entry['design']}'. Every output file downstream is "
                "named after it, so one design's prediction would overwrite the other's "
                "and the report would describe a design that does not exist (818954e9). "
                "Rename one of the inputs."
            )
        seen[entry["design"]] = entry["source_file"]
        entries.append(entry)

    sources = sorted({str(e["designed_chain_source"]) for e in entries})
    fused = sum(len(t["gaps"]) for e in entries for t in e["target_chains"])
    inventory = {
        "n_designs": len(entries),
        "n_skipped": sum(1 for e in entries if e["skipped"]),
        # Named here, once, so nodes 04 and 05 report the provenance of the chain choice
        # without re-deriving it or re-wording it.
        "designed_chain_sources": sources,
        "metrics_present": metrics_payload is not None,
        "n_target_numbering_gaps": fused,
        # Carried forward rather than looked up again: 037 is the only thing that knows
        # these, and a run whose provenance is unrecorded cannot be reproduced.
        "seed": (metrics_payload or {}).get("seed"),
        "hf_ckpt_revision": (metrics_payload or {}).get("hf_ckpt_revision"),
        "rfd4_proteina_version": (metrics_payload or {}).get("rfd4_proteina_version"),
        # 037's designs carry no SecureDNA screening, and passing validation makes
        # synthesis MORE likely, not less. The notice travels into every artifact.
        "securedna_screened": False,
        "designs": entries,
    }

    with open("./outputs/design_inventory.json", "w") as handle:
        json.dump(inventory, handle, indent=2)

    # Passed through unchanged, because nodes 04 and 05 merge 037's geometric metrics with
    # the confidence numbers computed here, and the merge has to be by design name against
    # the real upstream file rather than a re-derived copy of it.
    if metrics_payload is not None:
        shutil.copyfile(metrics_path, "./outputs/design_metrics.json")

    for entry in entries:
        targets = ",".join(t["chain"] for t in entry["target_chains"]) or "none"
        print(
            f"  {entry['design']}: binder chain {entry['binder_chain']} "
            f"len={entry['n_binder_residues']} (source: {entry['designed_chain_source']}), "
            f"target chains {targets}"
            + (f", {entry['n_unknown']} unknown residue(s)" if entry["n_unknown"] else "")
            + (f"  [SKIPPED: {entry['skip_reason']}]" if entry["skipped"] else ""),
            flush=True,
        )

    if fused:
        print(
            f"WARNING: the target chains carry {fused} internal numbering gap(s) in total. "
            "ColabFold cannot express a gap inside a chain, so node 03 folds each target "
            "chain WHOLE and the junction across every break is fused -- a real defect, "
            "kept knowingly. Trimming around the gap instead was measured at ~42 pLDDT "
            "where the whole chain reaches 81, and moved an interface PAE floor from 7.6 "
            "to 15.4 A against a threshold of 10 (8e4772a7, §8.12).",
            flush=True,
        )
    if "assumed-first-chain" in sources:
        n = sum(1 for e in entries if e["designed_chain_source"] == "assumed-first-chain")
        print(
            f"WARNING: {n} of {len(entries)} design(s) had their designed chain ASSUMED "
            "to be the first chain in the file rather than resolved. Every number this "
            "workflow reports for those designs describes whichever chain came first. "
            "Check that workflow-037's design_metrics.json reached this run.",
            flush=True,
        )
    skipped = [e["design"] for e in entries if e["skipped"]]
    if skipped:
        print(
            f"WARNING: {len(skipped)} design(s) have nothing foldable and will not be "
            f"predicted: {', '.join(skipped)}. They stay in the inventory, marked, so the "
            "report accounts for every design that entered (§8.5).",
            flush=True,
        )

    print("Wrote ./outputs/design_inventory.json", flush=True)


if __name__ == "__main__":
    main()
