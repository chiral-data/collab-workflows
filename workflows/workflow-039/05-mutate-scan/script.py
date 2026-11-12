#!/usr/bin/env python3
"""Node 05 — Mutation Scan.

Propose point mutants with AbLang2, then re-screen every proposal so none trades one risk for
another.

How the model is used: `ablang2.pretrained(model_to_use='ablang2-paired')` called with
`mode='probability'` returns, per clone, an array of shape [len(sequence), vocab] of softmax
probabilities over the amino acid at each position. The input string is `"VH|VL"` with the
heavy chain first, and `res_to_list` trims the output to the input length, so **row i is
character i of that string** and the separator occupies row len(VH). Vocabulary indices come
from `ablang.tokenizer.aa_to_token`, not from a hardcoded order.

With `scan_stepwise_masking=true` each scanned position is masked before scoring, which is the
more standard way to ask a masked language model what belongs at a position, at the cost of one
forward pass per position instead of one per clone.

A proposal is kept when:
  * the alternative residue scores higher than the wild-type residue at that position, and
  * the mutated sequence introduces no new liability motif, and
  * the substitution does not worsen the sequence-level charge and hydrophobicity picture in
    the direction the clone is already flagged for.

The second and third rules are the point of the node. A language model will happily suggest an
asparagine that completes an N-glycosylation sequon, or a lysine that deepens a positive patch
on a clone already flagged for PPC. Dropped proposals are written out with the reason, so the
screen is auditable rather than invisible.

**These are likelihood-based proposals, not affinity predictions.** No antigen is considered
anywhere in this workflow, so nothing here can speak to binding. A high-probability mutant is
one that looks like a natural antibody. Treat every row as a bench hypothesis.

The motif definitions come from node 01's annotations.json rather than being redeclared here,
so a custom motif set configured in node 01 governs the re-screen too.

Supports --dry-run: resolve and print every setting, then exit 0 without loading weights.
"""
import csv
import json
import os
import re
import sys
import time

os.environ.setdefault("HF_HUB_OFFLINE", "1")

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
SCAN_REGIONS = {"cdr", "cdrh3", "all"}

# Kyte & Doolittle comes from biopython at runtime; this is only the charge table, which the
# TAP paper specifies directly.
RESIDUE_CHARGE = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1}


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def warn(message):
    print(f"WARNING: {message}", flush=True)


def read_json(path, label, required=True):
    if not os.path.exists(path):
        if required:
            fail(f"{label} not found: {path}. The upstream node must run first.")
        return {}
    try:
        with open(path) as handle:
            return json.load(handle)
    except Exception as exc:
        fail(f"Could not read {label} ('{path}'): {exc}")


def package_version(name):
    try:
        from importlib.metadata import version

        return version(name)
    except Exception:
        return "unknown"


def write_empty(settings, reason):
    """Every declared output is always written: an empty table plus an explicit disabled flag."""
    with open("./outputs/proposals.csv", "w", newline="") as handle:
        csv.DictWriter(handle, fieldnames=PROPOSAL_FIELDS).writeheader()
    with open("./outputs/scan_manifest.json", "w") as handle:
        json.dump(
            {"enabled": False, "reason": reason, "settings": settings, "n_proposals": 0},
            handle,
            indent=2,
        )
    print(f"Scan not run: {reason}", flush=True)


PROPOSAL_FIELDS = [
    "clone_id",
    "chain",
    "position",
    "region",
    "wild_type",
    "mutant",
    "wt_probability",
    "mutant_probability",
    "probability_gain",
    "rank",
    "verdict",
    "reason",
]


def find_motif_hits(sequence, motifs):
    """Set of (motif name, start offset) hits, for before/after comparison."""
    hits = set()
    for motif in motifs:
        for match in re.finditer(f"(?={motif['pattern']})", sequence):
            hits.add((motif["name"], match.start()))
    return hits


def introduced_motifs(original, mutated, motifs):
    """Motif names present after the substitution that were not present before."""
    before = {name for name, _ in find_motif_hits(original, motifs)}
    after = find_motif_hits(mutated, motifs)
    return sorted({name for name, _ in after if name not in before})


def cysteine_parity_broken(original, mutated):
    return (original.count("C") % 2) != (mutated.count("C") % 2)


def developability_direction(wild_type, mutant, flags, hydrophobicity):
    """Reject a substitution that pushes a clone further in a direction it is flagged for.

    This is a sequence-level guard, not a recomputation of the structural metrics: node 03's
    metrics need a structure, and the mutant has none. It uses node 03's flags to decide which
    direction is dangerous for this particular clone.
    """
    wt_charge = RESIDUE_CHARGE.get(wild_type, 0.0)
    mutant_charge = RESIDUE_CHARGE.get(mutant, 0.0)
    charge_delta = mutant_charge - wt_charge
    hydrophobicity_delta = hydrophobicity.get(mutant, 0.0) - hydrophobicity.get(wild_type, 0.0)

    if flags.get("ppc_cdr_vicinity") in {"amber", "red"} and charge_delta > 0:
        return (
            f"clone is flagged {flags['ppc_cdr_vicinity']} for positive-charge patches and this "
            "substitution adds positive charge in a CDR"
        )
    if flags.get("pnc_cdr_vicinity") in {"amber", "red"} and charge_delta < 0:
        return (
            f"clone is flagged {flags['pnc_cdr_vicinity']} for negative-charge patches and this "
            "substitution adds negative charge in a CDR"
        )
    if flags.get("psh_cdr_vicinity") == "red" and hydrophobicity_delta > 0:
        return (
            "clone is flagged red for surface hydrophobicity and this substitution is more "
            "hydrophobic"
        )
    return None


def scannable_positions(chain_label, residues, scan_region):
    """Indices into the chain's domain sequence that the scan may vary."""
    selected = []
    for index, row in enumerate(residues):
        region = row.get("region", "")
        if scan_region == "all":
            keep = True
        elif scan_region == "cdr":
            keep = region.startswith("CDR")
        else:  # cdrh3
            keep = region == "CDR3" and chain_label == "VH"
        if keep:
            selected.append(index)
    return selected


def main():
    dry_run = "--dry-run" in sys.argv
    os.makedirs("./outputs", exist_ok=True)

    enabled = os.environ.get("PARAM_MUTATE_SCAN", "true").lower() == "true"
    scan_region = os.environ.get("PARAM_SCAN_REGION", "cdr")
    if scan_region not in SCAN_REGIONS:
        fail(f"scan_region must be one of {sorted(SCAN_REGIONS)}, got '{scan_region}'")
    stepwise = os.environ.get("PARAM_SCAN_STEPWISE_MASKING", "false").lower() == "true"
    max_mutants = int(os.environ.get("PARAM_MAX_MUTANTS_PER_CLONE", "20"))
    if max_mutants < 1:
        fail(f"max_mutants_per_clone must be at least 1, got {max_mutants}")
    seed = int(os.environ.get("PARAM_SEED", "42"))

    annotations = read_json("./inputs/annotations.json", "annotations.json")
    clones = annotations.get("clones") or []
    motifs = annotations.get("motifs") or []
    if not clones:
        fail("annotations.json contains no annotated clones.")
    if not motifs:
        warn(
            "annotations.json carries no motif definitions, so proposals cannot be re-screened "
            "for liabilities. Only the cysteine-parity and developability-direction rules will "
            "apply."
        )

    developability = read_json("./inputs/developability.json", "developability.json", required=False)
    flags_by_clone = {
        clone_id: (entry.get("flags") or {})
        for clone_id, entry in (developability.get("clones") or {}).items()
    }

    settings = {
        "mutate_scan": enabled,
        "scan_region": scan_region,
        "scan_stepwise_masking": stepwise,
        "max_mutants_per_clone": max_mutants,
        "seed": seed,
        "model": "ablang2-paired",
        "ablang2_version": package_version("ablang2"),
        "n_clones": len(clones),
        "n_motifs": len(motifs),
    }
    print(f"Resolved settings: {json.dumps(settings, indent=2)}", flush=True)

    if not enabled:
        write_empty(settings, "mutate_scan=false")
        return
    if dry_run:
        print("Dry run requested — not loading model weights.", flush=True)
        return

    from Bio.SeqUtils.ProtParamData import kd

    import ablang2

    try:
        model = ablang2.pretrained(
            model_to_use="ablang2-paired", random_init=False, ncpu=1, device="cpu"
        )
    except Exception as exc:
        fail(
            f"Could not load AbLang2: {exc}. The weights are baked into the image at build "
            "time, so this is an image-build problem rather than a network one."
        )
    vocabulary = dict(model.tokenizer.aa_to_token)

    rows, failures = [], []
    counts = {"kept": 0, "dropped": 0}
    started = time.time()

    for record in clones:
        clone_id = record["clone_id"]
        chains = record.get("chains") or {}
        heavy = (chains.get("VH") or {}).get("domain_sequence") or ""
        light = (chains.get("VL") or {}).get("domain_sequence") or ""
        if not heavy:
            failures.append({"clone_id": clone_id, "reason": "no heavy chain"})
            continue

        # AbLang2's paired model always takes "VH|VL"; a nanobody has an empty light chain.
        query = f"{heavy}|{light}"
        offsets = {"VH": 0, "VL": len(heavy) + 1}

        positions = []
        for chain_label in ("VH", "VL"):
            chain = chains.get(chain_label)
            if not chain:
                continue
            for index in scannable_positions(
                chain_label, chain.get("residues") or [], scan_region
            ):
                positions.append((chain_label, index))
        if not positions:
            failures.append(
                {
                    "clone_id": clone_id,
                    "reason": f"no positions in scan_region='{scan_region}'",
                }
            )
            continue

        try:
            if stepwise:
                probabilities = model(
                    [query], mode="probability", stepwise_masking=True
                )[0]
            else:
                probabilities = model([query], mode="probability")[0]
        except Exception as exc:
            failures.append({"clone_id": clone_id, "reason": f"{type(exc).__name__}: {exc}"})
            print(f"  {clone_id}: FAILED — {exc}", flush=True)
            continue

        if len(probabilities) < len(query):
            failures.append(
                {
                    "clone_id": clone_id,
                    "reason": (
                        f"AbLang2 returned {len(probabilities)} rows for a {len(query)}-"
                        "character query; the row-to-residue mapping cannot be trusted"
                    ),
                }
            )
            continue

        clone_flags = flags_by_clone.get(clone_id, {})
        candidates = []

        for chain_label, index in positions:
            chain = chains[chain_label]
            residues = chain.get("residues") or []
            sequence = chain.get("domain_sequence") or ""
            row_index = offsets[chain_label] + index
            if row_index >= len(probabilities) or index >= len(sequence):
                continue

            row = probabilities[row_index]
            wild_type = sequence[index]
            wt_token = vocabulary.get(wild_type)
            if wt_token is None or wt_token >= len(row):
                continue
            wt_probability = float(row[wt_token])

            for mutant in AMINO_ACIDS:
                if mutant == wild_type:
                    continue
                token = vocabulary.get(mutant)
                if token is None or token >= len(row):
                    continue
                mutant_probability = float(row[token])
                if mutant_probability <= wt_probability:
                    continue

                mutated = sequence[:index] + mutant + sequence[index + 1 :]
                verdict, reason = "kept", ""

                new_motifs = introduced_motifs(sequence, mutated, motifs)
                if new_motifs:
                    verdict = "dropped"
                    reason = f"introduces liability motif(s): {', '.join(new_motifs)}"
                elif cysteine_parity_broken(sequence, mutated):
                    verdict = "dropped"
                    reason = "changes cysteine parity, risking an unpaired thiol"
                else:
                    direction = developability_direction(
                        wild_type, mutant, clone_flags, kd
                    )
                    if direction:
                        verdict, reason = "dropped", direction

                residue_row = residues[index] if index < len(residues) else {}
                candidates.append(
                    {
                        "clone_id": clone_id,
                        "chain": chain_label,
                        "position": (
                            f"{residue_row.get('position', index + 1)}"
                            f"{residue_row.get('insertion', '')}"
                        ),
                        "region": residue_row.get("region", ""),
                        "wild_type": wild_type,
                        "mutant": mutant,
                        "wt_probability": round(wt_probability, 4),
                        "mutant_probability": round(mutant_probability, 4),
                        "probability_gain": round(mutant_probability - wt_probability, 4),
                        "verdict": verdict,
                        "reason": reason,
                    }
                )

        # Rank by probability gain; ties broken deterministically by position and residue so a
        # rerun with the same seed gives the same table.
        candidates.sort(
            key=lambda item: (
                -item["probability_gain"],
                item["chain"],
                str(item["position"]),
                item["mutant"],
            )
        )

        kept = [item for item in candidates if item["verdict"] == "kept"][:max_mutants]
        for rank, item in enumerate(kept, start=1):
            item["rank"] = rank
        dropped = [item for item in candidates if item["verdict"] == "dropped"]
        for item in dropped:
            item["rank"] = ""

        rows.extend(kept)
        rows.extend(dropped)
        counts["kept"] += len(kept)
        counts["dropped"] += len(dropped)

        print(
            f"  {clone_id}: {len(kept)} proposal(s) kept, {len(dropped)} dropped by the "
            f"re-screen, over {len(positions)} scanned position(s)",
            flush=True,
        )
        if candidates and not kept:
            print(
                f"    every proposal for {clone_id} was dropped — the re-screen is working; "
                "read the reasons in proposals.csv",
                flush=True,
            )

    with open("./outputs/proposals.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=PROPOSAL_FIELDS, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    manifest = {
        "enabled": True,
        "settings": settings,
        "n_proposals": counts["kept"],
        "n_dropped": counts["dropped"],
        "n_clones_scanned": len(clones) - len(failures),
        "failures": failures,
        "elapsed_seconds": round(time.time() - started, 1),
        "screen_rules": [
            "mutant probability must exceed the wild-type probability at that position",
            "no new liability motif may be introduced (motif set from node 01)",
            "cysteine parity must not change",
            "the substitution must not push a flagged clone further in its flagged direction",
        ],
        "caveat": (
            "Proposals are AbLang2 probabilities: they say a residue looks natural at that "
            "position, not that it binds better. No antigen is considered anywhere in this "
            "workflow, so these are bench hypotheses rather than affinity predictions."
        ),
    }
    with open("./outputs/scan_manifest.json", "w") as handle:
        json.dump(manifest, handle, indent=2)

    print(
        f"Proposed {counts['kept']} mutant(s) across "
        f"{manifest['n_clones_scanned']} clone(s); {counts['dropped']} dropped; "
        f"{manifest['elapsed_seconds']}s",
        flush=True,
    )
    if failures:
        warn(f"{len(failures)} clone(s) could not be scanned; see scan_manifest.json.")
    if not flags_by_clone:
        warn(
            "No developability flags were available, so the direction rule could not consider "
            "which risk each clone already carries. Run node 03 first for the full re-screen."
        )


if __name__ == "__main__":
    main()
