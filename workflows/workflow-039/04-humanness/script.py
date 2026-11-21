#!/usr/bin/env python3
"""Node 04 — Humanness.

How human each variable domain looks, by two independent measures:

  1. **Peptide content** — the fraction of the domain's 9-mers that appear in the reference
     database. With `human-oas` (the default) the reference is antibody peptides observed in at
     least 10% of human subjects, which makes this an OASis-like score.
  2. **Germline identity** — ANARCI's closest-germline V and J identities, carried from node 01.

The two disagree usefully: a germline-identical framework with an engineered CDR scores high on
identity and lower on peptide content, and that is exactly the clone a humanization campaign
wants to see.

**Humanness is not immunogenicity.** No MHC-II binding or T-cell epitope prediction happens
anywhere in this workflow, because no open, commercially usable class-II predictor exists:
NetMHCIIpan is academic-use-only and mhcflurry is class I. See "the MHC-II hole" in the
pipeline doc. Treating a humanness score as an immunogenicity prediction is the single most
likely way to misread this node's output.

This node depends on node 01 and not on node 02, deliberately: it still produces a score for
clones whose structure prediction failed.

promb's reference databases are baked into the image, so nothing here touches the network.
"""
import csv
import json
import os
import sys

DATABASES = {"human-oas", "human-reference", "human-swissprot"}

# Only `human-oas` has a fixed peptide length, set by the database itself.
FIXED_LENGTH_DATABASES = {"human-oas"}


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def warn(message):
    print(f"WARNING: {message}", flush=True)


def read_json(path, label):
    if not os.path.exists(path):
        fail(f"{label} not found: {path}. Node 01 must run first.")
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


def open_database(name, peptide_length):
    from promb import init_db

    if name in FIXED_LENGTH_DATABASES:
        print(f"Opening promb database '{name}' (peptide length fixed at 9)", flush=True)
        return init_db(name)
    print(f"Opening promb database '{name}' with peptide length {peptide_length}", flush=True)
    return init_db(name, peptide_length)


def main():
    os.makedirs("./outputs", exist_ok=True)

    database_name = os.environ.get("PARAM_HUMANNESS_DB", "human-oas")
    if database_name not in DATABASES:
        fail(f"humanness_db must be one of {sorted(DATABASES)}, got '{database_name}'")
    peptide_length = int(os.environ.get("PARAM_PEPTIDE_LENGTH", "9"))
    if peptide_length < 3:
        fail(f"peptide_length must be at least 3, got {peptide_length}")
    if database_name in FIXED_LENGTH_DATABASES and peptide_length != 9:
        print(
            f"NOTE: peptide_length={peptide_length} is ignored for '{database_name}', whose "
            "peptide length is fixed at 9 by the database itself.",
            flush=True,
        )
        peptide_length = 9

    annotations = read_json("./inputs/annotations.json", "annotations.json")
    clones = annotations.get("clones") or []
    if not clones:
        fail("annotations.json contains no annotated clones.")

    settings = {
        "humanness_db": database_name,
        "peptide_length": peptide_length,
        "promb_version": package_version("promb"),
        "n_clones": len(clones),
    }
    print(f"Resolved settings: {json.dumps(settings, indent=2)}", flush=True)

    try:
        database = open_database(database_name, peptide_length)
    except Exception as exc:
        fail(
            f"Could not open the promb database '{database_name}': {exc}. The predefined "
            "databases are baked into the image at build time, so this is an image-build "
            "problem rather than a network one."
        )

    rows, detail, peptide_rows, failures = [], {}, [], []

    for record in clones:
        clone_id = record["clone_id"]
        clone_detail = {"chains": {}, "receptor_type": record.get("receptor_type")}
        chain_scores = {}

        for chain_label, chain in (record.get("chains") or {}).items():
            sequence = chain.get("domain_sequence") or ""
            if not sequence:
                continue
            try:
                content = float(database.compute_peptide_content(sequence))
                peptides = list(database.chop_seq_peptides(sequence))
                human_flags = [bool(database.contains(peptide)) for peptide in peptides]
            except Exception as exc:
                failures.append(
                    {
                        "clone_id": clone_id,
                        "chain": chain_label,
                        "reason": f"{type(exc).__name__}: {exc}",
                    }
                )
                print(f"  {clone_id} {chain_label}: FAILED — {exc}", flush=True)
                continue

            non_human = [
                peptide for peptide, is_human in zip(peptides, human_flags) if not is_human
            ]

            # The nearest human peptide for each non-human one is what tells a humanization
            # campaign where to start. It is an exhaustive search, so it is the slow path here.
            nearest_by_peptide = {}
            if non_human:
                try:
                    nearest = database.find_nearest_peptides(non_human, 1)
                    for peptide, hits in zip(non_human, nearest):
                        nearest_by_peptide[peptide] = (hits or [""])[0]
                except Exception as exc:
                    warn(
                        f"{clone_id} {chain_label}: nearest-peptide search failed ({exc}). "
                        "Peptide content is unaffected; peptide_detail.csv will not carry "
                        "nearest hits for this chain."
                    )

            offset_lookup = {}
            for index, peptide in enumerate(peptides):
                offset_lookup.setdefault(peptide, index)

            for peptide, is_human in zip(peptides, human_flags):
                if is_human:
                    continue
                start = offset_lookup.get(peptide, 0)
                covered = (chain.get("residues") or [])[start : start + peptide_length]
                regions = sorted({row.get("region", "?") for row in covered}) if covered else []
                peptide_rows.append(
                    {
                        "clone_id": clone_id,
                        "chain": chain_label,
                        "peptide": peptide,
                        "start_offset": start,
                        "position": (
                            f"{covered[0]['position']}{covered[0].get('insertion', '')}"
                            if covered
                            else ""
                        ),
                        "regions": ",".join(regions),
                        "in_cdr": any(region.startswith("CDR") for region in regions),
                        "nearest_human_peptide": nearest_by_peptide.get(peptide, ""),
                    }
                )

            chain_scores[chain_label] = content
            clone_detail["chains"][chain_label] = {
                "peptide_content": round(content, 4),
                "n_peptides": len(peptides),
                "n_human_peptides": sum(human_flags),
                "n_non_human_peptides": len(non_human),
                "v_gene": chain.get("v_gene"),
                "v_identity": chain.get("v_identity"),
                "j_gene": chain.get("j_gene"),
                "j_identity": chain.get("j_identity"),
                "species": chain.get("species"),
            }

        if not chain_scores:
            rows.append(
                {
                    "clone_id": clone_id,
                    "status": "failed",
                    "reason": "no chain could be scored",
                }
            )
            continue

        heavy = clone_detail["chains"].get("VH", {})
        light = clone_detail["chains"].get("VL", {})
        # Mean over the chains present, so a nanobody is not penalised for having one.
        mean_content = sum(chain_scores.values()) / len(chain_scores)
        identities = [
            chain.get("v_identity")
            for chain in clone_detail["chains"].values()
            if chain.get("v_identity") is not None
        ]
        mean_v_identity = sum(identities) / len(identities) if identities else None
        species = sorted(
            {
                chain.get("species")
                for chain in clone_detail["chains"].values()
                if chain.get("species")
            }
        )

        clone_detail["mean_peptide_content"] = round(mean_content, 4)
        clone_detail["mean_v_identity"] = (
            round(mean_v_identity, 4) if mean_v_identity is not None else None
        )
        clone_detail["germline_species"] = species
        detail[clone_id] = clone_detail

        rows.append(
            {
                "clone_id": clone_id,
                "receptor_type": record.get("receptor_type"),
                "mean_peptide_content": round(mean_content, 4),
                "vh_peptide_content": (
                    round(heavy["peptide_content"], 4) if heavy else ""
                ),
                "vl_peptide_content": (
                    round(light["peptide_content"], 4) if light else ""
                ),
                "vh_v_gene": heavy.get("v_gene") or "",
                "vh_v_identity": heavy.get("v_identity") if heavy.get("v_identity") is not None else "",
                "vl_v_gene": light.get("v_gene") or "",
                "vl_v_identity": light.get("v_identity") if light.get("v_identity") is not None else "",
                "mean_v_identity": (
                    round(mean_v_identity, 4) if mean_v_identity is not None else ""
                ),
                "germline_species": ",".join(species),
                "n_non_human_peptides": sum(
                    chain.get("n_non_human_peptides", 0)
                    for chain in clone_detail["chains"].values()
                ),
                "status": "ok",
                "reason": "",
            }
        )
        print(
            f"  {clone_id}: peptide content {mean_content:.3f}, V identity "
            f"{f'{mean_v_identity:.3f}' if mean_v_identity is not None else 'n/a'}, "
            f"species {','.join(species) or '?'}",
            flush=True,
        )

    fields = [
        "clone_id",
        "receptor_type",
        "mean_peptide_content",
        "vh_peptide_content",
        "vl_peptide_content",
        "vh_v_gene",
        "vh_v_identity",
        "vl_v_gene",
        "vl_v_identity",
        "mean_v_identity",
        "germline_species",
        "n_non_human_peptides",
        "status",
        "reason",
    ]
    with open("./outputs/humanness.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    peptide_fields = [
        "clone_id",
        "chain",
        "peptide",
        "position",
        "start_offset",
        "regions",
        "in_cdr",
        "nearest_human_peptide",
    ]
    with open("./outputs/peptide_detail.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=peptide_fields, extrasaction="ignore")
        writer.writeheader()
        for row in peptide_rows:
            writer.writerow(row)

    payload = {
        "settings": settings,
        "clones": detail,
        "failures": failures,
        "n_scored": len(detail),
        "score_definition": (
            f"peptide_content is the fraction of the domain's {peptide_length}-mers present in "
            f"the '{database_name}' reference. With human-oas the reference is antibody peptides "
            "observed in at least 10% of human subjects, making this an OASis-like score."
        ),
        "caveat": (
            "Humanness is not immunogenicity: no MHC-II binding or T-cell epitope prediction is "
            "performed anywhere in this workflow. The score is also germline-biased, because "
            "both the reference repertoire and the identity baseline are built from observed "
            "antibodies, so an unusual but perfectly good clone scores worse than a common "
            "mediocre one."
        ),
    }
    with open("./outputs/humanness.json", "w") as handle:
        json.dump(payload, handle, indent=2)

    if not detail:
        fail("No clone could be scored for humanness. See humanness.csv for the reasons.")

    print(
        f"Scored {len(detail)} clone(s); {len(peptide_rows)} non-human peptide(s) recorded",
        flush=True,
    )
    print(
        "Humanness is a proxy for immunogenicity risk, not a T-cell epitope prediction.",
        flush=True,
    )


if __name__ == "__main__":
    main()
