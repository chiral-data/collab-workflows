#!/usr/bin/env python3
"""Node 02 — Fetch MSA.

Retrieve the multiple sequence alignment for every protein chain in the system and write it
to ./outputs/msa/<chain>.a3m, so node 03 can load the alignment from disk and run with no
network access at all.

Why this is a node of its own: the MSA is the only step that needs the public internet, and
it is the step most likely to fail for reasons unrelated to the model — a rate limit, a
server outage, a sequence the server cannot align. Doing it here means such a failure costs
seconds of CPU instead of a GPU job, and it makes the alignment an inspectable artifact:
msa_summary.json carries the depth per chain, which is the first thing to check when an
ensemble comes out suspiciously narrow, because shallow MSA is upstream's stated first cause
of degraded quality.

If the user drops their own .a3m files into ./inputs/, they are used verbatim and no server
is contacted. That is also the polite path for repeated runs — upstream explicitly asks users
not to send bulk workloads to the public ColabFold server.
"""
import glob
import json
import os
import shutil
import sys
from pathlib import Path

A3M_DIR = "./outputs/msa"


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def unpaired_a3m_paths(directory):
    """The per-chain unpaired A3M files, excluding the paired ones.

    A plain `*.a3m` glob also matches `<chain>.paired.a3m`, whose stem is `<chain>.paired` and
    is not a chain id. bioemu2's `load_a3m_msa` requires the mapping's keys to equal the
    system's protein chain ids exactly, so including the paired files would make it raise.
    This one helper is the only place that distinction is encoded.
    """
    return sorted(
        path
        for path in glob.glob(os.path.join(directory, "*.a3m"))
        if not path.endswith(".paired.a3m")
    )


def count_a3m_rows(path):
    """Number of records in an A3M file: every line beginning with '>'."""
    rows = 0
    try:
        with open(path) as handle:
            for line in handle:
                if line.startswith(">"):
                    rows += 1
    except Exception as exc:
        fail(f"Could not read A3M '{path}': {exc}")
    return rows


def load_inventory():
    path = "./inputs/system_inventory.json"
    if not os.path.exists(path):
        return {}
    try:
        with open(path) as handle:
            return json.load(handle)
    except Exception as exc:
        print(f"WARNING: could not read system inventory: {exc}", flush=True)
        return {}


def user_supplied_a3m(protein_chains):
    """Map chain id -> uploaded A3M path, if the user supplied any.

    A single uploaded A3M is taken to be the alignment for a single-protein-chain system.
    For several chains, files must be named <chain>.a3m.
    """
    uploads = unpaired_a3m_paths("./inputs")
    if not uploads:
        return {}
    if len(uploads) == 1 and len(protein_chains) == 1:
        return {protein_chains[0]: uploads[0]}

    mapping = {}
    for path in uploads:
        chain = os.path.splitext(os.path.basename(path))[0]
        if chain in protein_chains:
            mapping[chain] = path
    missing = [chain for chain in protein_chains if chain not in mapping]
    if missing:
        fail(
            f"Uploaded A3M files do not cover every protein chain. Missing: {missing}. "
            f"Name each file <chain_id>.a3m (chains present: {protein_chains}), or upload "
            "none and let the ColabFold server be queried."
        )
    return mapping


def main():
    os.makedirs("./outputs", exist_ok=True)
    os.makedirs(A3M_DIR, exist_ok=True)

    msa_depth = int(os.environ.get("PARAM_MSA_DEPTH", "0"))
    use_env = os.environ.get("PARAM_USE_ENV", "true").lower() == "true"

    system_path = "./inputs/system.yaml"
    if not os.path.exists(system_path):
        fail(f"System YAML not found: {system_path}")

    # Pass the system through unchanged so node 03 gets it from one upstream node.
    shutil.copyfile(system_path, "./outputs/system.yaml")
    # The inventory travels with it. Node 03 needs the token count to advise on batch size and
    # the report needs the chain list to say what was sampled; both would otherwise be
    # recomputed downstream from the same file, which is the defect class CLAUDE.md names.
    if os.path.exists("./inputs/system_inventory.json"):
        shutil.copyfile("./inputs/system_inventory.json", "./outputs/system_inventory.json")
    else:
        print(
            "WARNING: node 01's system_inventory.json was not staged here, so the token count "
            "and chain list will be absent downstream.",
            flush=True,
        )

    if msa_depth < 0:
        fail(f"msa_depth must be 0 (full alignment) or positive, got {msa_depth}")

    from bioemu2 import System

    try:
        system = System.from_yaml(Path(system_path))
    except Exception as exc:
        fail(f"Could not load system YAML '{system_path}': {exc}")

    inventory = load_inventory()
    protein_chains = list(inventory.get("protein_chains") or [])
    if not protein_chains:
        # Fall back to the system itself rather than trusting the inventory alone.
        for entity in getattr(system, "entities", []) or []:
            entity_type = str(
                getattr(entity, "entity_type", None) or type(entity).__name__
            ).lower()
            if "protein" in entity_type:
                protein_chains.extend(list(getattr(entity, "asym_ids", []) or []))

    summary = {
        "protein_chains": protein_chains,
        "msa_depth_param": msa_depth,
        "use_env": use_env,
        "source": None,
        "chains": {},
        "note": (
            "BioEmu2 ensemble quality tracks alignment depth; upstream names shallow or poor "
            "MSAs as the first cause of degraded quality. Check these depths before "
            "believing a narrow ensemble."
        ),
    }

    if not protein_chains:
        print(
            "No protein chain in this system — BioEmu2 computes no MSA. Writing an empty "
            "alignment set.",
            flush=True,
        )
        summary["source"] = "none (no protein chain)"
        # msa/ is declared as a directory and was created above, so it ships empty rather
        # than missing. The empty-but-present directory is not a fallback signal on its own:
        # node 03 decides whether an alignment was due from protein_chains in this summary,
        # which is why that key is written on every path including this one.
        with open("./outputs/msa_summary.json", "w") as handle:
            json.dump(summary, handle, indent=2)
        return

    uploaded = user_supplied_a3m(protein_chains)

    if uploaded:
        print(f"Using {len(uploaded)} user-supplied alignment(s); no server contacted.", flush=True)
        summary["source"] = "user-supplied A3M"
        # msa_depth subsamples what comes back from the server. An uploaded alignment is used
        # verbatim, so the parameter does nothing here — say so rather than leaving a run that
        # looks subsampled and is not.
        if msa_depth > 0:
            print(
                f"NOTE: msa_depth={msa_depth} is ignored for a user-supplied alignment, which "
                "is used exactly as given. Subsample the file before uploading it if you want "
                "a shallower alignment.",
                flush=True,
            )
            summary["msa_depth_applied"] = False
        for chain, path in uploaded.items():
            dest = os.path.join(A3M_DIR, f"{chain}.a3m")
            shutil.copyfile(path, dest)
            summary["chains"][chain] = {
                "a3m": f"msa/{chain}.a3m",
                "rows": count_a3m_rows(dest),
                "subsampled": False,
                # An uploaded set is unpaired only. For a heteromer this gives up the paired
                # feature the model would otherwise get, which is a real cost of the upload
                # path and is recorded rather than hidden.
                "paired_a3m": None,
            }
        summary["paired"] = False
    else:
        from bioemu2 import get_msa

        print(
            f"Querying the public ColabFold server for {len(protein_chains)} protein "
            f"chain(s) (use_env={use_env}). This needs network access.",
            flush=True,
        )
        try:
            msa = get_msa(system, use_env=use_env)
        except Exception as exc:
            fail(
                f"MSA retrieval failed: {exc}. The public ColabFold server may be "
                "unreachable or rate-limiting. Retry later, or supply your own alignment by "
                "uploading <chain>.a3m files into the workflow inputs."
            )
        summary["source"] = "ColabFold server"
        # A heteromer gets two alignments per chain from ColabFold: the independent (unpaired)
        # one and a positionally paired one. The model uses both — bioemu2's pairmixer is built
        # with use_paired_feature=True, and prepare_structure() substitutes zeros when the
        # paired side is None — so writing only the unpaired half silently degrades every
        # multi-sequence complex. get_msa returns paired alignments only when the system has
        # more than one unique protein sequence, so a monomer or homomer legitimately has none.
        paired_available = getattr(msa, "paired_a3m", None) is not None
        summary["paired"] = paired_available
        if paired_available:
            print(
                "Heteromer: writing the positionally paired alignment alongside the unpaired "
                "one, so node 03 can hand the model both.",
                flush=True,
            )
        else:
            print(
                "No paired alignment: the system has one unique protein sequence, so BioEmu2 "
                "needs only the unpaired alignment.",
                flush=True,
            )

        subsampled = False
        if msa_depth > 0:
            try:
                msa = msa[:msa_depth]
                subsampled = True
                print(f"Subsampled the alignment to {msa_depth} row(s)", flush=True)
            except Exception as exc:
                print(f"WARNING: could not subsample the alignment: {exc}", flush=True)

        for chain in protein_chains:
            try:
                a3m_text = msa.a3m(chain)
            except Exception as exc:
                fail(f"Could not extract the alignment for chain '{chain}': {exc}")
            dest = os.path.join(A3M_DIR, f"{chain}.a3m")
            with open(dest, "w") as handle:
                handle.write(a3m_text)
            summary["chains"][chain] = {
                "a3m": f"msa/{chain}.a3m",
                "rows": count_a3m_rows(dest),
                "subsampled": subsampled,
                "paired_a3m": None,
            }
            if paired_available:
                # Subsampling happened on the MSA object above, which applies one set of row
                # indices to every chain and to both alignments, so the paired rows stay
                # positionally aligned with each other. Extracting the text per chain
                # afterwards cannot break that; extracting first and subsampling per file
                # would have.
                try:
                    paired_text = msa.a3m(chain, paired=True)
                except Exception as exc:
                    fail(
                        f"The alignment reported a paired set but chain '{chain}' could not be "
                        f"extracted from it: {exc}"
                    )
                paired_dest = os.path.join(A3M_DIR, f"{chain}.paired.a3m")
                with open(paired_dest, "w") as handle:
                    handle.write(paired_text)
                summary["chains"][chain]["paired_a3m"] = f"msa/{chain}.paired.a3m"
                summary["chains"][chain]["paired_rows"] = count_a3m_rows(paired_dest)
        summary["msa_depth_applied"] = bool(msa_depth > 0 and subsampled)

    written = sorted(unpaired_a3m_paths(A3M_DIR))
    if not written:
        fail("No alignment files were written, but the system has protein chains.")
    if len(written) != len(protein_chains):
        fail(
            f"Wrote {len(written)} unpaired alignment(s) for {len(protein_chains)} protein "
            "chain(s). bioemu2's load_a3m_msa requires exactly one per chain and rejects a "
            "mismatch, so this would fail in node 03 instead."
        )

    depths = {chain: info["rows"] for chain, info in summary["chains"].items()}
    summary["min_depth"] = min(depths.values()) if depths else 0
    summary["max_depth"] = max(depths.values()) if depths else 0

    with open("./outputs/msa_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)

    print(f"Source: {summary['source']}", flush=True)
    for chain, rows in sorted(depths.items()):
        print(f"  chain {chain}: {rows} alignment row(s)", flush=True)
    if summary["min_depth"] < 32:
        print(
            f"WARNING: shallowest alignment has only {summary['min_depth']} row(s). "
            "BioEmu2 quality degrades on shallow alignments (orphan and some viral "
            "proteins); treat the resulting ensemble with caution.",
            flush=True,
        )


if __name__ == "__main__":
    main()
