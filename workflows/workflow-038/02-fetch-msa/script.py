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
    uploads = sorted(glob.glob("./inputs/*.a3m"))
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
        # Keep the declared msa/*.a3m glob honest: there is genuinely nothing to write, and
        # node 03 treats an empty directory as "let BioEmu2 decide".
        with open("./outputs/msa_summary.json", "w") as handle:
            json.dump(summary, handle, indent=2)
        return

    uploaded = user_supplied_a3m(protein_chains)

    if uploaded:
        print(f"Using {len(uploaded)} user-supplied alignment(s); no server contacted.", flush=True)
        summary["source"] = "user-supplied A3M"
        for chain, path in uploaded.items():
            dest = os.path.join(A3M_DIR, f"{chain}.a3m")
            shutil.copyfile(path, dest)
            summary["chains"][chain] = {
                "a3m": f"msa/{chain}.a3m",
                "rows": count_a3m_rows(dest),
                "subsampled": False,
            }
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
            }

    written = sorted(glob.glob(os.path.join(A3M_DIR, "*.a3m")))
    if not written:
        fail("No alignment files were written, but the system has protein chains.")

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
