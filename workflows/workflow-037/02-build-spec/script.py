#!/usr/bin/env python3
"""Node 02 — Build Condition Spec.

Render the task-specific RFD4 condition spec (input / contig / conditions, with
C_CRDGRP, C_SEQ, C_HOT, C_CRD and sel() selectors) from the workflow parameters, then
validate it on CPU before any GPU time is spent.

This node exists because sel() and contig errors are the most common RFD4 failure mode.
Upstream's scripts/write_condition_spec_cifs.py renders a spec to an annotated CIF without
running the model, which turns a 20-minute GPU failure into a 10-second CPU failure — and
spec_preview.cif is itself useful, since the user can open it and see exactly which atoms
are conditioned. This node always writes its own annotation of the conditioned residues,
and additionally runs the upstream renderer over the spec when it can be located in the
image, letting its richer output replace ours.

All residue ids written into the spec are mmCIF label_seq_id, per the RFD4 docs. Parameters
given in author numbering (the default) are translated here using node 01's inventory.
"""
import json
import os
import re
import shutil
import subprocess
import sys

TASKS = ("binder", "motif_scaffold", "monomer", "binder_redesign")
LENGTH_RE = re.compile(r"^\d+(-\d+)?$")
# Candidate locations for the upstream spec renderer inside the image.
UPSTREAM_RENDERER_CANDIDATES = (
    "/opt/rfd4/scripts/write_condition_spec_cifs.py",
    "/opt/rfd4-proteina/scripts/write_condition_spec_cifs.py",
    "/usr/local/scripts/write_condition_spec_cifs.py",
)


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def load_inventory(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except Exception as exc:
        fail(f"Could not read inventory '{path}': {exc}")


def chain_by_id(inventory, chain_id):
    """Resolve a user-supplied chain id against either numbering."""
    for chain in inventory["chains"]:
        if chain_id in (chain["label_asym_id"], chain["auth_asym_id"]):
            return chain
    present = ", ".join(
        f"{c['label_asym_id']}(auth {c['auth_asym_id']})" for c in inventory["chains"]
    )
    fail(f"Chain '{chain_id}' is not in the target structure. Chains present: {present}")


def parse_residue_tokens(spec, default_chain):
    """Parse 'A45,A46-50,82' into [(chain, residue_id), ...] as written by the user."""
    out = []
    for raw in spec.split(","):
        token = raw.strip()
        if not token:
            continue
        match = re.match(r"^([A-Za-z]*)(\d+)(?:-(\d+))?$", token)
        if not match:
            fail(
                f"Could not parse residue token '{token}'. Expected forms: A45, A45-50, 45."
            )
        chain = match.group(1) or default_chain
        start = int(match.group(2))
        end = int(match.group(3)) if match.group(3) else start
        if end < start:
            fail(f"Residue range '{token}' runs backwards.")
        for residue in range(start, end + 1):
            out.append((chain, residue))
    return out


def to_label_ids(residues, inventory, numbering):
    """Translate (chain, residue) pairs to (label_asym_id, label_seq_id)."""
    translated = []
    for chain_id, residue_id in residues:
        chain = chain_by_id(inventory, chain_id)
        if numbering == "label":
            if str(residue_id) not in chain["label_to_auth"]:
                fail(
                    f"Residue {chain_id}{residue_id} (label numbering) is not in the "
                    f"structure. Chain {chain['label_asym_id']} covers label "
                    f"{chain['label_seq_range'][0]}-{chain['label_seq_range'][1]}."
                )
            label_seq = int(residue_id)
        else:
            label_seq = chain["auth_to_label"].get(str(residue_id))
            if label_seq is None:
                fail(
                    f"Residue {chain_id}{residue_id} (author numbering) is not in the "
                    f"structure. Chain {chain['label_asym_id']} covers auth "
                    f"{chain['auth_seq_range'][0]}-{chain['auth_seq_range'][1]}. "
                    "If you meant mmCIF numbering, set residue_numbering=label."
                )
        translated.append((chain["label_asym_id"], label_seq))
    return translated


def contiguous_runs(pairs):
    """Collapse [(chain, seq), ...] into [(chain, start, end), ...]."""
    runs = []
    for chain, seq in sorted(set(pairs)):
        if runs and runs[-1][0] == chain and seq == runs[-1][2] + 1:
            runs[-1] = (chain, runs[-1][1], seq)
        else:
            runs.append((chain, seq, seq))
    return runs


def sel(chain="", comp="", res="", atom=""):
    """Build an RFD4 selector. Documented field order: chain/res_name/res_id/atom_name.

    Empty fields are wildcards. The exact wildcard spelling is the one upstream-documented
    unknown in this node — confirm against the RFD4 docs on the first in-image run.
    """
    return f"sel('{chain}/{comp}/{res}/{atom}')"


def chain_segment(chain):
    """Contig segment spanning a chain, warning if it spans residues absent from the file.

    An unresolved loop leaves a hole in label_seq_id, so a span like A1-129 can cover
    residues that have no coordinates. RFD4 is conditioned on the structure, so what it
    does with a hole is its business — but this is worth seeing before the GPU runs, not
    after a design comes out with a gap in the middle of the target.
    """
    lo, hi = chain["label_seq_range"]
    present = {residue["label_seq_id"] for residue in chain["residues"]}
    missing = [seq for seq in range(lo, hi + 1) if seq not in present]
    if missing:
        runs = []
        for seq in missing:
            if runs and seq == runs[-1][1] + 1:
                runs[-1] = (runs[-1][0], seq)
            else:
                runs.append((seq, seq))
        gaps = ", ".join(f"{a}" if a == b else f"{a}-{b}" for a, b in runs)
        print(
            f"WARNING: chain {chain['label_asym_id']} is missing {len(missing)} residue(s) "
            f"inside its span (label {gaps}) — unresolved loops or non-standard residues. "
            f"The contig segment {chain['label_asym_id']}{lo}-{hi} covers them.",
            flush=True,
        )
    return f"{chain['label_asym_id']}{lo}-{hi}"


def build_binder(params, inventory):
    target = chain_by_id(inventory, params["target_chain"])
    contig = f"{params['binder_length']}, /0, {chain_segment(target)}"
    conditions = {}
    if params["hotspot_residues"]:
        hotspots = to_label_ids(
            parse_residue_tokens(params["hotspot_residues"], target["label_asym_id"]),
            inventory,
            params["residue_numbering"],
        )
        conditions["C_HOT"] = [
            sel(chain=chain, res=",".join(str(s) for _, s in group))
            for chain, group in group_by_chain(hotspots)
        ]
    return contig, conditions, [target["label_asym_id"]]


def build_monomer(params, inventory):
    return params["binder_length"], {}, []


def build_motif_scaffold(params, inventory):
    if not params["motif_residues"]:
        fail("task=motif_scaffold requires motif_residues.")
    target = chain_by_id(inventory, params["target_chain"])
    motif = to_label_ids(
        parse_residue_tokens(params["motif_residues"], target["label_asym_id"]),
        inventory,
        params["residue_numbering"],
    )
    runs = contiguous_runs(motif)
    flank = params["binder_length"]
    print(
        f"NOTE: motif_scaffold uses binder_length='{flank}' as the length range of each "
        f"scaffold flank around {len(runs)} motif segment(s).",
        flush=True,
    )
    segments = [flank]
    for chain, start, end in runs:
        segments.append(f"{chain}{start}-{end}")
        segments.append(flank)
    contig = ", ".join(segments)
    res_by_chain = group_by_chain(motif)
    conditions = {
        # Hold the motif's geometry as one rigid group, and its identity as sequence.
        "C_CRDGRP": [
            sel(chain=chain, res=",".join(str(s) for _, s in group))
            for chain, group in res_by_chain
        ],
        "C_SEQ": [
            sel(chain=chain, res=",".join(str(s) for _, s in group))
            for chain, group in res_by_chain
        ],
    }
    return contig, conditions, [chain for chain, _ in res_by_chain]


def build_binder_redesign(params, inventory):
    target = chain_by_id(inventory, params["target_chain"])
    redesign = chain_by_id(inventory, params["redesign_chain"])
    if redesign["label_asym_id"] == target["label_asym_id"]:
        fail("redesign_chain and target_chain must differ for task=binder_redesign.")
    print(
        f"NOTE: regenerating chain {redesign['label_asym_id']} "
        f"({redesign['n_residues']} residues in the input) at length "
        f"'{params['binder_length']}' — set binder_length to match the existing binder "
        "unless you intend to change its length.",
        flush=True,
    )
    contig = f"{params['binder_length']}, /0, {chain_segment(target)}"
    conditions = {"C_CRD": [sel(chain=target["label_asym_id"])]}
    if params["hotspot_residues"]:
        hotspots = to_label_ids(
            parse_residue_tokens(params["hotspot_residues"], target["label_asym_id"]),
            inventory,
            params["residue_numbering"],
        )
        conditions["C_HOT"] = [
            sel(chain=chain, res=",".join(str(s) for _, s in group))
            for chain, group in group_by_chain(hotspots)
        ]
    return contig, conditions, [target["label_asym_id"]]


def group_by_chain(pairs):
    grouped = {}
    for chain, seq in sorted(set(pairs)):
        grouped.setdefault(chain, []).append((chain, seq))
    return list(grouped.items())


BUILDERS = {
    "binder": build_binder,
    "monomer": build_monomer,
    "motif_scaffold": build_motif_scaffold,
    "binder_redesign": build_binder_redesign,
}


def render_spec(task, context):
    from jinja2 import Environment, FileSystemLoader, StrictUndefined

    env = Environment(
        loader=FileSystemLoader("templates"),
        undefined=StrictUndefined,
        keep_trailing_newline=True,
        trim_blocks=False,
        lstrip_blocks=False,
    )
    template = env.get_template(f"{task}.yaml.j2")
    return template.render(**context)


def write_annotation(target_cif, conditioned, dest):
    """Copy the target and append our own annotation of the conditioned residues.

    A whole-chain condition such as C_CRD sel('A///') is recorded with label_seq_id '.',
    so chain-level conditioning is visible rather than missing. An empty loop_ is invalid
    mmCIF, so with nothing to annotate we write a comment instead of an empty table.
    """
    shutil.copyfile(target_cif, dest)
    with open(dest, "rb") as handle:
        handle.seek(-1, os.SEEK_END)
        needs_newline = handle.read(1) != b"\n"
    with open(dest, "a") as handle:
        if needs_newline:
            handle.write("\n")
        handle.write("#\n# Conditioned residues (annotation written by workflow-037 node 02)\n")
        if not conditioned:
            handle.write("# none: this spec sets no conditions\n#\n")
            return
        handle.write("loop_\n_rfd4_spec_preview.condition\n")
        handle.write("_rfd4_spec_preview.label_asym_id\n_rfd4_spec_preview.label_seq_id\n")
        for condition, chain, seq in conditioned:
            handle.write(f"{condition} {chain or '.'} {seq or '.'}\n")
        handle.write("#\n")


def run_upstream_renderer(spec_path, dest):
    """If the upstream renderer is in the image, let its output replace ours."""
    script = os.environ.get("RFD4_SPEC_CIF_SCRIPT", "")
    candidates = [script] if script else list(UPSTREAM_RENDERER_CANDIDATES)
    found = next((c for c in candidates if c and os.path.exists(c)), None)
    if not found:
        print(
            "WARNING: upstream write_condition_spec_cifs.py not found in this image; "
            "spec_preview.cif carries this node's own annotation only. Set "
            "RFD4_SPEC_CIF_SCRIPT to its path to enable the upstream validation pass.",
            flush=True,
        )
        return False
    print(f"Validating spec with {found}", flush=True)
    result = subprocess.run(
        [sys.executable, found, spec_path, "--out", dest],
        capture_output=True,
        text=True,
    )
    if result.stdout:
        print(result.stdout.strip(), flush=True)
    if result.returncode != 0:
        print(result.stderr.strip(), flush=True)
        fail(
            "Condition spec failed upstream validation (see above). Fix the contig or "
            "the hotspot/motif residues; no GPU time has been spent."
        )
    return True


def main():
    os.makedirs("./outputs", exist_ok=True)

    params = {
        "task": os.environ.get("PARAM_TASK", "binder"),
        "target_chain": os.environ.get("PARAM_TARGET_CHAIN", "A"),
        "binder_length": os.environ.get("PARAM_BINDER_LENGTH", "50-120"),
        "contig": os.environ.get("PARAM_CONTIG", ""),
        "hotspot_residues": os.environ.get("PARAM_HOTSPOT_RESIDUES", ""),
        "motif_residues": os.environ.get("PARAM_MOTIF_RESIDUES", ""),
        "redesign_chain": os.environ.get("PARAM_REDESIGN_CHAIN", "B"),
        "residue_numbering": os.environ.get("PARAM_RESIDUE_NUMBERING", "auth"),
    }

    if params["task"] not in TASKS:
        fail(f"Unknown task '{params['task']}'. Expected one of: {', '.join(TASKS)}")
    if params["residue_numbering"] not in ("auth", "label"):
        fail(
            f"residue_numbering must be 'auth' or 'label', got "
            f"'{params['residue_numbering']}'"
        )
    if not LENGTH_RE.match(params["binder_length"]):
        fail(
            f"binder_length '{params['binder_length']}' is not a length or range, "
            "e.g. 80 or 50-120."
        )

    target_cif = "./inputs/target.cif"
    inventory_path = "./inputs/target_inventory.json"
    for path in (target_cif, inventory_path):
        if not os.path.exists(path):
            fail(f"Required input not found: {path}")
    inventory = load_inventory(inventory_path)

    contig, conditions, conditioned_chains = BUILDERS[params["task"]](params, inventory)

    if params["contig"]:
        print(
            f"NOTE: contig override in use; templated contig '{contig}' discarded.",
            flush=True,
        )
        contig = params["contig"]

    needs_input = params["task"] != "monomer"
    spec_text = render_spec(
        params["task"],
        {
            "input_cif": "target.cif" if needs_input else "",
            "contig": contig,
            "conditions": conditions,
            "task": params["task"],
            "residue_numbering": params["residue_numbering"],
        },
    )
    spec_path = "./outputs/design_spec.yaml"
    with open(spec_path, "w") as handle:
        handle.write(spec_text)

    print(f"Task: {params['task']} (numbering: {params['residue_numbering']})", flush=True)
    print(f"Contig: {contig}", flush=True)
    for name, selectors in conditions.items():
        print(f"Condition {name}: {' '.join(selectors)}", flush=True)
    if not conditions:
        print("Conditions: none", flush=True)

    conditioned = []
    for name, selectors in conditions.items():
        for selector in selectors:
            match = re.match(r"^sel\('([^/]*)/[^/]*/([^/]*)/", selector)
            if not match:
                continue
            chain = match.group(1)
            residues = [s for s in match.group(2).split(",") if s]
            if not residues:
                conditioned.append((name, chain, "."))  # whole-chain condition
                continue
            for seq in residues:
                conditioned.append((name, chain, seq))
    preview_path = "./outputs/spec_preview.cif"
    if needs_input:
        write_annotation(target_cif, conditioned, preview_path)
        run_upstream_renderer(spec_path, preview_path)
    else:
        with open(preview_path, "w") as handle:
            handle.write("# task=monomer: no input structure, nothing to annotate\n")
    # Always passed through, so node 03's declared inputs hold for every task.
    shutil.copyfile(target_cif, "./outputs/target.cif")

    print(f"Wrote {spec_path} and {preview_path}", flush=True)


if __name__ == "__main__":
    main()
