#!/usr/bin/env python3
"""Node 02 — Build Condition Spec.

Render the task-specific RFD4 condition spec from the workflow parameters, then validate it
on CPU before any GPU time is spent.

The spec format is upstream's, confirmed against the RFD4-Proteina 0.0.3 README, the
condition-spec guide and the cookbook:

    design:                                     # an example id, not a bare mapping
      input: target.cif
      contig: "50-120, /0, A17-131"
      conditions:
        C_CRDGRP: {1: [{select: "sel('A/*/*/[N,CA,C,O]')"}]}
        C_SEQ:    {true: [{select: "sel('A')"}]}
        C_HOT:    {true: [{select: "sel('A/*/56')"}]}

Conditions are keyed by the condition *value* (a group id for C_CRDGRP, true for the
boolean conditions), each holding a list of {select: <query>} entries. Selections are
atomworks query strings, where sel() takes chain/res_name/res_id/atom_name with `*` as the
wildcard and [A,B] for lists.

Per the upstream guide, res_id in both contig and selections resolves against the input's
mmCIF label_seq_id, which is why node 01 emits both numberings and this node translates.

This node exists because selection and contig errors are the most common RFD4 failure mode.
It validates by calling resolve_condition_spec_inputs(), the same code path generation uses
to turn a spec into an annotated CIF, so a bad spec fails here in seconds instead of 20
minutes into a GPU job.
"""
import json
import os
import re
import shutil
import sys
import tempfile

TASKS = ("binder", "motif_scaffold", "monomer", "binder_redesign")
LENGTH_RE = re.compile(r"^\d+(-\d+)?$")
EXAMPLE_ID = "design"
BACKBONE_ATOMS = "[N,CA,C,O]"


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
            fail(f"Could not parse residue token '{token}'. Expected forms: A45, A45-50, 45.")
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
            label_seq = int(label_seq)
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


def sel(chain, res_name=None, res_id=None, atoms=None):
    """Build an atomworks selection string in sel() form.

    Fields are chain/res_name/res_id/atom_name with `*` as the wildcard, as in upstream's
    sel('A/*/*/[N,CA,C,O]'). Trailing fields are omitted rather than left empty, matching
    upstream's sel('A').
    """
    parts = [chain]
    for value in (res_name, res_id, atoms):
        if value is None:
            break
        parts.append(str(value))
    return f"sel('{'/'.join(parts)}')"


def select_entries(queries):
    return [{"select": query} for query in queries]


def hold_target(chain_id):
    """Hold a chain's backbone coordinates as one rigid group and freeze its sequence.

    This is upstream's binder recipe: without it the "target" is not actually held.
    """
    return {
        "C_CRDGRP": {1: select_entries([sel(chain_id, "*", "*", BACKBONE_ATOMS)])},
        "C_SEQ": {True: select_entries([sel(chain_id)])},
    }


def hotspot_conditions(hotspots):
    """One {select: ...} entry per hotspot residue, as in upstream's cookbook."""
    return {
        "C_HOT": {
            True: select_entries([sel(chain, "*", seq) for chain, seq in hotspots])
        }
    }


def chain_segment(chain):
    """Contig segment spanning a chain, warning if it spans residues absent from the file."""
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


def resolve_hotspots(params, inventory, target):
    if not params["hotspot_residues"]:
        return {}
    hotspots = to_label_ids(
        parse_residue_tokens(params["hotspot_residues"], target["label_asym_id"]),
        inventory,
        params["residue_numbering"],
    )
    return hotspot_conditions(hotspots)


def build_binder(params, inventory):
    target = chain_by_id(inventory, params["target_chain"])
    contig = f"{params['binder_length']}, /0, {chain_segment(target)}"
    conditions = hold_target(target["label_asym_id"])
    conditions.update(resolve_hotspots(params, inventory, target))
    return contig, conditions, False


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
    conditions = hold_target(target["label_asym_id"])
    conditions.update(resolve_hotspots(params, inventory, target))
    return contig, conditions, False


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
    residues = [f"{seq}" for _, seq in sorted(set(motif))]
    by_chain = {}
    for chain, seq in sorted(set(motif)):
        by_chain.setdefault(chain, []).append(seq)
    conditions = {
        # Hold the motif's backbone as one rigid group and freeze its identity.
        "C_CRDGRP": {
            1: select_entries(
                [
                    sel(chain, "*", f"[{','.join(str(s) for s in seqs)}]", BACKBONE_ATOMS)
                    for chain, seqs in by_chain.items()
                ]
            )
        },
        "C_SEQ": {
            True: select_entries(
                [
                    sel(chain, "*", f"[{','.join(str(s) for s in seqs)}]")
                    for chain, seqs in by_chain.items()
                ]
            )
        },
    }
    return contig, conditions, False


MONOMER_SCAFFOLD_NAME = "target.cif"  # the scaffold replaces the target for task=monomer

MONOMER_SCAFFOLD_TEMPLATE = """data_{name}
#
_entry.id       {name}
_entry.author   workflow-037
#
_entity.id                 0
_entity.type               polymer
_entity.pdbx_description   .
#
_entity_poly.entity_id                      0
_entity_poly.type                           polypeptide(L)
_entity_poly.nstd_linkage                   no
_entity_poly.nstd_monomer                   yes
_entity_poly.pdbx_seq_one_letter_code       X
_entity_poly.pdbx_seq_one_letter_code_can   X
_entity_poly.pdbx_strand_id                 A
_entity_poly.pdbx_target_identifier         ?
#
_entity_poly_seq.entity_id   0
_entity_poly_seq.hetero      n
_entity_poly_seq.mon_id      {comp}
_entity_poly_seq.num         1
#
_chem_comp.id                        {comp}
_chem_comp.name                      {comp}
_chem_comp.type                      NON-POLYMER
_chem_comp.formula                   ?
_chem_comp.mon_nstd_parent_comp_id   ?
_chem_comp.pdbx_synonyms             ?
_chem_comp.formula_weight            ?
#
_atom_site.group_PDB                         HETATM
_atom_site.type_symbol                       X
_atom_site.label_atom_id                     UNK
_atom_site.label_alt_id                      .
_atom_site.label_comp_id                     {comp}
_atom_site.label_asym_id                     A
_atom_site.label_entity_id                   0
_atom_site.label_seq_id                      1
_atom_site.pdbx_PDB_ins_code                 .
_atom_site.auth_seq_id                       1
_atom_site.auth_comp_id                      {comp}
_atom_site.auth_asym_id                      A
_atom_site.auth_atom_id                      UNK
_atom_site.B_iso_or_equiv                    nan
_atom_site.occupancy                         0.0
_atom_site.pdbx_formal_charge                ?
_atom_site.annotation_expsegmax_1_residue    {high}
_atom_site.mask_expsegmin_1_residue          True
_atom_site.mask_expsegmax_1_residue          True
_atom_site.annotation_expsegmin_1_residue    {low}
_atom_site.condition_hbond_1_atom            0
_atom_site.condition_index-group_1_token     0
_atom_site.mask_index_1_residue              True
_atom_site.mask_coordinate_1_atom            False
_atom_site.mask_sequence_1_residue           False
_atom_site.annotation_src-component_1_atom   .
_atom_site.atom_array_index                  0
_atom_site.Cartn_x                           nan
_atom_site.Cartn_y                           nan
_atom_site.Cartn_z                           nan
_atom_site.pdbx_PDB_model_num                1
_atom_site.id                                1
#
_rfd_design.sampled_contig   {comp}
#
_chem_comp_atom.comp_id                    {comp}
_chem_comp_atom.atom_id                    UNK
_chem_comp_atom.type_symbol                X
_chem_comp_atom.charge                     0
_chem_comp_atom.pdbx_aromatic_flag         N
_chem_comp_atom.pdbx_leaving_atom_flag     N
_chem_comp_atom.pdbx_model_Cartn_x_ideal   nan
_chem_comp_atom.pdbx_model_Cartn_y_ideal   nan
_chem_comp_atom.pdbx_model_Cartn_z_ideal   nan
#
"""


def write_monomer_scaffold(low, high, dest):
    """Write an annotated CIF holding one expandable chain and no real atoms.

    A contig of nothing but expandable segments is rejected upstream ("input exists to be
    referenced by a preserved segment"), so a monomer cannot be asked for with a contig.
    Upstream's answer is an annotated CIF carrying a single expandable segment and no
    coordinates: examples/monomer.cif in the model repository is 3 KB of annotations, one
    sentinel HETATM with NaN coordinates, and a component named for its length.

    This writes the same structure at the requested length. The format is copied from that
    file: the component id is "<low>-<high>", and the bounds appear again as
    annotation_expsegmin/max_1_residue. Upstream documents none of this, so node 02
    validates the result through RFD4's own parser before any GPU time is spent.
    """
    with open(dest, "w") as handle:
        handle.write(
            MONOMER_SCAFFOLD_TEMPLATE.format(
                name=MONOMER_SCAFFOLD_NAME, comp=f"{low}-{high}", low=low, high=high
            )
        )


def build_monomer(params, inventory):
    """Unconditioned single chain, designed from a scaffold written in place of the target.

    The scaffold replaces target.cif rather than travelling as a file of its own, so every
    task hands node 03 the same two files and the spec's `input:` never varies. main()
    writes it, because it owns the output copy step.
    """
    low, _, high = params["binder_length"].partition("-")
    high = high or low
    print(
        f"NOTE: task=monomer ignores the target structure. It designs a standalone chain "
        f"of {low}-{high} residues from a generated scaffold, written as target.cif, "
        "because a contig built only from expandable segments is rejected upstream.",
        flush=True,
    )
    return None, {}, True


BUILDERS = {
    "binder": build_binder,
    "motif_scaffold": build_motif_scaffold,
    "monomer": build_monomer,
    "binder_redesign": build_binder_redesign,
}


def check_contig(contig):
    """Upstream rejects a contig whose segments are all expandable."""
    segments = [s.strip() for s in contig.split(",") if s.strip() and s.strip() != "/0"]
    preserved = [s for s in segments if not LENGTH_RE.match(s)]
    if not preserved:
        fail(
            f"Contig '{contig}' has no preserved segment. Upstream rejects a contig built "
            "only from expandable lengths, since the input exists to be referenced."
        )


def dump_spec(spec):
    import yaml

    return yaml.safe_dump(spec, sort_keys=False, default_flow_style=False, width=100)


def validate_with_upstream(spec_path, dest):
    """Render the spec to an annotated CIF through the code generation itself uses."""
    try:
        from rfproteina.inference_utils.condition_spec import resolve_condition_spec_inputs
    except ImportError:
        # Package absent: the CPU-only development path. Degrade to our own annotation.
        print(
            "WARNING: rfproteina is not importable here, so the spec was not validated "
            "against upstream. spec_preview.cif carries this node's own annotation only.",
            flush=True,
        )
        return False
    except Exception as exc:
        # Package present but broken — e.g. Open Babel failing to load its plugins takes
        # the whole import chain down. Generation would fail the same way, so stop here.
        fail(
            f"rfproteina is installed but failed to import: {type(exc).__name__}: {exc}\n"
            "The image is broken, not the spec. Generation would fail the same way."
        )

    # Rendered CIFs go to a scratch directory, not ./outputs, so a validation failure
    # leaves no half-written directory among the declared outputs.
    out_dir = tempfile.mkdtemp(prefix="spec_cifs_")
    try:
        resolve_condition_spec_inputs(spec_path, out_dir, validate=True)
    except Exception as exc:
        fail(
            f"Condition spec rejected by RFD4: {type(exc).__name__}: {exc}\n"
            "Fix the contig or the selections; no GPU time has been spent."
        )
    rendered = sorted(
        os.path.join(out_dir, name) for name in os.listdir(out_dir) if name.endswith(".cif")
    )
    if not rendered:
        print("WARNING: upstream validation wrote no CIF; keeping this node's annotation.", flush=True)
        shutil.rmtree(out_dir, ignore_errors=True)
        return False
    shutil.copyfile(rendered[0], dest)
    shutil.rmtree(out_dir, ignore_errors=True)
    print(f"Spec validated by RFD4 and rendered to {dest}", flush=True)
    return True


def write_annotation(target_cif, conditions, dest):
    """Fallback annotation: the target plus a record of what each condition selects."""
    shutil.copyfile(target_cif, dest)
    with open(dest, "rb") as handle:
        handle.seek(-1, os.SEEK_END)
        needs_newline = handle.read(1) != b"\n"
    rows = [
        (name, value, entry["select"])
        for name, by_value in conditions.items()
        for value, entries in by_value.items()
        for entry in entries
    ]
    with open(dest, "a") as handle:
        if needs_newline:
            handle.write("\n")
        handle.write("#\n# Conditions (annotation written by workflow-037 node 02)\n")
        if not rows:
            handle.write("# none: this spec sets no conditions\n#\n")
            return
        for name, value, query in rows:
            handle.write(f"# {name} [{value}]: {query}\n")
        handle.write("#\n")


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

    contig, conditions, as_annotated_cif = BUILDERS[params["task"]](params, inventory)

    if params["contig"]:
        print(
            f"NOTE: contig override in use; templated contig '{contig}' discarded.",
            flush=True,
        )
        contig = params["contig"]
        as_annotated_cif = False
    if contig:
        check_contig(contig)

    example = {"input": "target.cif"}
    if as_annotated_cif:
        example["input_as_annotated_cif"] = True
    if contig:
        example["contig"] = contig
    if conditions:
        example["conditions"] = conditions

    spec_path = "./outputs/design_spec.yaml"
    with open(spec_path, "w") as handle:
        handle.write(f"# RFD4 condition spec — task: {params['task']}\n")
        handle.write(
            f"# res_id values are mmCIF label_seq_id "
            f"(parameters were given in: {params['residue_numbering']})\n"
        )
        handle.write(dump_spec({EXAMPLE_ID: example}))

    print(f"Task: {params['task']} (numbering: {params['residue_numbering']})", flush=True)
    print(f"Contig: {contig or 'none (annotated CIF passthrough)'}", flush=True)
    for name, by_value in conditions.items():
        for value, entries in by_value.items():
            queries = " ".join(entry["select"] for entry in entries)
            print(f"Condition {name}[{value}]: {queries}", flush=True)
    if not conditions:
        print("Conditions: none", flush=True)

    # The spec names its input as a bare relative filename, which upstream resolves against
    # the spec's own directory — so the input has to sit beside the spec before validation,
    # not after. Node 03 relies on the same side-by-side layout.
    if as_annotated_cif:
        low, _, high = params["binder_length"].partition("-")
        write_monomer_scaffold(int(low), int(high or low), "./outputs/target.cif")
    else:
        shutil.copyfile(target_cif, "./outputs/target.cif")

    preview_path = "./outputs/spec_preview.cif"
    write_annotation("./outputs/target.cif", conditions, preview_path)
    validate_with_upstream(spec_path, preview_path)

    print(f"Wrote {spec_path} and {preview_path}", flush=True)


if __name__ == "__main__":
    main()
