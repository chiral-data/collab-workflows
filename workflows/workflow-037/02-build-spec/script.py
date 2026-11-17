#!/usr/bin/env python3
"""Node 02 — Build Condition Spec.

Render the task-specific RFD4 condition spec from the workflow parameters, then validate it
on CPU before any GPU time is spent.

The spec format is upstream's, as given in the RFD4-Proteina 0.0.3 README, the
condition-spec guide and the cookbook:

    design:                                     # an example id, not a bare mapping
      input: target.cif
      contig: "50-120, /0, A17-131"
      conditions:
        C_CRDGRP: {1: [{select: "sel('A/*/*/[N,CA,C,O]')"}]}
        C_SEQ:    {true: [{select: "sel('A')"}]}
        C_HOT:    {true: [{select: "sel('A/*/56/[CG,OH]')"}]}

C_HOT is atom-level in all three of those sources — the cookbook's own example is
sel('A/*/56/[CG,OH]'), annotated "the PD-L1 side chain atoms to engage". An earlier
version of this node emitted sel('A/*/56') and claimed upstream confirmed it; that was
wrong, and the residue form marks every atom of the residue instead of a couple.

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
import math
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


def numeric_param(env_name, display_name, default, cast):
    """Read a numeric parameter, failing the same way every other bad parameter does.

    Casting inline at the point of use turns a typo into a raw ValueError traceback, which
    is the one error in this node that does not look like the others.
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
    # are caught by a range check: every comparison against nan is false, so a nan
    # cutoff passes "must be positive" and "outside the expected range" alike and then
    # makes every distance test fail, producing an empty interface rather than an error.
    # inf does the opposite and takes the whole chain.
    if not math.isfinite(value):
        fail(f"{display_name} must be a finite number, got '{raw}'.")
    return value


def load_inventory(path):
    try:
        with open(path) as handle:
            return json.load(handle)
    except Exception as exc:
        fail(f"Could not read inventory '{path}': {exc}")


def chain_by_id(inventory, chain_id):
    """Resolve a user-supplied chain id against either numbering, label first.

    An id can be the label_asym_id of one chain and the auth_asym_id of a different one.
    Node 01 resolves that case to the label match and says so; this has to agree with it,
    or the two nodes name different chains for the same parameter and the one the user is
    shown is not the one the design is built against.
    """
    for chain in inventory["chains"]:
        if chain_id == chain["label_asym_id"]:
            return chain
    for chain in inventory["chains"]:
        if chain_id == chain["auth_asym_id"]:
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


def insertion_code_siblings(chain, auth_residue_id):
    """Residues sharing this author number through an insertion code, label id and all.

    Node 01 already warns that author->label translation takes the first occurrence, and
    lists the label numbers of the residues it cannot name. The warning is repeated here
    because this is where the loss actually happens. Node 01 warns about the structure;
    the user names the residue in a different node, against a different parameter, and
    nothing in between mentions it again. On an antibody the shadowed residues are
    100A/100B/100C in CDR-H3 -- exactly the ones someone would name as hotspots.
    """
    prefix = str(auth_residue_id)
    siblings = []
    for label_seq, auth in chain.get("label_to_auth", {}).items():
        # Suffix must be alphabetic, or auth 1000 would read as a sibling of auth 100.
        suffix = auth[len(prefix):]
        if auth.startswith(prefix) and suffix and suffix.isalpha():
            siblings.append((auth, int(label_seq)))
    siblings.sort(key=lambda item: item[1])
    return siblings


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
            shadowed = insertion_code_siblings(chain, residue_id)
            if shadowed:
                names = ", ".join(name for name, _ in shadowed)
                labels = ", ".join(str(seq) for _, seq in shadowed)
                print(
                    f"WARNING: {chain_id}{residue_id} (author numbering) names label "
                    f"{label_seq} and nothing else. Chain {chain['label_asym_id']} also "
                    f"carries {names} at label {labels}, which author numbering cannot "
                    "address, so those residues are not conditioned here. Set "
                    "residue_numbering=label and name the label numbers to include them.",
                    flush=True,
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
    """One {select: ...} per hotspot residue, naming atoms rather than the whole residue.

    C_HOT is condition_hotspot_1_atom — an atom-level annotation, and upstream's examples
    are all atom-level: sel('A/*/56/[CG,OH]') with the comment "the PD-L1 side chain atoms
    to engage". Leaving the atom field off selects every atom of the residue, backbone
    included, which marks 19 atoms on a tyrosine where upstream would mark one or two.

    Takes [(chain, label_seq, [atom_name, ...]), ...].
    """
    return {
        "C_HOT": {
            True: select_entries(
                [
                    sel(chain, "*", label_seq, f"[{','.join(atom_names)}]")
                    for chain, label_seq, atom_names in hotspots
                ]
            )
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


PRESERVED_SEGMENT_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]*?)(\d+)(?:-(\d+))?$")


def contig_chain_ids(contig):
    """Chain ids of the preserved segments in a contig, or None if it cannot be read.

    A preserved segment names a chain and a residue range, "A1-115"; an expandable one is
    a bare length, "50-120"; "/0" is a chain break. Returning None rather than an empty
    set for anything unrecognised matters: the caller refuses selections on the strength
    of this, and a contig it cannot parse is not evidence that a chain is absent.
    """
    chains = set()
    for raw in contig.split(","):
        segment = raw.strip()
        if not segment or segment.startswith("/"):
            continue
        if LENGTH_RE.match(segment):  # expandable: a length or a range of lengths
            continue
        match = PRESERVED_SEGMENT_RE.match(segment)
        if not match:
            return None
        chains.add(match.group(1))
    return chains


def reject_offtarget(residues, target, param_name, contig=""):
    """Refuse residues naming a chain the contig does not preserve.

    A selection naming an absent chain is resolved against a system that does not contain
    it, so upstream's parser accepts the query, matches nothing, and conditions nothing —
    with no error anywhere. The run then spends full GPU time on an unaimed design while
    every log line says it was aimed.

    Naming the input's second chain is the likely mistake, not an exotic one: on a complex
    the binder usually is chain B, redesign_chain defaults to "B", and "aim at the binder's
    residues" is a reasonable thing to try. It is also the most confusing failure, because
    the designed segment is itself labelled with a free chain id, so sel('B/...') can look
    like it addresses the new binder when it addresses nothing at all.

    Which chains count as present depends on the contig actually used. The derived contigs
    are "<length>, /0, <target segment>", where the target is the only input chain. A
    contig override can preserve others, and this must read it rather than assume: refusing
    a selection with "only the target chain is preserved" when the user's own contig
    plainly preserves the chain they named is a worse failure than the one being prevented,
    because the message is false and the fix it asks for is wrong.
    """
    target_id = target["label_asym_id"]
    preserved = contig_chain_ids(contig) if contig else {target_id}
    if preserved is None:
        # An unreadable contig is the user's to debug: check_contig() and then upstream's
        # own parser will speak to it. Guessing here could only refuse something valid.
        return
    stray = sorted({chain for chain, _ in residues if chain not in preserved})
    if not stray:
        return
    if preserved == {target_id}:
        detail = (
            f"but only the target chain ({target_id}) is preserved by the contig.\n"
            f"Hotspots name residues on the TARGET that the new chain should engage, so "
            f"give them on chain {target_id}. They cannot name the binder: the binder is "
            "what is being designed and does not exist yet."
        )
    else:
        detail = (
            f"but the contig in use preserves only {', '.join(sorted(preserved))}.\n"
            "Name residues on a chain the contig carries, or widen the contig to carry "
            "the one you want."
        )
    fail(
        f"{param_name} names chain(s) {', '.join(stray)}, {detail}\n"
        "Residues on a chain the contig does not carry are absent from the designed "
        "system, so they would condition nothing and the run would produce an unaimed "
        "design."
    )


def resolve_hotspots(params, inventory, target, atoms):
    """Turn user-named hotspot residues into atom-level selections."""
    if not params["hotspot_residues"]:
        return {}
    residues = to_label_ids(
        parse_residue_tokens(params["hotspot_residues"], target["label_asym_id"]),
        inventory,
        params["residue_numbering"],
    )
    reject_offtarget(residues, target, "hotspot_residues", params["contig"])

    structure = atoms()
    hotspots = []
    for chain, label_seq in sorted(set(residues)):
        bucket = structure.get((chain, label_seq))
        if not bucket:
            fail(
                f"Residue {chain}{label_seq} (label numbering) is in the inventory but has "
                "no heavy atoms in the structure, so there is nothing to mark as a hotspot."
            )
        hotspots.append((chain, label_seq, outermost_sidechain_atoms(bucket)))

    marked = sum(len(names) for _, _, names in hotspots)
    print(
        f"Marking {marked} atom(s) across {len(hotspots)} hotspot residue(s), using the "
        "side-chain atoms reaching furthest from CA.",
        flush=True,
    )
    if marked > params["max_hotspot_atoms"]:
        print(
            f"WARNING: {marked} hotspot atoms is above max_hotspot_atoms="
            f"{params['max_hotspot_atoms']}. Upstream trained with about 5 and never more "
            "than 10% of the interface, so conditioning this hard is outside the regime "
            "the model saw. Naming fewer residues aims it more like upstream does.",
            flush=True,
        )
    return hotspot_conditions(hotspots)


def build_binder(params, inventory, atoms):
    target = chain_by_id(inventory, params["target_chain"])
    contig = f"{params['binder_length']}, /0, {chain_segment(target)}"
    conditions = hold_target(target["label_asym_id"])
    conditions.update(resolve_hotspots(params, inventory, target, atoms))
    return contig, conditions, False


# Upstream's trained definition of a hotspot, from docs/reference/annotations.md: "a
# target-chain atom within 4.5 A of any binder-chain atom", count sampled geometric(p=0.2)
# and capped at 10% of eligible atoms. Three numbers matter and all three are atom-level:
# the distance, the atom (not the residue), and the count.
#
# This is NOT node 04's 8.0 A CA-CA contact count, and the two are deliberately different.
# See the note on CONTACT_CUTOFF_A in node 04: that one describes a generated complex
# coarsely and cheaply; this one reproduces the conditioning convention the model was
# trained with. Using 8 A here would be wrong twice over -- 8 A is the Cb-Cb convention,
# not the CA-CA one -- and on 5O45 it misses 8 of the 17 real heavy-atom contacts.
DEFAULT_EPITOPE_CONTACT_CUTOFF_A = 4.5
HOTSPOT_FRACTION_OF_POOL = 0.1
DEFAULT_MAX_HOTSPOT_ATOMS = 6
MIN_HOTSPOT_SEPARATION_A = 5.0
BACKBONE = {"N", "CA", "C", "O", "OXT"}


def load_heavy_atoms(cif_path):
    """Heavy atoms grouped by (label_asym_id, label_seq_id).

    Hydrogens are excluded because upstream excludes them when it marks hotspots, and
    alternate locations are collapsed to one so the same atom is not counted twice.

    The collapse keeps the lowest-labelled altloc present for each atom rather than
    keeping altloc A specifically. Those are the same thing on most files and not on
    some: a residue modelled only at B and C -- which happens where A was rejected in
    refinement -- carries no atom labelled A at all, and whitelisting A dropped that
    residue out of the interface entirely. It could then never be picked as a hotspot
    and nothing would say why. Taking the lowest label present also makes the choice
    independent of the order the atoms appear in the file.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    cif = MMCIF2Dict(cif_path)
    asym = cif["_atom_site.label_asym_id"]
    seq = cif["_atom_site.label_seq_id"]
    name = cif["_atom_site.label_atom_id"]
    element = cif.get("_atom_site.type_symbol", ["C"] * len(asym))
    altloc = cif.get("_atom_site.label_alt_id", ["."] * len(asym))
    xs, ys, zs = (
        cif["_atom_site.Cartn_x"],
        cif["_atom_site.Cartn_y"],
        cif["_atom_site.Cartn_z"],
    )

    atoms = {}
    kept_altloc = {}
    for i in range(len(asym)):
        if seq[i] in (".", "?", ""):
            continue
        if element[i].strip().upper() in ("H", "D"):
            continue
        try:
            point = (float(xs[i]), float(ys[i]), float(zs[i]))
        except ValueError:
            continue
        key = (asym[i], int(seq[i]))
        atom_name = name[i]
        alt = "" if altloc[i] in (".", "?") else altloc[i]
        previous = kept_altloc.get((key, atom_name))
        if previous is not None and previous <= alt:
            continue
        atoms.setdefault(key, {})[atom_name] = point
        kept_altloc[(key, atom_name)] = alt
    return {key: sorted(bucket.items()) for key, bucket in atoms.items()}


def _sq_distance(a, b):
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2


def interface_atom_pool(atoms, target_chain, binder_chain, cutoff):
    """Target heavy atoms within the cutoff of any binder heavy atom.

    Returns [(distance, label_seq, atom_name, point), ...] sorted nearest first, which is
    upstream's eligible set: the atoms a hotspot may be drawn from.
    """
    binder = [
        point
        for (chain, _), bucket in atoms.items()
        if chain == binder_chain
        for _, point in bucket
    ]
    if not binder:
        return []
    cutoff_sq = cutoff ** 2
    pool = []
    for (chain, label_seq), bucket in atoms.items():
        if chain != target_chain:
            continue
        for atom_name, point in bucket:
            nearest_sq = min(_sq_distance(point, other) for other in binder)
            if nearest_sq <= cutoff_sq:
                pool.append((nearest_sq ** 0.5, label_seq, atom_name, point))
    pool.sort(key=lambda row: (row[0], row[1], row[2]))
    return pool


def choose_hotspot_atoms(pool, budget):
    """Thin the eligible pool down to the number upstream would have marked.

    Upstream samples the count from a geometric distribution (mean about 5) and caps it at
    10% of the eligible atoms. Marking the whole pool instead puts the conditioning an
    order of magnitude outside the regime the model was trained in: on 5O45 the pool is 69
    atoms and upstream's cap is 6.

    Two knowing departures from upstream, both in how the atoms are picked rather than how
    many. The pick is deterministic where upstream samples at random, because a workflow
    with a seed parameter should give the same spec for the same inputs. And it is spread
    out -- nearest contact first, one atom per residue, no two picks within
    MIN_HOTSPOT_SEPARATION_A -- so the marks cover the patch rather than piling into
    whichever side chain happens to be closest. Upstream has no such rule: it samples
    uniformly over the eligible atoms with no spacing or per-residue constraint, so the
    spread is this node's own and not a reproduction of upstream's.
    """
    if not pool:
        return []
    allowed = max(1, int(len(pool) * HOTSPOT_FRACTION_OF_POOL))
    allowed = min(allowed, budget)
    chosen = []
    for distance, label_seq, atom_name, point in pool:
        if len(chosen) >= allowed:
            break
        if any(seen_seq == label_seq for _, seen_seq, _, _ in chosen):
            continue
        if any(
            _sq_distance(point, seen_point) < MIN_HOTSPOT_SEPARATION_A ** 2
            for _, _, _, seen_point in chosen
        ):
            continue
        chosen.append((distance, label_seq, atom_name, point))
    if len(chosen) < allowed:
        # Under-delivering is legitimate -- a small or tightly packed epitope genuinely
        # has nowhere to put another well-separated mark -- but it is not obvious from
        # the outside, and the printed cap would otherwise be the last number the user
        # sees before a smaller one appears with no explanation.
        print(
            f"NOTE: {allowed} hotspot atom(s) were allowed here but only {len(chosen)} "
            f"could be placed. Picks are one per residue and at least "
            f"{MIN_HOTSPOT_SEPARATION_A:g} A apart, and this epitope has no room for "
            "more; the design is conditioned on the smaller number.",
            flush=True,
        )
    return chosen


def outermost_sidechain_atoms(bucket, count=2):
    """The side-chain heavy atoms reaching furthest from CA, or CA itself if there are none.

    Used when the user names hotspot residues directly: there is no binder yet to measure
    against, so "which atoms of this residue would a binder engage" has to be answered from
    the residue alone. The atoms at the end of the side chain are the ones that reach out
    of the surface, which is why upstream's own examples pick things like [CG,OH] on a
    tyrosine rather than the whole residue. Glycine, and anything with no side chain
    resolved, falls back to CA.
    """
    by_name = dict(bucket)
    ca = by_name.get("CA")
    sidechain = [(name, point) for name, point in bucket if name not in BACKBONE]
    if not sidechain or ca is None:
        return ["CA"] if ca is not None else [name for name, _ in bucket[:1]]
    sidechain.sort(key=lambda item: _sq_distance(item[1], ca), reverse=True)
    return [name for name, _ in sidechain[:count]]


def build_binder_redesign(params, inventory, atoms):
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
    if params["hotspot_residues"]:
        conditions.update(resolve_hotspots(params, inventory, target, atoms))
    else:
        # What distinguishes a redesign from a fresh binder: the existing binder's own
        # epitope becomes the hotspot set, so the new chain is aimed at the same patch.
        # Without this the task would be byte-identical to task=binder.
        cutoff = params["epitope_contact_cutoff_a"]
        pool = interface_atom_pool(
            atoms(),
            target["label_asym_id"],
            redesign["label_asym_id"],
            cutoff,
        )
        picks = choose_hotspot_atoms(pool, params["max_hotspot_atoms"])
        if picks:
            cap = max(1, int(len(pool) * HOTSPOT_FRACTION_OF_POOL))
            detail = ", ".join(
                f"{label_seq}/{atom_name}" for _, label_seq, atom_name, _ in picks
            )
            print(
                f"Derived {len(picks)} hotspot atom(s) from chain "
                f"{redesign['label_asym_id']}'s contacts with "
                f"{target['label_asym_id']} (label {detail}), drawn from "
                f"{len(pool)} interface atom(s) within {cutoff} A heavy-atom "
                "(epitope_contact_cutoff_a).",
                flush=True,
            )
            print(
                f"Upstream marks about 5 and never more than "
                f"{int(HOTSPOT_FRACTION_OF_POOL * 100)}% of the interface, which is "
                f"{cap} atom(s) here; max_hotspot_atoms={params['max_hotspot_atoms']}. "
                "Set hotspot_residues to name them yourself.",
                flush=True,
            )
            conditions.update(
                hotspot_conditions(
                    [
                        (target["label_asym_id"], label_seq, [atom_name])
                        for _, label_seq, atom_name, _ in picks
                    ]
                )
            )
        else:
            print(
                f"WARNING: no heavy atom of chain {target['label_asym_id']} lies within "
                f"{cutoff} A of chain {redesign['label_asym_id']}, so no epitope could be "
                "derived. This run is equivalent to task=binder. Give hotspot_residues to "
                "aim it, or raise epitope_contact_cutoff_a if the chains really are in "
                "contact.",
                flush=True,
            )
    return contig, conditions, False


def build_motif_scaffold(params, inventory, atoms):
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
    chains_used = {chain for chain, _, _ in runs}
    if len(chains_used) > 1:
        print(
            f"WARNING: motif_residues spans {len(chains_used)} chains "
            f"({', '.join(sorted(chains_used))}). They are grafted into a SINGLE designed "
            "chain with no chain break between them, which is unlikely to be what you "
            "want. Scaffold one chain's motif at a time, or supply your own contig with "
            "explicit /0 breaks.",
            flush=True,
        )
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


def build_monomer(params, inventory, atoms):
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


def report_dropped_chains(contig, inventory):
    """Name the chains the contig does not carry into the designed system.

    Everything the contig does not preserve is absent from what RFD4 is given. That is
    the right default -- the point of a contig is to choose -- but it is a different
    structure from the one that was uploaded, and nothing said so. 1BRS holds three
    copies of the barnase-barstar complex; five of its six chains vanish here without
    a word, and every interface number downstream is computed against what is left.
    """
    preserved = contig_chain_ids(contig)
    if preserved is None:
        # An unreadable contig is reported by check_contig and upstream; guessing which
        # chains it drops would be worse than saying nothing.
        return
    dropped = [c for c in inventory["chains"] if c["label_asym_id"] not in preserved]
    if not dropped:
        return
    named = ", ".join(
        f"{c['label_asym_id']} (auth {c['auth_asym_id']}, {c['n_residues']} residues)"
        for c in dropped
    )
    print(
        f"WARNING: the contig preserves {', '.join(sorted(preserved))}. "
        f"{len(dropped)} other chain(s) in target.cif are not carried into the designed "
        f"system: {named}. They are absent from the model's view and from every interface "
        "number computed downstream. Widen the contig if one of them is part of the "
        "interface you are designing against.",
        flush=True,
    )


def check_contig(contig):
    """Upstream rejects a contig whose segments are all expandable.

    This also enforces designed-segment-first, which is not upstream's rule but ours.
    Node 04 identifies the design as the first chain of each generated structure, and
    RFD4 writes chains in contig order, so a contig that names the target first makes
    every per-design metric describe the target instead. Nothing downstream can notice:
    the chain count is unchanged, the numbers are all plausible, and node 04's
    several-chains warning does not fire on the ordinary two-chain case.

    The contigs this node derives always put the length first, so this only ever bites a
    hand-written `contig` override.
    """
    segments = [s.strip() for s in contig.split(",") if s.strip() and s.strip() != "/0"]
    preserved = [s for s in segments if not LENGTH_RE.match(s)]
    if not preserved:
        fail(
            f"Contig '{contig}' has no preserved segment. Upstream rejects a contig built "
            "only from expandable lengths, since the input exists to be referenced."
        )
    if segments and not LENGTH_RE.match(segments[0]):
        designed = next((s for s in segments if LENGTH_RE.match(s)), None)
        rewritten = (
            f"{designed}, /0, " + ", ".join(s for s in segments if s != designed)
            if designed
            else None
        )
        fail(
            f"Contig '{contig}' puts the preserved segment '{segments[0]}' first. The "
            "designed length has to come first: node 04 reads the first chain of each "
            "generated structure as the design, so a target-first contig makes "
            "designed_chain, designed_length, the designed sequence, the hydrophobic "
            "fraction, the net charge and the radius of gyration all describe the target "
            "-- with no error and no warning, because the chain count is unchanged.\n"
            + (f"Write it as '{rewritten}'." if rewritten else
               "Put the designed length first.")
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


HOTSPOT_COLUMN = "_atom_site.condition_hotspot_1_atom"


def count_flagged(cif_path, column):
    """Atoms the rendered spec actually marks in a boolean condition column.

    Returns None when the column is absent or the file cannot be read, so a parsing
    problem here never turns into a spurious failure of an otherwise valid spec.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    try:
        cif = MMCIF2Dict(cif_path)
    except Exception:
        return None
    values = cif.get(column)
    if values is None:
        return None
    if isinstance(values, str):
        # Measured against the pinned Biopython: every category comes back as a list,
        # a one-row loop and a bare key/value pair included, so this does not fire
        # today. It stays because the failure if it ever did would be silent and
        # wrong-way-round -- a string is iterable, so the count would run over its
        # characters, reach 0, and report a valid spec as marking no hotspot at all.
        values = [values]
    return sum(1 for value in values if str(value).strip() == "True")


def hotspot_selections(conditions):
    return [
        entry["select"]
        for entries in conditions["C_HOT"].values()
        for entry in entries
    ]


def expected_hotspot_atoms(conditions):
    """How many atoms the C_HOT selections ask for, or None if that cannot be counted.

    The selections are the ones sel() just produced, so the atom list is always the
    bracketed final field: sel('A/*/38/[CD1,CG1]'). Anything that does not match that
    shape returns None and the caller falls back to the total-miss check rather than
    failing a spec on a counting rule it cannot apply.
    """
    total = 0
    for query in hotspot_selections(conditions):
        match = re.search(r"/\[([^\]]*)\]'\)$", query)
        if not match:
            return None
        names = [name for name in match.group(1).split(",") if name.strip()]
        if not names:
            return None
        total += len(names)
    return total


def verify_hotspots_landed(preview_path, conditions):
    """Check the rendered system really carries the hotspots the spec asked for.

    reject_offtarget() catches the known way to select nothing. This catches the rest:
    it reads what upstream actually produced rather than trusting that a query which
    parsed must also have matched. A selection that matches no atom is not an error to
    upstream — an empty match is simply an empty set — so this node has to be the one
    that objects, which is the whole reason it renders the spec on CPU first.
    """
    if "C_HOT" not in conditions:
        return
    expected = expected_hotspot_atoms(conditions)
    flagged = count_flagged(preview_path, HOTSPOT_COLUMN)
    if flagged is None:
        print(
            "WARNING: could not read the hotspot annotation back from the rendered spec, "
            "so it was not confirmed that the hotspots match any atom.",
            flush=True,
        )
        return
    queries = " ".join(hotspot_selections(conditions))
    if flagged == 0:
        fail(
            f"The spec asks for hotspots ({queries}) but the rendered system flags no "
            "atoms at all, so C_HOT conditions nothing and generation would produce an "
            "unaimed design. The selection parsed but matched nothing — check the chain "
            "and residue ids against the inventory printed by node 01."
        )
    # A partial miss is the more dangerous shape: most of the aim lands, one selection
    # quietly matches nothing, and the design is aimed at a patch the user did not ask
    # for. Comparing against the atoms actually requested is what catches it; a bare
    # non-zero check passes it at rc=0.
    if expected is not None and flagged < expected:
        fail(
            f"The spec asks for {expected} hotspot atom(s) ({queries}) but the rendered "
            f"system flags only {flagged}. At least one selection matched nothing, so the "
            "design would be aimed at only part of the patch you asked for. Check every "
            "chain, residue id and atom name against the inventory printed by node 01, "
            "and confirm each residue is inside the contig."
        )
    # Over-flagging is not a mis-aim and is not worth refusing a run over: a residue kept
    # at more than one altloc can carry the same atom name twice, and both copies are
    # correctly inside the patch the user asked for. Say so rather than fail.
    if expected is not None and flagged > expected:
        print(
            f"NOTE: the spec names {expected} hotspot atom(s) and the rendered system "
            f"flags {flagged}. Every requested atom landed; the extras are additional "
            "copies of the same atom names (alternate conformations).",
            flush=True,
        )
    print(f"C_HOT marks {flagged} atom(s) in the rendered system.", flush=True)


def write_annotation(target_cif, conditions, dest):
    """Fallback annotation: the target plus a record of what each condition selects."""
    shutil.copyfile(target_cif, dest)
    if os.path.getsize(dest) == 0:
        fail(f"'{target_cif}' is empty; node 01 produced no usable structure.")
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

    # epitope_cutoff_a was a CA-CA distance defaulting to 8.0. Its replacement is a
    # heavy-atom distance defaulting to 4.5, so the same number means something different
    # and a value carried over would silently widen the epitope. Refuse it rather than
    # reinterpret it.
    if os.environ.get("PARAM_EPITOPE_CUTOFF_A"):
        fail(
            "epitope_cutoff_a has been replaced by epitope_contact_cutoff_a, which is a "
            "heavy-atom distance rather than a CA-CA one and defaults to 4.5 A, matching "
            "the criterion RFD4 was trained with. The same number does not mean the same "
            "thing in both, so the old parameter is not reused. Remove it, and set "
            "epitope_contact_cutoff_a only if you need something other than 4.5."
        )

    params["epitope_contact_cutoff_a"] = numeric_param(
        "PARAM_EPITOPE_CONTACT_CUTOFF_A",
        "epitope_contact_cutoff_a",
        DEFAULT_EPITOPE_CONTACT_CUTOFF_A,
        float,
    )
    params["max_hotspot_atoms"] = numeric_param(
        "PARAM_MAX_HOTSPOT_ATOMS", "max_hotspot_atoms", DEFAULT_MAX_HOTSPOT_ATOMS, int
    )

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
    cutoff = params["epitope_contact_cutoff_a"]
    if cutoff <= 0:
        fail(f"epitope_contact_cutoff_a must be positive, got {cutoff}.")
    if not 3.0 <= cutoff <= 8.0:
        print(
            f"WARNING: epitope_contact_cutoff_a={cutoff} A is outside the range that reads "
            "as heavy-atom contact. Below about 3 A almost nothing qualifies, since that is "
            "shorter than a hydrogen bond; much above 5 A it stops being contact, and by "
            "8 A it is the whole face of the target. Upstream trains at 4.5 A.",
            flush=True,
        )
    if params["max_hotspot_atoms"] < 1:
        fail(
            f"max_hotspot_atoms must be at least 1, got {params['max_hotspot_atoms']}."
        )
    if params["contig"] and params["task"] == "monomer":
        fail(
            "contig cannot be combined with task=monomer: a monomer is designed from a "
            "generated scaffold with no contig at all. Use task=binder (or another "
            "conditioned task) if you want to supply a contig."
        )
    low_s, _, high_s = params["binder_length"].partition("-")
    low_n, high_n = int(low_s), int(high_s or low_s)
    if low_n < 1:
        fail(f"binder_length '{params['binder_length']}' must be at least 1 residue.")
    if low_n > high_n:
        fail(
            f"binder_length '{params['binder_length']}' runs backwards; write it as "
            f"{high_n}-{low_n}."
        )

    target_cif = "./inputs/target.cif"
    inventory_path = "./inputs/target_inventory.json"
    for path in (target_cif, inventory_path):
        if not os.path.exists(path):
            fail(f"Required input not found: {path}")
    inventory = load_inventory(inventory_path)

    consumes = {
        "binder": {"hotspot_residues"},
        "binder_redesign": {"hotspot_residues"},
        "motif_scaffold": {"motif_residues"},
        "monomer": set(),
    }[params["task"]]
    for name in ("hotspot_residues", "motif_residues"):
        if params[name] and name not in consumes:
            print(
                f"WARNING: {name}='{params[name]}' is ignored by task={params['task']} "
                "and conditions nothing. Check you meant the other parameter.",
                flush=True,
            )

    # Loaded at most once, and only by the tasks that need coordinates: monomer ignores
    # the target entirely and motif_scaffold works from the inventory.
    cache = {}

    def atoms():
        if "value" not in cache:
            cache["value"] = load_heavy_atoms(target_cif)
        return cache["value"]

    contig, conditions, as_annotated_cif = BUILDERS[params["task"]](
        params, inventory, atoms
    )

    if params["contig"]:
        print(
            f"NOTE: contig override in use; derived contig '{contig}' discarded.",
            flush=True,
        )
        contig = params["contig"]
        if not conditions:
            # Every templated task holds its target; an override must not silently drop
            # that, or the model is handed a "target" it is free to move and rewrite.
            fail(
                "contig override produced a spec with no conditions, so the target would "
                "not be held. This is a bug in the task builder, not in your input."
            )
    if contig:
        check_contig(contig)
        report_dropped_chains(contig, inventory)

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
    if validate_with_upstream(spec_path, preview_path):
        verify_hotspots_landed(preview_path, conditions)

    print(f"Wrote {spec_path} and {preview_path}", flush=True)


if __name__ == "__main__":
    main()
