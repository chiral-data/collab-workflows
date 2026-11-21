#!/usr/bin/env python3
"""Check a GPU campaign's output against a measurement of its input.

Exit codes say a run finished, not that it was right. Every check here compares something
a node *reported* against something measured independently: `tests/probe_entry.py` for the
input, a scaling law or a recomputation for the output. The detail string carries the
numbers, so a failure reads as evidence rather than as an assertion.

Run:  python3 -I tests/check_campaign.py <results-dir> <downloads-dir>

Both directories are produced by tests/run_gpu_campaign.sh.
"""
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from probe_entry import probe  # noqa: E402

PASS, FAIL, GAP, SKIP = "PASS", "FAIL", "GAP", "SKIP"

# Below this the "ensemble" is one structure plus sampling noise. It is the scale of
# ordinary backbone fluctuation, and it is a floor for spotting collapse -- not a
# validation of the ensemble. "Non-zero spread" is not a test: two independently sampled
# frames always differ.
COLLAPSE_FLOOR_NM = 0.1

# Rg of a compact globular protein, nm, from the residue count. The same relation
# workflow-037 uses in angstrom (2.2 N^0.38). A factor of 1.5 either way is generous
# enough for a real fold and still catches an unfolded or exploded structure.
#
# Generous on purpose: the relation is for a *globular* fold, and a fold that is elongated by
# design sits well above it without being wrong. Measured here: 1UBQ 1.00x, 5PTI 1.08x, and
# 2ZTA — the GCN4 leucine zipper, a parallel coiled coil — 1.39x. Tightening this limit would
# flag the coiled coil, so it stays loose and the check keeps its real job: catching a chain
# that has unfolded or exploded, which lands far outside any fold's range.
RG_TOLERANCE = 1.5


def rg_expected_nm(n_residues):
    return 0.22 * (n_residues ** 0.38)


class Case:
    """Everything one system's run produced, loaded once."""

    def __init__(self, results, downloads, entry):
        self.id = entry
        base = os.path.join(results, entry)
        # A case named "<ID>-propagated" is that entry run with propagation on. The DAG hands
        # node 04 only what node 03 wrote, so the propagated case has no node 01 or node 02
        # directory of its own -- but it is the same system, so those facts are the source
        # case's. Reading them from there is the difference between checking the propagated
        # run against its input and reporting four blanks.
        self.source = entry[: -len("-propagated")] if entry.endswith("-propagated") else entry
        upstream = os.path.join(results, self.source)

        def _either(node, name):
            """The file from this case if it ran that node, else from the source case."""
            own = os.path.join(base, node, "outputs", name)
            return own if os.path.exists(own) else os.path.join(upstream, node, "outputs", name)

        self.inventory = _json(_either("01", "system_inventory.json"))
        self.msa_summary = _json(_either("02", "msa_summary.json"))
        self.sampling = _json(os.path.join(base, "03", "outputs", "sampling_manifest.json"))
        self.physicality = _json(os.path.join(base, "03", "outputs", "physicality.json"))
        self.propagation = _json(os.path.join(base, "04", "outputs", "propagation_manifest.json"))
        self.propagated_physicality = _json(
            os.path.join(base, "04", "outputs", "propagated_physicality.json")
        )
        self.summary = _json(os.path.join(base, "05", "outputs", "ensemble_summary.json"))
        self.report = _text(os.path.join(base, "06", "outputs", "report.html"))
        self.top_pdb = _text(os.path.join(base, "05", "outputs", "ensemble_top.pdb"))
        self.logs = {
            name: _text(os.path.join(base, f"{name}.log"))
            for name in ("01", "02", "03", "04", "05", "06")
        }
        # The exit status each node actually returned, which the driver records next to the
        # log. A node that was never run has no file and no status; that is not the same thing
        # as a node that ran and returned 0, and the checks below need to tell them apart.
        self.rc = {}
        for name in ("01", "02", "03", "04", "05", "06"):
            text = _text(os.path.join(base, f"{name}.rc"))
            if text and text.strip().isdigit():
                self.rc[name] = int(text.strip())

        src = os.path.join(downloads, self.source + ".cif")
        self.truth = probe(src) if os.path.exists(src) else {}

    @property
    def sampled(self):
        return ((self.summary or {}).get("ensembles") or {}).get("sampled") or {}

    @property
    def is_protein(self):
        return bool(self.truth.get("protein_chains"))

    @property
    def modelled_protein_residues(self):
        """(standard, total) protein residues in the system BioEmu2 was given.

        From node 01's inventory, not from the file: `System.from_cif` builds the entity
        sequence, which upstream documents as including unresolved residues, so 1CLL is 148
        modelled against 144 resolved. The molecule in the frames is the modelled one, so the
        C-alpha count and the Rg scaling both have to use these numbers.

        Two of them, because `X` in an entity sequence is a non-standard component and whether
        it carries a C-alpha depends on what it is. 2ZTA's `X` is an N-terminal acetyl cap,
        which has none, so its 34-residue chains give 33 C-alphas each; a modified amino acid
        such as selenomethionine would give one. The C-alpha count therefore has to fall
        between the standard count and the total, and the interesting check is not which end
        it lands on but whether an ion has been counted among them.
        """
        standard = total = 0
        for chain in ((self.inventory or {}).get("chains") or []):
            if "protein" not in str(chain.get("entity_type") or "").lower():
                continue
            sequence = chain.get("sequence")
            length = int(chain.get("length") or 0)
            total += length
            standard += length - (sequence.count("X") if isinstance(sequence, str) else 0)
        return (standard, total) if total else (None, None)


def _json(path, default=None):
    if not os.path.exists(path):
        return default
    try:
        with open(path) as handle:
            return json.load(handle)
    except Exception:
        return default


def _text(path):
    if not os.path.exists(path):
        return ""
    with open(path, errors="replace") as handle:
        return handle.read()


# ---------------------------------------------------------------- checks

def check_token_count(case):
    """Node 01's token total must be a measured number, not null.

    This read the wrong attribute (`n_tokens` for `num_tokens`) and so was null on every
    input ever run, which also made the >10k GPU-fit warning unreachable. One token per
    polymer residue plus one per ligand atom, so the total is at least the residue count.
    """
    if not case.inventory:
        return FAIL, "node 01 wrote no inventory"
    tokens = case.inventory.get("n_tokens")
    source = case.inventory.get("n_tokens_from")
    residues = case.truth.get("n_polymer_residues")
    if not isinstance(tokens, int):
        return FAIL, f"n_tokens={tokens!r}, so the GPU-fit warning cannot fire"
    if not source:
        return FAIL, f"n_tokens={tokens} but n_tokens_from is empty: provenance unrecorded"
    if residues and tokens < residues:
        return FAIL, f"n_tokens={tokens} < {residues} polymer residues measured in the file"
    return PASS, f"n_tokens={tokens} from {source}, against {residues} polymer residue(s)"


def check_chains_match_probe(case):
    """Node 01's chain list must be a subset of what the file holds, and name the proteins."""
    if not case.inventory:
        return FAIL, "node 01 wrote no inventory"
    reported = set(case.inventory.get("chain_ids") or [])
    measured = set(case.truth.get("chain_ids") or [])
    if not measured:
        return SKIP, "input was not probed"
    extra = reported - measured
    if extra:
        return FAIL, f"inventory names chains {sorted(extra)} that the file does not contain"
    proteins = set(case.inventory.get("protein_chains") or [])
    measured_proteins = set(case.truth.get("protein_chains") or [])
    if proteins - measured:
        return FAIL, f"protein_chains {sorted(proteins - measured)} are not in the file"
    return PASS, (
        f"{len(reported)} chain(s) reported of {len(measured)} in the file; "
        f"protein {sorted(proteins)} against measured {sorted(measured_proteins)}"
    )


def check_paired_alignment(case):
    """A paired alignment is due exactly when the system has >1 unique protein sequence.

    ColabFold returns a positionally paired alignment only for a heteromer. Node 02 wrote
    only the unpaired half, and `load_a3m_msa` always returns `MSA(paired_a3m=None)`, so
    the paired feature the model expects was replaced with zeros on every complex. The
    counter-case matters as much: a homodimer has two chains and one sequence, and must
    NOT be reported as paired.
    """
    unique = case.truth.get("unique_protein_sequences")
    if unique is None:
        return SKIP, "input was not probed"
    if not case.is_protein:
        return SKIP, "no protein chain, so no alignment of either kind is due"
    due = unique > 1
    node02 = (case.msa_summary or {}).get("paired")
    if node02 is None:
        return FAIL, "node 02 recorded no `paired` field"
    if bool(node02) != due:
        return FAIL, (
            f"{unique} unique protein sequence(s), so paired={due} is due, but node 02 "
            f"recorded paired={node02}"
        )
    if not case.sampling:
        # Node 03 did not get as far as writing a manifest, so what reached the model cannot
        # be judged. Node 02's half is still a result worth reporting.
        return FAIL, (
            f"{unique} unique sequence(s) and node 02 recorded paired={node02} correctly, but "
            "node 03 wrote no manifest, so what reached the model is unknown"
        )
    node03 = case.sampling.get("msa_paired")
    if bool(node03) != due:
        return FAIL, (
            f"node 02 recorded paired={node02} but node 03 used msa_paired={node03}: the "
            "paired alignment was written and then dropped on the way to the model"
        )
    return PASS, (
        f"{unique} unique sequence(s) -> paired={due}; node 02 and node 03 agree "
        f"({case.sampling.get('msa_how')})"
    )


def check_msa_reached_the_model(case):
    """The alignment the model used, not the one node 02 fetched. These disagreed silently."""
    if not case.sampling:
        return FAIL, "node 03 wrote no manifest"
    loaded = case.sampling.get("msa_loaded")
    how = case.sampling.get("msa_how")
    if case.is_protein:
        if not loaded:
            return FAIL, f"protein system sampled with no alignment: {how}"
        files = case.sampling.get("msa_files") or []
        if not files:
            return FAIL, "msa_loaded is true but no files were recorded"
        return PASS, f"{len(files)} file(s): {how}"
    if loaded:
        return FAIL, f"no protein chain, yet an alignment was loaded: {how}"
    return PASS, f"no protein chain, no alignment: {how}"


def check_ca_selection_excludes_ions(case):
    """The C-alpha set must be the protein's C-alphas and nothing else.

    A calcium ion is an atom named CA in a residue named CA. mdtraj matches atom names, so
    `name CA` folds structural calcium into the C-alpha set and shifts every RMSD, every
    RMSF and the medoid frame. The per-residue RMSF list is one entry per selected atom,
    so its length IS the size of the selection -- the one number that makes this visible.
    """
    if not case.sampled:
        return FAIL, "node 05 wrote no summary"
    selected = len(case.sampled.get("residue_rmsf") or [])
    standard, total = case.modelled_protein_residues
    resolved = case.truth.get("ca_atoms")
    calcium = case.truth.get("calcium_atoms")
    if not case.is_protein:
        return SKIP, "no protein chain; the all-atom fallback is the correct path here"
    if standard is None:
        return FAIL, "node 01's inventory does not state the protein chains' sequences"
    for ion_count, label in ((calcium, "calcium ion"),):
        if ion_count and selected in (standard + ion_count, total + ion_count):
            return FAIL, (
                f"{selected} atoms selected = {selected - ion_count} protein residue(s) + "
                f"{ion_count} {label}(s). The ions are being counted as C-alphas, which moves "
                "every RMSD, the RMSF and the medoid"
            )
    if not (standard <= selected <= total):
        return FAIL, (
            f"{selected} atoms selected, outside the {standard}-{total} protein residue(s) the "
            f"system models ({resolved} C-alpha resolved in the file)"
        )
    detail = f"{selected} C-alpha selected, within the {standard}-{total} modelled"
    if resolved is not None and resolved != selected:
        detail += f" ({resolved} resolved in the file)"
    if calcium:
        detail += f"; the {calcium} calcium ion(s) were excluded, as they must be"
    return PASS, detail


def check_selection_recorded(case):
    """Node 05 must say which selection it used; 'C-alpha' alone was not falsifiable."""
    selection = case.sampled.get("rmsd_atom_selection")
    if not selection:
        return FAIL, "no rmsd_atom_selection recorded"
    allowed = (
        "protein and name CA",
        "name CA in residues with more than one atom",
        "all atoms",
    )
    if selection not in allowed:
        return FAIL, f"unrecognised selection {selection!r}"
    if case.is_protein and selection == "all atoms":
        return FAIL, "a protein system fell back to all atoms; the C-alpha selection failed"
    return PASS, f"selection recorded as {selection!r}"


def check_frames_and_physicality(case):
    """One physicality record per frame, zero-based and contiguous, in frame order."""
    n_frames = case.sampled.get("n_frames")
    records = case.physicality
    if n_frames is None:
        return FAIL, "node 05 reported no frame count"
    if not isinstance(records, list):
        return FAIL, f"physicality.json is {type(records).__name__}, not a list"
    frames = [r.get("frame") for r in records if isinstance(r, dict)]
    if frames != list(range(n_frames)):
        return FAIL, (
            f"{len(frames)} record(s) for {n_frames} frame(s); frames are not a contiguous "
            f"zero-based range (first few: {frames[:5]})"
        )
    statuses = sorted({r.get("status") for r in records})
    return PASS, f"{n_frames} frame(s), one record each, status(es) {statuses}"


def check_not_collapsed(case):
    """The ensemble must show more spread than ordinary backbone noise."""
    if not case.is_protein:
        return SKIP, "collapse floor is stated for folded proteins"
    pairwise = case.sampled.get("pairwise_rmsd") or {}
    maximum = pairwise.get("max_nm")
    if maximum is None:
        return FAIL, "no pairwise RMSD recorded"
    if maximum <= COLLAPSE_FLOOR_NM:
        return FAIL, (
            f"max pairwise C-alpha RMSD {maximum} nm is at or below the {COLLAPSE_FLOOR_NM} nm "
            "floor: this is one structure plus sampling noise, not an ensemble"
        )
    return PASS, (
        f"max pairwise C-alpha RMSD {maximum} nm over {pairwise.get('n_pairs')} pair(s), "
        f"above the {COLLAPSE_FLOOR_NM} nm floor"
    )


def check_radius_of_gyration(case):
    """Rg must sit near the globular scaling 0.22 N^0.38 nm.

    This is the one check here that can catch a model that has stopped working. Everything
    else -- frame counts, physicality records, cluster counts -- is computable on an
    unfolded or exploded chain and looks ordinary. Rg is not.
    """
    if not case.is_protein:
        return SKIP, "the scaling relation is for folded proteins"
    rg = (case.sampled.get("radius_of_gyration_nm") or {}).get("mean")
    # The modelled residue count, since the frames hold the modelled molecule.
    residues = case.modelled_protein_residues[0] or case.truth.get("n_polymer_residues")
    if rg is None or not residues:
        return SKIP, "no Rg or no residue count"
    expected = rg_expected_nm(residues)
    ratio = rg / expected
    detail = (
        f"Rg {rg:.3f} nm against {expected:.3f} nm expected for {residues} residues "
        f"(0.22 N^0.38), ratio {ratio:.2f}"
    )
    if ratio > RG_TOLERANCE:
        return FAIL, detail + " - too expanded for a folded chain"
    if ratio < 1.0 / RG_TOLERANCE:
        return FAIL, detail + " - too compact to be a real fold"
    return PASS, detail


def check_medoid_is_argmin(case):
    """Recompute the medoid from the matrix the run published and compare.

    Node 05 set `medoid = 0` above 200 frames while keeping the name `medoid_frame`. Where
    the matrix is published this is checkable directly rather than by trusting the label.
    """
    matrix = (case.summary or {}).get("rmsd_matrix")
    medoid = case.sampled.get("medoid_frame")
    if not matrix:
        n_frames = case.sampled.get("n_frames")
        pairs = (case.sampled.get("pairwise_rmsd") or {}).get("n_pairs")
        if n_frames and pairs == n_frames * (n_frames - 1) // 2:
            return PASS, (
                f"matrix not published, but the spread covers all {pairs} pair(s), so the "
                f"medoid ({medoid}) was measured rather than defaulted"
            )
        return FAIL, f"no matrix and only {pairs} pair(s) measured for {n_frames} frame(s)"
    sums = [sum(row) for row in matrix]
    best = min(range(len(sums)), key=lambda i: sums[i])
    if medoid != best:
        return FAIL, (
            f"medoid_frame={medoid} but the published matrix's row sums are minimised at "
            f"{best} ({sums[best]:.4f} against {sums[medoid]:.4f})"
        )
    return PASS, f"medoid_frame={medoid} is the argmin of the published matrix's row sums"


def check_viewer_frames_superposed(case):
    """The viewer's models must be mutually superposed, measured from the file.

    `ensemble_top.pdb` was written from the un-superposed frames while the RMSF came from
    superposed ones, so the 3D viewer animated the whole molecule tumbling. Two frames that
    differ only by a rigid rotation have a large *raw* coordinate RMS difference and a small
    superposed RMSD, so comparing the raw difference against the reported superposed spread
    separates the two cases without needing a superposition here.
    """
    models = _pdb_ca_models(case.top_pdb)
    if len(models) < 2:
        return SKIP, f"{len(models)} model(s) in ensemble_top.pdb; need two to compare"
    reference = models[0]
    worst = 0.0
    for other in models[1:]:
        if len(other) != len(reference):
            return FAIL, "models have different atom counts"
        total = sum(
            (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2
            for a, b in zip(reference, other)
        )
        worst = max(worst, math.sqrt(total / len(reference)))
    worst_nm = worst / 10.0  # the PDB is in angstrom; everything else here is nm
    reported = (case.sampled.get("pairwise_rmsd") or {}).get("max_nm")
    if reported is None:
        return SKIP, "no reported spread to compare against"
    limit = max(3.0 * reported, reported + 0.05)
    detail = (
        f"worst raw model-to-model C-alpha RMS difference {worst_nm:.3f} nm against a "
        f"reported superposed spread of {reported:.3f} nm"
    )
    if worst_nm > limit:
        return FAIL, detail + " - the models are not superposed; the viewer will show tumbling"
    return PASS, detail


def check_report(case):
    """The report must exist and carry the claims the summary makes."""
    if not case.report:
        return FAIL, "node 06 wrote no report"
    required = [
        ("the selection", case.sampled.get("rmsd_atom_selection") or ""),
        ("the caveat block", "Frame counts are not populations"),
        ("the data object", "const DATA ="),
    ]
    missing = [label for label, text in required if text and text not in case.report]
    if missing:
        return FAIL, f"report is missing {', '.join(missing)}"
    return PASS, f"{round(len(case.report) / 1024)} KB, carrying the selection and the caveats"


def check_propagation_disabled_cleanly(case):
    """Propagation off must still write every declared output, and say it was an identity."""
    if not case.propagation:
        return FAIL, "node 04 wrote no manifest"
    if case.propagation.get("enabled"):
        return SKIP, "propagation ran in this case"
    if (case.summary or {}).get("propagation_enabled"):
        return FAIL, "node 04 recorded enabled=false but node 05 reported propagation enabled"
    if "propagated" in ((case.summary or {}).get("ensembles") or {}):
        return FAIL, "an identity pass-through was analysed as a second, time-lagged ensemble"
    return PASS, f"identity pass-through recorded: {case.propagation.get('reason')}"


def check_propagation_ran(case):
    """When propagation is on, it must produce a second ensemble of the same size.

    This is the only place the time-lagged path is executed with a model: the offline suite
    reaches node 04's identity pass-through and nothing else. `propagate_with_embeddings`
    requires byte-identical starting frames inside a resampling group, which is why the
    resampling samplers are refused outright; the frames it returns are the same conformations
    advanced by `lag_time_ns`, so the count must match and the molecule must not have changed.
    """
    manifest = case.propagation
    if not manifest:
        return FAIL, "node 04 wrote no manifest"
    if not manifest.get("enabled"):
        return SKIP, "propagation was disabled in this case"
    ensembles = (case.summary or {}).get("ensembles") or {}
    if "propagated" not in ensembles:
        return FAIL, f"node 04 propagated and node 05 analysed only {sorted(ensembles)}"
    sampled, propagated = ensembles.get("sampled") or {}, ensembles["propagated"]
    if propagated.get("n_frames") != sampled.get("n_frames"):
        return FAIL, (
            f"{propagated.get('n_frames')} propagated frame(s) from "
            f"{sampled.get('n_frames')} sampled"
        )
    if propagated.get("n_atoms") != sampled.get("n_atoms"):
        return FAIL, (
            f"{propagated.get('n_atoms')} atoms after propagation against "
            f"{sampled.get('n_atoms')} before: the molecule changed"
        )
    records = case.propagated_physicality
    if not isinstance(records, list) or len(records) != propagated["n_frames"]:
        return FAIL, (
            f"propagated_physicality.json holds "
            f"{len(records) if isinstance(records, list) else records!r} record(s) for "
            f"{propagated['n_frames']} frame(s)"
        )
    moved = (propagated.get("pairwise_rmsd") or {}).get("max_nm")
    return PASS, (
        f"{propagated['n_frames']} frame(s) at lag {manifest.get('lag_time_ns')} ns, "
        f"{propagated['n_atoms']} atoms unchanged, {len(records)} physicality record(s); "
        f"spread within the propagated ensemble {moved} nm"
    )


def check_propagated_is_physical(case):
    """And the propagated ensemble must still be a folded protein.

    The failure mode of a time-lagged model is not a crash: it is frames that have drifted
    somewhere unphysical while every count stays right. Rg against the globular scaling is the
    check that notices.
    """
    manifest = case.propagation
    if not manifest or not manifest.get("enabled"):
        return SKIP, "propagation was disabled in this case"
    propagated = ((case.summary or {}).get("ensembles") or {}).get("propagated") or {}
    rg = (propagated.get("radius_of_gyration_nm") or {}).get("mean")
    residues = case.modelled_protein_residues[0] or case.truth.get("n_polymer_residues")
    if rg is None or not residues:
        return SKIP, "no propagated Rg or no residue count"
    expected = rg_expected_nm(residues)
    ratio = rg / expected
    sampled_rg = ((case.sampled or {}).get("radius_of_gyration_nm") or {}).get("mean")
    detail = (
        f"propagated Rg {rg:.3f} nm against {expected:.3f} nm expected and {sampled_rg} nm "
        f"before propagation, ratio {ratio:.2f}"
    )
    if ratio > RG_TOLERANCE or ratio < 1.0 / RG_TOLERANCE:
        return FAIL, detail + " - the propagated frames are not a folded chain"
    return PASS, detail


def _pdb_ca_models(text):
    """CA coordinates per MODEL, parsed by column as the PDB format fixes them."""
    models, current = [], None
    for line in text.splitlines():
        if line.startswith("MODEL"):
            current = []
            continue
        if line.startswith("ENDMDL"):
            if current:
                models.append(current)
            current = None
            continue
        if line.startswith(("ATOM", "HETATM")) and line[12:16].strip() == "CA":
            if current is None:
                current = []
            try:
                current.append(
                    (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                )
            except ValueError:
                continue
    if current:
        models.append(current)
    return models


NODE_ORDER = ("01", "02", "03", "04", "05", "06")


def short_reason(line, limit=110):
    """The readable front of a node's ERROR line: no prefix, cut at a sentence end.

    CUDA's own message runs to several hundred characters of allocator advice, and a
    mid-word truncation of it is the kind of detail string that reads as noise rather than
    as evidence.
    """
    text = re.sub(r"^ERROR:?\s*", "", line.strip())
    if len(text) <= limit:
        return text
    cut = text.rfind(". ", 0, limit)
    return (text[: cut + 1] if cut > 40 else text[:limit].rstrip() + "...")


# A failure that means "this hardware is too small", rather than "this code is wrong". The
# distinction is what lets the campaign be a gate at all: an 8 GB card cannot hold a 542-token
# system at any batch size, and a run that stops there has behaved correctly. Everything else
# that stops a node is a defect until shown otherwise.
RESOURCE_LIMITS = ("out of memory", "outofmemory", "no cuda-capable device", "no gpu")


def first_aborted(case):
    """(node, reason) for the earliest node that exited non-zero, or (None, None).

    The reason is the node's own last ERROR line, or None if it wrote none -- stopping
    *quietly* being the thing this whole suite exists to catch.
    """
    for node in NODE_ORDER:
        if case.rc.get(node, 0) == 0:
            continue
        reported = [
            line.strip()
            for line in (case.logs.get(node) or "").splitlines()
            if line.strip().startswith("ERROR")
        ]
        return node, (reported[-1] if reported else None)
    return None, None


def check_failure_was_reported(case):
    """A node that stops must say why, and the why must be a resource limit.

    Three outcomes, and the middle one is the point. A node that stops without writing an
    ERROR line fails, because a silent stop is the defect class this suite is built around.
    A node that stops on a resource limit passes: the card is too small, the node said so,
    and the checks below it report as skipped rather than as failures. A node that stops for
    any *other* stated reason fails here even though it reported itself, so that the campaign
    exits non-zero -- otherwise a real defect that the node happens to describe politely
    would leave a green run behind a wall of skips.
    """
    if not case.rc:
        return SKIP, "the driver recorded no exit statuses for this case"
    node, reason = first_aborted(case)
    if node is None:
        return PASS, f"every node exited 0 ({len(case.rc)} recorded)"
    if reason is None:
        return FAIL, (
            f"node {node} exited {case.rc[node]} without writing an ERROR line: "
            "the run stopped and the log does not say why"
        )
    if any(limit in reason.lower() for limit in RESOURCE_LIMITS):
        return PASS, (
            f"node {node} exited {case.rc[node]} on a resource limit, reported: "
            f"{short_reason(reason)}"
        )
    return FAIL, (
        f"node {node} exited {case.rc[node]} for a reason that is not a resource limit, so "
        f"the checks below it are skipped against a defect: {short_reason(reason)}"
    )


# (label, function, the node whose output the check reads). The node matters because a run
# that stopped at node 03 cannot be judged on node 05's numbers: those checks report as
# skipped-with-a-reason rather than as failures, so a FAIL always means a node that ran.
CHECKS = [
    ("stopped only on a resource limit", check_failure_was_reported, "01"),
    ("token count measured", check_token_count, "01"),
    ("chains match the file", check_chains_match_probe, "01"),
    ("paired alignment due iff heteromer", check_paired_alignment, "03"),
    ("alignment reached the model", check_msa_reached_the_model, "03"),
    ("C-alpha set excludes ions", check_ca_selection_excludes_ions, "05"),
    ("selection recorded", check_selection_recorded, "05"),
    ("one physicality record per frame", check_frames_and_physicality, "05"),
    ("ensemble not collapsed", check_not_collapsed, "05"),
    ("Rg near the globular scaling", check_radius_of_gyration, "05"),
    ("medoid is the argmin", check_medoid_is_argmin, "05"),
    ("viewer models superposed", check_viewer_frames_superposed, "05"),
    ("report carries its claims", check_report, "06"),
    ("identity propagation is labelled", check_propagation_disabled_cleanly, "04"),
    ("propagation produced a second ensemble", check_propagation_ran, "04"),
    ("propagated frames are still folded", check_propagated_is_physical, "05"),
]


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    results, downloads = sys.argv[1], sys.argv[2]
    entries = sorted(
        name for name in os.listdir(results) if os.path.isdir(os.path.join(results, name))
    )
    if not entries:
        print(f"no case directories in {results}")
        return 1

    rows, tally = [], {PASS: 0, FAIL: 0, GAP: 0, SKIP: 0}
    for entry in entries:
        case = Case(results, downloads, entry)
        aborted_at, aborted_why = first_aborted(case)
        for label, runner, needs in CHECKS:
            try:
                status, detail = runner(case)
            except Exception as exc:  # a check must not take the run down with it
                status, detail = FAIL, f"check raised {type(exc).__name__}: {exc}"
            # A check whose input never got written because an earlier node stopped is not
            # evidence of anything. `check_failure_was_reported` is what judges that stop,
            # and it reads node 01, so it is never converted here.
            if (
                status == FAIL
                and aborted_at is not None
                and aborted_why is not None
                and NODE_ORDER.index(needs) >= NODE_ORDER.index(aborted_at)
            ):
                status = SKIP
                detail = (
                    f"node {aborted_at} exited {case.rc[aborted_at]}, so node {needs} "
                    f"produced nothing to check: {short_reason(aborted_why)}"
                )
            rows.append((entry, label, status, detail))
            tally[status] += 1

    width = max(len(r[1]) for r in rows)
    current = None
    for entry, label, status, detail in rows:
        if entry != current:
            print(f"\n{entry}")
            current = entry
        print(f"  {status:<4} {label:<{width}}  {detail}")

    print(
        f"\n{tally[PASS]} pass, {tally[FAIL]} fail, {tally[GAP]} gap, {tally[SKIP]} skipped"
    )
    return 1 if tally[FAIL] else 0


if __name__ == "__main__":
    sys.exit(main())
