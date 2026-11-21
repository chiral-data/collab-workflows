#!/usr/bin/env python3
"""Check a fixture run against a measurement of its input.

    python3 -I tests/check_fixtures.py <results-dir> <downloads-dir>

`tests/run_fixtures.sh` produces both directories; this script does the judging, so it can
be re-run over an existing run without paying for the run again.

Every check compares something a node *reported* against something measured independently:
`probe_entry.py` parses the mmCIF by hand, so the measurement does not share a parser -- or
a bug -- with the code under test. Each check returns `(status, detail)` and the detail
carries the numbers, so a failure reads as evidence rather than as an assertion.

An expectation marked `"status": "open"` in fixtures.json reports as GAP and does not fail
the run: a reproduced but unfixed defect, recorded rather than hidden.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from probe_entry import probe  # noqa: E402

PASS, FAIL, GAP, SKIP = "PASS", "FAIL", "GAP", "SKIP"

# Node 05's own collapse threshold (05-analyze/script.py: COLLAPSE_NM). Kept in step with it
# deliberately: this suite checks that the node's test fires, so it has to test the same thing.
COLLAPSE_NM = 0.01


def read_json(path, default=None):
    if not os.path.exists(path):
        return default
    try:
        with open(path) as handle:
            return json.load(handle)
    except Exception:
        return default


def read_text(path):
    if not os.path.exists(path):
        return ""
    with open(path, errors="replace") as handle:
        return handle.read()


def read_csv_rows(path):
    import csv

    if not os.path.exists(path):
        return []
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


class Case:
    """One fixture's run: its directories, its logs, and the measurement of its input."""

    def __init__(self, section, entry, results, downloads):
        self.section = section
        self.entry = entry
        self.results = results
        # A rejection case carries both a name and the entry it is built from; the name is
        # what identifies it, and what the driver named its directory after.
        self.name = entry.get("name") or entry.get("id")
        self.dir = os.path.join(results, section, slug(self.name))
        self.downloads = downloads
        self._truth = None

    def out(self, node, name):
        return os.path.join(self.dir, node, "outputs", name)

    def json(self, node, name, default=None):
        return read_json(self.out(node, name), default)

    def log(self, node):
        return read_text(os.path.join(self.dir, f"{node}.log"))

    def rc(self, node):
        text = read_text(os.path.join(self.dir, f"{node}.rc")).strip()
        return int(text) if text.isdigit() else None

    @property
    def truth(self):
        """The measured shape of this case's input.

        For a fixture that is the downloaded entry. For a downstream case it is the
        topology actually handed to the node -- which is what the node's own numbers have
        to match, and which also records whether the shape being tested for is present.
        """
        if self._truth is None:
            source = self.entry.get("source") or self.name
            candidates = [
                os.path.join(self.dir, "base", "topology.cif"),
                os.path.join(self.downloads, f"{source}.cif"),
            ]
            self._truth = {}
            for path in candidates:
                if os.path.exists(path):
                    self._truth = probe(path)
                    self._truth["measured_from"] = path
                    break
        return self._truth


def slug(name):
    return "".join(c if c.isalnum() or c in "-_." else "-" for c in name).strip("-")


# ----------------------------------------------------------------- unit section

def unit_result(case, which):
    results = read_json(os.path.join(case.results, "unit", "unit_results.json"))
    if not results:
        return FAIL, "tests/unit_checks.py wrote no results; it did not run"
    record = results.get(which)
    if not record:
        return FAIL, f"{which} is absent from unit_results.json"
    status = PASS if record.get("status") == "pass" else FAIL
    return status, record.get("detail", "")


def unit_a3m_path_helpers(case):
    return unit_result(case, "unit_a3m_path_helpers")


def unit_pairwise_pass(case):
    return unit_result(case, "unit_pairwise_pass")


def unit_select_frames(case):
    return unit_result(case, "unit_select_frames")


# -------------------------------------------------------------- fixtures section

def inventory_matches_probe(case):
    """Node 01 must describe the file it was given, not a file like it."""
    inventory = case.json("01", "system_inventory.json")
    if not inventory:
        return FAIL, f"node 01 wrote no inventory (exit {case.rc('01')})"
    truth = case.truth
    if not truth:
        return SKIP, "input was not probed"
    problems = []
    reported = set(inventory.get("chain_ids") or [])
    measured = set(truth.get("chain_ids") or [])
    if reported - measured:
        problems.append(f"chains {sorted(reported - measured)} are not in the file")
    proteins = set(inventory.get("protein_chains") or [])
    measured_proteins = set(truth.get("protein_chains") or [])
    if proteins - measured_proteins:
        problems.append(
            f"protein chains {sorted(proteins - measured_proteins)} carry no amino acids "
            "in the file"
        )
    if bool(inventory.get("has_nucleic_acid")) != bool(truth.get("has_nucleic_acid")):
        problems.append(
            f"has_nucleic_acid={inventory.get('has_nucleic_acid')} against "
            f"nucleic chains {truth.get('nucleic_chains')} measured"
        )
    # Node 01 counts the entity sequence; the probe counts residues the file resolves. The
    # entity is the longer of the two whenever a terminus or a loop is disordered -- 1CLL is
    # a 148-residue protein with 144 resolved -- so the entity count must be at least the
    # resolved count, and never below it. Reading these two numbers as the same quantity is
    # how "148 residues" ends up in a report about a 144-residue crystal structure, which is
    # why node 01 records where its number comes from.
    residues = inventory.get("n_polymer_residues")
    resolved = truth.get("n_polymer_residues") or 0
    if isinstance(residues, int) and residues < resolved:
        problems.append(
            f"{residues} polymer residues reported, {resolved} resolved in the file: the "
            "entity sequence cannot be shorter than the residues built from it"
        )
    if isinstance(residues, int) and residues > resolved * 2 and resolved:
        problems.append(
            f"{residues} polymer residues reported against {resolved} resolved: too large a "
            "gap to be disorder"
        )
    if problems:
        return FAIL, "; ".join(problems)
    unresolved = (residues - resolved) if isinstance(residues, int) else None
    return PASS, (
        f"{len(reported)} chain(s), protein {sorted(proteins)}, {residues} residue(s) in "
        f"the entity sequence against {resolved} resolved in the file"
        + (f" ({unresolved} unresolved)" if unresolved else "")
    )


def inventory_tokens_measured(case):
    """`n_tokens` must be a number. It was null on every run, from reading `n_tokens` for
    `num_tokens`, which also made the >10k GPU-fit warning unreachable."""
    inventory = case.json("01", "system_inventory.json")
    if not inventory:
        return FAIL, "node 01 wrote no inventory"
    tokens, source = inventory.get("n_tokens"), inventory.get("n_tokens_from")
    residues = (case.truth or {}).get("n_polymer_residues")
    if not isinstance(tokens, int):
        return FAIL, f"n_tokens={tokens!r}: the batch-size and GPU-fit advice cannot fire"
    if not source:
        return FAIL, f"n_tokens={tokens} with no n_tokens_from: provenance unrecorded"
    if residues and tokens < residues:
        return FAIL, f"n_tokens={tokens} is below the {residues} polymer residues measured"
    return PASS, f"n_tokens={tokens} from {source}, against {residues} polymer residue(s)"


def inventory_reports_ligands(case):
    """A non-polymer component must be reported as present.

    These are what bioemu2 samples as ligand tokens, and the report says whether the
    system had any. A haem, an AP5A or four calcium ions going unmentioned is how a
    reader concludes an apo ensemble was sampled holo, or the reverse.
    """
    inventory = case.json("01", "system_inventory.json")
    if not inventory:
        return FAIL, "node 01 wrote no inventory"
    truth = case.truth or {}
    hetero = truth.get("hetero_rows") or 0
    if not hetero:
        return SKIP, "no non-polymer rows in the file"
    if not inventory.get("has_ligands"):
        return FAIL, (
            f"{hetero} non-polymer atom row(s) in the file and has_ligands is "
            f"{inventory.get('has_ligands')!r}"
        )
    return PASS, f"has_ligands=true, against {hetero} non-polymer atom row(s) measured"


def node02_no_protein_chain(case):
    """No protein chain means no alignment, said out loud, with `msa/` shipped empty."""
    summary = case.json("02", "msa_summary.json")
    if not summary:
        return FAIL, f"node 02 wrote no summary (exit {case.rc('02')})"
    if summary.get("protein_chains"):
        return FAIL, f"protein_chains={summary['protein_chains']} for a nucleic-only entry"
    source = summary.get("source") or ""
    if "no protein chain" not in source:
        return FAIL, f"source={source!r} does not say why there is no alignment"
    msa_dir = case.out("02", "msa")
    if not os.path.isdir(msa_dir):
        return FAIL, "msa/ was not shipped; node 03 cannot tell empty from missing"
    staged = os.listdir(msa_dir)
    if staged:
        return FAIL, f"msa/ holds {staged} for a system with no protein chain"
    return PASS, f"source={source!r}, msa/ present and empty"


def node02_paired_expectation(case):
    """A paired alignment is due exactly when the system has >1 unique protein sequence.

    Only checkable on the server path: an uploaded alignment is unpaired by construction.
    `run_gpu_campaign.sh` runs that path for real; here it reports what it can.
    """
    summary = case.json("02", "msa_summary.json")
    if not summary:
        return FAIL, "node 02 wrote no summary"
    unique = (case.truth or {}).get("unique_protein_sequences")
    if unique is None:
        return SKIP, "input was not probed"
    if (summary.get("source") or "") != "ColabFold server":
        return SKIP, (
            f"ran on the {summary.get('source')!r} path, which is unpaired by construction; "
            "the server path is checked by run_gpu_campaign.sh"
        )
    due = unique > 1
    if bool(summary.get("paired")) != due:
        return FAIL, (
            f"{unique} unique protein sequence(s), so paired={due} is due, and node 02 "
            f"recorded paired={summary.get('paired')}"
        )
    if due:
        missing = [c for c, i in (summary.get("chains") or {}).items() if not i.get("paired_a3m")]
        if missing:
            return FAIL, f"paired=true but chain(s) {missing} have no paired file"
    return PASS, f"{unique} unique sequence(s) -> paired={summary.get('paired')}"


def node02_upload_path_records_no_pairing(case):
    """The upload path's cost must be recorded, not silent.

    An uploaded alignment set is unpaired, and bioemu2's pairmixer substitutes zeros for
    the missing paired feature. On a heteromer that is a real degradation of the run, so
    node 02 has to say `paired: false` and leave `paired_a3m: null` per chain rather than
    leaving a reader to assume the model got everything ColabFold would have given it.
    """
    summary = case.json("02", "msa_summary.json")
    if not summary:
        return FAIL, "node 02 wrote no summary"
    if (summary.get("source") or "") != "user-supplied A3M":
        return SKIP, f"ran on the {summary.get('source')!r} path"
    if summary.get("paired") is not False:
        return FAIL, f"upload path recorded paired={summary.get('paired')!r}, not false"
    chains = summary.get("chains") or {}
    wrong = [c for c, info in chains.items() if "paired_a3m" not in info or info["paired_a3m"]]
    if wrong:
        return FAIL, f"chain(s) {wrong} do not record paired_a3m: null"
    unique = (case.truth or {}).get("unique_protein_sequences")
    return PASS, (
        f"paired=false and paired_a3m null for all {len(chains)} chain(s); the file has "
        f"{unique} unique protein sequence(s), so this is a recorded cost, not an absence"
    )


def node03_dryrun_accepts_staged_alignment(case):
    """Node 03's dry run must find one alignment file per protein chain and say so.

    The check that matters is the pairing of the two numbers. The staging failure this
    suite exists for left the directory absent while node 02's summary named chains, and
    nothing downstream compared them.
    """
    rc = case.rc("03")
    log = case.log("03")
    summary = case.json("02", "msa_summary.json") or {}
    expected = list(summary.get("protein_chains") or [])
    if rc != 0:
        return FAIL, f"node 03 --dry-run exited {rc}: {last_error(log)}"
    if "Alignment staged" not in log:
        return FAIL, "the dry run did not report what it staged"
    for chain in expected:
        if f"{chain}.a3m" not in log:
            return FAIL, f"node 02 named chain {chain} and the dry run staged no {chain}.a3m"
    if not expected:
        return SKIP, "no protein chain in this entry"
    return PASS, f"{len(expected)} protein chain(s), one .a3m each: {staged_line(log)}"


def node03_dryrun_accepts_absent_alignment(case):
    """And must not demand one when none was due. A guard that refuses every empty `msa/`
    passes the rejection test above and breaks every nucleic-acid or ligand-only system."""
    rc = case.rc("03")
    if rc != 0:
        return FAIL, (
            f"node 03 --dry-run exited {rc} for a system with no protein chain: "
            f"{last_error(case.log('03'))}"
        )
    summary = case.json("02", "msa_summary.json") or {}
    if summary.get("protein_chains"):
        return SKIP, "this entry does have a protein chain"
    return PASS, "exit 0 with no alignment, because none was due"


def no_warnings(case):
    """The control: an ordinary monomer must raise nothing, so a warning elsewhere means
    something about that input, not about the pipeline."""
    noisy = {}
    for node in ("01", "02", "03"):
        lines = [l.strip() for l in case.log(node).splitlines() if "WARNING" in l]
        if lines:
            noisy[node] = lines
    if noisy:
        first = next(iter(noisy.items()))
        return FAIL, f"node {first[0]} warned: {first[1][0][:160]}"
    return PASS, "nodes 01-03 ran clean"


def staged_line(log):
    for line in log.splitlines():
        if "Alignment staged" in line:
            return line.strip()[:120]
    return ""


def last_error(log):
    lines = [l.strip() for l in log.splitlines() if l.strip()]
    for line in reversed(lines):
        if "ERROR" in line or "Error" in line:
            return line[:200]
    return lines[-1][:200] if lines else "(no output)"


# ------------------------------------------------------------ downstream section

def node05_ca_selection_excludes_ions(case):
    """The C-alpha set must be the protein's C-alphas and nothing else.

    A calcium ion is an atom named CA in a residue named CA, and mdtraj's `name CA`
    matches atom names. The per-residue RMSF list holds one entry per selected atom, so
    its length is the size of the selection -- the number that makes this visible.
    """
    sampled = sampled_summary(case)
    if not sampled:
        return FAIL, f"node 05 wrote no summary (exit {case.rc('05')})"
    selected = len(sampled.get("residue_rmsf") or [])
    truth = case.truth or {}
    ca, calcium = truth.get("ca_atoms"), truth.get("calcium_atoms")
    if ca is None:
        return SKIP, "the topology handed to node 05 was not probed"
    if not calcium:
        entry = probe(os.path.join(case.downloads, f"{case.entry.get('source') or case.name}.cif"))
        return SKIP, (
            f"the topology node 05 was given carries no calcium ({ca} C-alpha), so this "
            f"case cannot discriminate; the entry itself has "
            f"{entry.get('calcium_atoms', '?')} -- see inventory_reports_ligands"
        )
    if selected == ca + calcium:
        return FAIL, (
            f"{selected} atoms selected = {ca} C-alpha + {calcium} calcium ion(s): the ions "
            "are counted as C-alphas, which moves every RMSD, the RMSF and the medoid"
        )
    if selected != ca:
        return FAIL, f"{selected} atoms selected against {ca} C-alpha measured"
    return PASS, (
        f"{selected} C-alpha selected, matching the {ca} measured; the {calcium} calcium "
        "ion(s) were excluded"
    )


def node05_selection_recorded(case):
    """Which selection was used, in the output. 'C-alpha' alone was not falsifiable."""
    sampled = sampled_summary(case)
    if not sampled:
        return FAIL, "node 05 wrote no summary"
    selection = sampled.get("rmsd_atom_selection")
    allowed = (
        "protein and name CA",
        "name CA in residues with more than one atom",
        "all atoms",
    )
    if not selection:
        return FAIL, "no rmsd_atom_selection recorded"
    if selection not in allowed:
        return FAIL, f"unrecognised selection {selection!r}"
    if (case.truth or {}).get("protein_chains") and selection == "all atoms":
        return FAIL, "a protein system fell back to all atoms: the C-alpha selection failed"
    excluded = sampled.get("excluded_chains")
    return PASS, f"{selection!r}, excluded chains {excluded}"


def node05_rmsf_has_multiple_chains(case):
    """Per-residue RMSF must carry a chain label per entry, or the report cannot split it."""
    sampled = sampled_summary(case)
    if not sampled:
        return FAIL, "node 05 wrote no summary"
    rmsf = sampled.get("residue_rmsf") or []
    if not rmsf:
        return FAIL, "no per-residue RMSF written"
    if any("chain" not in entry for entry in rmsf):
        return FAIL, "some RMSF entries carry no chain label"
    chains = sorted({entry["chain"] for entry in rmsf})
    measured = (case.truth or {}).get("protein_chains") or []
    if len(chains) < 2:
        return FAIL, (
            f"RMSF spans chain(s) {chains} for a topology with protein chains {measured}: "
            "a homodimer drawn as one series strokes back to the left edge where the second "
            "chain restarts, which reads as a flexible region"
        )
    return PASS, f"{len(rmsf)} entries over chain(s) {chains}, each labelled"


def node06_rmsf_traces_per_chain(case):
    """And the report must actually split it, by reading the labels out of its own DATA."""
    report = read_text(case.out("06", "report.html"))
    if not report:
        return FAIL, "node 06 wrote no report"
    data = report_data(report)
    if data is None:
        return FAIL, "could not parse the report's DATA object"
    ensembles = data.get("ensembles") or []
    if not ensembles:
        return FAIL, "DATA carries no ensembles"
    rmsf = ensembles[0].get("rmsf") or []
    chains = sorted({entry.get("chain") for entry in rmsf if isinstance(entry, dict)})
    if len(chains) < 2:
        return FAIL, f"the report's RMSF data spans chain(s) {chains}"
    if "chains.forEach" not in report or "DASHES" not in report:
        return FAIL, "the report has per-chain data but builds a single trace from it"
    return PASS, f"{len(rmsf)} RMSF points over chain(s) {chains}, one trace each"


def node05_matrix_dropped_above_cap(case):
    """Above the cap the matrix is not retained -- that part was never the bug."""
    summary = case.json("05", "ensemble_summary.json")
    if not summary:
        return FAIL, "node 05 wrote no summary"
    frames = (sampled_summary(case) or {}).get("n_frames")
    if summary.get("rmsd_matrix") is not None:
        return FAIL, f"an n^2 matrix was written for {frames} frames"
    return PASS, f"{frames} frames: matrix not retained, as designed"


def node05_medoid_measured_over_all_pairs(case):
    """The medoid must stay measured where the matrix is dropped.

    Node 05 used to set `medoid = 0` above 200 frames while keeping the column name
    `rmsd_to_medoid_nm`, so the ensemble's reference frame silently became "whichever frame
    was sampled first". Two things prove it is real now: the spread covers every pair, and
    the per-frame RMSD-to-medoid column is zero exactly at the frame named as the medoid.
    """
    sampled = sampled_summary(case)
    if not sampled:
        return FAIL, "node 05 wrote no summary"
    frames = sampled.get("n_frames") or 0
    medoid = sampled.get("medoid_frame")
    pairs = (sampled.get("pairwise_rmsd") or {}).get("n_pairs")
    expected_pairs = frames * (frames - 1) // 2
    if pairs != expected_pairs:
        return FAIL, f"{pairs} pair(s) measured for {frames} frames, expected {expected_pairs}"
    rows = [r for r in read_csv_rows(case.out("05", "ensemble_metrics.csv"))
            if r.get("ensemble") == "sampled"]
    if not rows:
        return FAIL, "no per-frame rows to check the medoid against"
    try:
        distances = {int(r["frame"]): float(r["rmsd_to_medoid_nm"]) for r in rows}
    except (KeyError, ValueError) as exc:
        return FAIL, f"could not read rmsd_to_medoid_nm: {exc}"
    zeros = [f for f, d in distances.items() if d < 1e-6]
    if zeros != [medoid]:
        return FAIL, (
            f"medoid_frame={medoid} but rmsd_to_medoid_nm is zero at {zeros}"
        )
    argmin = min(distances, key=lambda f: distances[f])
    return PASS, (
        f"{frames} frames, all {pairs} pair(s) measured; rmsd_to_medoid is 0 at frame "
        f"{argmin} == medoid_frame, so the reference frame is real"
    )


def node05_collapse_warned(case):
    """Frames that differ only by a rigid rotation are one conformation, and must be
    reported as one. The fixture rotates every frame randomly, so a node that forgot to
    superpose measures a large spread here and calls a collapsed ensemble healthy."""
    sampled = sampled_summary(case)
    if not sampled:
        return FAIL, "node 05 wrote no summary"
    maximum = (sampled.get("pairwise_rmsd") or {}).get("max_nm")
    log = case.log("05")
    if maximum is None:
        return FAIL, "no pairwise RMSD recorded"
    # The same threshold node 05 warns at. An XTC quantises coordinates to 0.001 nm, so one
    # conformation written twice comes back differing by about that much: measured 0.0018-0.0022
    # nm here, which is why the node's test is 0.01 nm and not float epsilon.
    if maximum > COLLAPSE_NM:
        return FAIL, (
            f"max pairwise RMSD {maximum} nm over identical-but-rotated frames, above the "
            f"{COLLAPSE_NM} nm collapse threshold: either the superposition did not happen, so "
            "rigid-body drift is being reported as conformational spread, or the frames are "
            "not the fixture's repeated conformation"
        )
    if "collapsed" not in log:
        return FAIL, f"max pairwise RMSD {maximum} nm and the log says nothing about it"
    return PASS, f"max pairwise RMSD {maximum} nm, and the log says the ensemble collapsed"


def node05_single_cluster(case):
    """One conformation must cluster as one, whatever the rotations."""
    sampled = sampled_summary(case)
    if not sampled:
        return FAIL, "node 05 wrote no summary"
    clusters = sampled.get("n_clusters")
    if clusters != 1:
        return FAIL, f"{clusters} clusters for an ensemble of one repeated conformation"
    return PASS, f"1 cluster at {sampled.get('rmsd_cluster_threshold_nm')} nm"


def node05_pass_rate_null_when_unchecked(case):
    """Unchecked is not clean. `n_checked = 0` must give `pass_rate: null`, never 0 or 1."""
    sampled = sampled_summary(case)
    if not sampled:
        return FAIL, "node 05 wrote no summary"
    physicality = sampled.get("physicality") or {}
    records = case.json("05", "physicality.json") or read_json(
        os.path.join(case.dir, "base", "physicality.json")
    ) or []
    statuses = sorted({r.get("status") for r in records if isinstance(r, dict)})
    if statuses and statuses != ["not_checked"]:
        return SKIP, f"this case's frames are {statuses}, not not_checked"
    if physicality.get("n_checked"):
        return FAIL, f"n_checked={physicality['n_checked']} where every frame is not_checked"
    if physicality.get("pass_rate") is not None:
        return FAIL, f"pass_rate={physicality['pass_rate']!r} with nothing checked"
    return PASS, f"n_checked=0, pass_rate=null over {len(records)} not_checked frame(s)"


def node06_pass_rate_not_shown_as_clean(case):
    """And the report must say so, rather than rendering a clean-looking number."""
    report = read_text(case.out("06", "report.html"))
    if not report:
        return FAIL, "node 06 wrote no report"
    sampled = sampled_summary(case) or {}
    physicality = sampled.get("physicality") or {}
    if physicality.get("pass_rate") is not None:
        return SKIP, f"pass_rate={physicality['pass_rate']} in this case"
    if "not checked" not in report:
        return FAIL, "pass_rate is null and the report does not say the frames were not checked"
    reason = (case.json("03", "sampling_manifest.json") or read_json(
        os.path.join(case.dir, "base", "sampling_manifest.json")
    ) or {}).get("physicality_not_checked_because")
    if reason and reason not in report:
        return FAIL, f"the report does not carry the recorded reason: {reason!r}"
    return PASS, f"shown as 'not checked'" + (f", because {reason}" if reason else "")


def node04_passthrough_outputs_complete(case):
    """Propagation off still has to write every declared output.

    Silva collects outputs with `nullglob`, so a declared output that is never written
    vanishes without an error and the next node sees a missing file instead of a disabled
    feature. Node 04's identity path is the one that has to be complete.
    """
    manifest = case.json("04", "propagation_manifest.json")
    if not manifest:
        return FAIL, f"node 04 wrote no manifest (exit {case.rc('04')})"
    if manifest.get("enabled"):
        return SKIP, "propagation ran in this case"
    declared = ("propagated_samples.xtc", "propagated_topology.cif",
                "propagated_physicality.json", "propagation_manifest.json")
    missing = [name for name in declared if not os.path.exists(case.out("04", name))]
    if missing:
        return FAIL, f"declared output(s) {missing} were not written on the identity path"
    return PASS, f"all {len(declared)} declared output(s) present: {manifest.get('reason')}"


def node04_passthrough_physicality_carried(case):
    """The identity path must carry node 03's per-frame physicality forward unchanged.

    Node 05 reads it to label frames. Re-deriving it here, or dropping it, is the defect
    class CLAUDE.md names: an upstream fact recomputed instead of read.
    """
    carried = case.json("04", "propagated_physicality.json")
    if carried is None:
        return FAIL, "no propagated_physicality.json"
    manifest = case.json("04", "propagation_manifest.json") or {}
    if manifest.get("enabled"):
        return SKIP, "propagation ran in this case"
    original = read_json(os.path.join(case.dir, "base", "physicality.json"))
    if original is None:
        return SKIP, "no input physicality to compare against"
    if carried != original:
        return FAIL, (
            f"{len(carried)} record(s) carried against {len(original)} in, and they differ: "
            "the identity path changed the frames' physicality"
        )
    return PASS, f"{len(carried)} record(s) carried forward byte-for-byte"


def node05_passthrough_not_treated_as_lagged(case):
    """An identity copy is not a second, time-lagged ensemble, and must not be charted as one."""
    summary = case.json("05", "ensemble_summary.json")
    if not summary:
        return FAIL, "node 05 wrote no summary"
    if summary.get("propagation_enabled"):
        return FAIL, "node 04 ran an identity pass and node 05 reported propagation enabled"
    labels = sorted((summary.get("ensembles") or {}).keys())
    if "propagated" in labels:
        return FAIL, f"ensembles {labels}: the identity copy was analysed as a second ensemble"
    return PASS, f"ensembles {labels}, propagation_enabled=false"


def node05_refuses_missing_manifest(case):
    """A missing propagation manifest is not disabled propagation.

    It means node 04 did not run, or its outputs were not staged -- the same staging
    failure this suite exists for, one node further down. Reading it as "propagation off"
    turns a broken DAG into a quiet half-run.
    """
    rc = case.rc("05")
    log = case.log("05")
    if rc == 0:
        return FAIL, "node 05 exited 0 with no propagation manifest staged"
    if "propagation_manifest.json" not in log:
        return FAIL, f"node 05 failed without naming the missing file: {last_error(log)}"
    return PASS, f"exit {rc}: {last_error(log)}"


def node05_refuses_lost_ensemble(case):
    """And propagation enabled with no propagated ensemble is not disabled propagation either."""
    rc = case.rc("05")
    log = case.log("05")
    if rc == 0:
        return FAIL, (
            "node 05 exited 0 where the manifest says propagation ran and the propagated "
            "ensemble is absent"
        )
    if "propagated" not in log:
        return FAIL, f"node 05 failed without saying what was missing: {last_error(log)}"
    return PASS, f"exit {rc}: {last_error(log)}"


def node05_filler_counted_honestly(case):
    """The viewer's filler frames must be counted as filler.

    Cluster centres come first and consecutive frames pad the rest. Reporting the pad as a
    distinct conformation overstates what the ensemble found.
    """
    summary = case.json("05", "ensemble_summary.json")
    if not summary:
        return FAIL, "node 05 wrote no summary"
    frames = summary.get("viewer_frames") or []
    from_clusters = summary.get("viewer_frames_from_clusters")
    sampled = sampled_summary(case) or {}
    clusters = sampled.get("n_clusters")
    if from_clusters is None:
        return FAIL, "viewer_frames_from_clusters was not written"
    if from_clusters > len(frames):
        return FAIL, f"{from_clusters} claimed from clusters out of {len(frames)} frames"
    if clusters is not None and from_clusters > clusters:
        return FAIL, f"{from_clusters} cluster centres claimed where {clusters} cluster(s) exist"
    if len(frames) != len(set(frames)):
        return FAIL, f"viewer_frames repeats a frame: {frames}"
    if clusters == 1 and from_clusters != 1:
        return FAIL, f"one cluster and {from_clusters} centres claimed"
    return PASS, (
        f"{len(frames)} viewer frame(s), {from_clusters} real centre(s) and "
        f"{len(frames) - from_clusters} filler, over {clusters} cluster(s)"
    )


def node06_caption_matches_selection(case):
    """The 3D viewer's caption must state what the viewer is showing.

    The frames are cluster centres plus filler, superposed on the medoid. A caption that
    says "cluster centres" alone is wrong about the filler, and one that says nothing about
    the superposition leaves a reader to read alignment into the animation.
    """
    report = read_text(case.out("06", "report.html"))
    summary = case.json("05", "ensemble_summary.json") or {}
    if not report:
        return FAIL, "node 06 wrote no report"
    frames = summary.get("viewer_frames") or []
    from_clusters = summary.get("viewer_frames_from_clusters")
    medoid = summary.get("viewer_frames_superposed_on_frame")
    needed = [str(len(frames))]
    if from_clusters is not None and from_clusters < len(frames):
        needed.append(str(len(frames) - from_clusters))
    missing = [n for n in needed if n not in report]
    if missing:
        return FAIL, f"the viewer caption does not carry {missing}"
    if medoid is not None and "superposed" not in report:
        return FAIL, "the caption does not say the frames are superposed"
    return PASS, (
        f"caption carries {len(frames)} frame(s), "
        f"{len(frames) - (from_clusters or 0)} filler, superposed on frame {medoid}"
    )


def node06_report_built(case):
    """One self-contained page, carrying the caveats that keep its numbers readable."""
    path = case.out("06", "report.html")
    report = read_text(path)
    if not report:
        return FAIL, f"no report.html (node 06 exited {case.rc('06')})"
    required = {
        "the DATA object": "const DATA =",
        "the caveat": "Frame counts are not populations",
    }
    missing = [label for label, text in required.items() if text not in report]
    if missing:
        return FAIL, f"report is missing {', '.join(missing)}"
    sampled = sampled_summary(case) or {}
    selection = sampled.get("rmsd_atom_selection")
    if selection and selection not in report:
        return FAIL, f"the report does not state the selection ({selection!r})"
    return PASS, f"{round(len(report) / 1024)} KB, with the selection and the caveats"


def sampled_summary(case):
    summary = case.json("05", "ensemble_summary.json")
    if not summary:
        return None
    return (summary.get("ensembles") or {}).get("sampled") or {}


def report_data(report):
    """The report's own DATA object, so a check reads what the page plots rather than
    searching the prose -- node 06 explains its charts in text, so every label appears in
    the file whether the chart is built or not."""
    marker = "const DATA = "
    start = report.find(marker)
    if start < 0:
        return None
    start += len(marker)
    depth, index = 0, start
    while index < len(report):
        if report[index] == "{":
            depth += 1
        elif report[index] == "}":
            depth -= 1
            if depth == 0:
                try:
                    return json.loads(report[start : index + 1])
                except Exception:
                    return None
        index += 1
    return None


# ------------------------------------------------------------ rejection section

def check_rejection(case):
    """An input the pipeline must refuse, or must accept -- the counter-cases matter as
    much, since a guard that refuses everything passes every refusal test."""
    entry = case.entry
    stage = entry.get("stage", "")
    node = "03" if "node03" in stage else "04"
    rc = case.rc(node)
    log = case.log(node)
    if rc is None:
        return FAIL, f"node {node} did not run"
    if entry.get("expect_ok"):
        if rc != 0:
            return FAIL, f"expected acceptance, got exit {rc}: {last_error(log)}"
        return PASS, f"accepted, exit 0"
    needle = entry.get("expect_error") or ""
    if rc == 0:
        return FAIL, f"expected a refusal naming {needle!r}, got exit 0"
    if needle and needle.lower() not in log.lower():
        return FAIL, f"refused with exit {rc} but not for the stated reason: {last_error(log)}"
    return PASS, f"exit {rc}, naming {needle!r}: {last_error(log)}"


CHECKS = {
    "unit_a3m_path_helpers": unit_a3m_path_helpers,
    "unit_pairwise_pass": unit_pairwise_pass,
    "unit_select_frames": unit_select_frames,
    "inventory_matches_probe": inventory_matches_probe,
    "inventory_tokens_measured": inventory_tokens_measured,
    "inventory_reports_ligands": inventory_reports_ligands,
    "node02_no_protein_chain": node02_no_protein_chain,
    "node02_paired_expectation": node02_paired_expectation,
    "node02_upload_path_records_no_pairing": node02_upload_path_records_no_pairing,
    "node03_dryrun_accepts_staged_alignment": node03_dryrun_accepts_staged_alignment,
    "node03_dryrun_accepts_absent_alignment": node03_dryrun_accepts_absent_alignment,
    "no_warnings": no_warnings,
    "node05_ca_selection_excludes_ions": node05_ca_selection_excludes_ions,
    "node05_selection_recorded": node05_selection_recorded,
    "node05_rmsf_has_multiple_chains": node05_rmsf_has_multiple_chains,
    "node06_rmsf_traces_per_chain": node06_rmsf_traces_per_chain,
    "node05_matrix_dropped_above_cap": node05_matrix_dropped_above_cap,
    "node05_medoid_measured_over_all_pairs": node05_medoid_measured_over_all_pairs,
    "node05_collapse_warned": node05_collapse_warned,
    "node05_single_cluster": node05_single_cluster,
    "node05_pass_rate_null_when_unchecked": node05_pass_rate_null_when_unchecked,
    "node06_pass_rate_not_shown_as_clean": node06_pass_rate_not_shown_as_clean,
    "node04_passthrough_outputs_complete": node04_passthrough_outputs_complete,
    "node04_passthrough_physicality_carried": node04_passthrough_physicality_carried,
    "node05_passthrough_not_treated_as_lagged": node05_passthrough_not_treated_as_lagged,
    "node05_refuses_missing_manifest": node05_refuses_missing_manifest,
    "node05_refuses_lost_ensemble": node05_refuses_lost_ensemble,
    "node05_filler_counted_honestly": node05_filler_counted_honestly,
    "node06_caption_matches_selection": node06_caption_matches_selection,
    "node06_report_built": node06_report_built,
}


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    results, downloads = sys.argv[1], sys.argv[2]
    spec = read_json(os.path.join(HERE, "fixtures.json"))
    if not spec:
        print("could not read tests/fixtures.json")
        return 2

    rows, tally = [], {PASS: 0, FAIL: 0, GAP: 0, SKIP: 0}

    def record(section, name, label, status, detail):
        rows.append((section, name, label, status, detail))
        tally[status] += 1

    for section in ("unit", "fixtures", "rejections", "downstream"):
        for entry in spec.get(section) or []:
            case = Case(section, entry, results, downloads)
            if section == "rejections":
                expectations = [{"check": "rejection", "status": entry.get("status")}]
            else:
                expectations = entry.get("expect") or []
            for expectation in expectations:
                label = expectation.get("check")
                runner = check_rejection if label == "rejection" else CHECKS.get(label)
                if runner is None:
                    record(section, case.name, label, FAIL,
                           "no such check in check_fixtures.py")
                    continue
                try:
                    status, detail = runner(case)
                except Exception as exc:
                    status, detail = FAIL, f"check raised {type(exc).__name__}: {exc}"
                if status == FAIL and expectation.get("status") == "open":
                    status, detail = GAP, f"known gap: {detail}"
                record(section, case.name, label, status, detail)

    width = max(len(r[2] or "") for r in rows)
    section, name = None, None
    for row_section, row_name, label, status, detail in rows:
        if row_section != section:
            print(f"\n=== {row_section}")
            section, name = row_section, None
        if row_name != name:
            print(f"  {row_name}")
            name = row_name
        print(f"    {status:<4} {label:<{width}}  {detail}")

    print(
        f"\n{tally[PASS]} pass, {tally[FAIL]} fail, {tally[GAP]} gap, {tally[SKIP]} skipped"
    )
    if tally[GAP]:
        print("Gaps are reproduced, unfixed defects recorded in fixtures.json; they do not fail.")
    return 1 if tally[FAIL] else 0


if __name__ == "__main__":
    sys.exit(main())
