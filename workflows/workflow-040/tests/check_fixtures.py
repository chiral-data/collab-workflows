#!/usr/bin/env python3
"""Judge a staged fixture run. Reads run.json and the tree beside it; runs nothing.

Separate from run_fixtures.py on purpose. A runner that also judged could only ever
report pass or fail, and when something failed there was nothing left on disk to look at
-- which is how 037's first suite went. Here the run leaves a tree and this file reads it,
so a failure can be re-judged, and re-judged with a changed check, without re-running a
single node.

Three statuses, and the difference between the last two is the whole point:

  PASS   the claim was measured and holds. The detail says what was measured.
  FAIL   the claim was measured and does not hold. Phase 1 CI must have none.
  GAP    the claim cannot be measured yet, because the behaviour it is about is Phase 2
         and Phase 2 is not built. Not a pass. Listed, counted and named every run.

A check NEVER returns PASS for an artifact it could not find. That is the failure mode
that matters most here: a suite that reports green because it looked in the wrong place is
worse than no suite, and 029's interface filters read as unreachable for exactly that kind
of reason. So every "not found" is a FAIL, and every PASS detail carries the numbers it
was decided on.

Non-vacuity
-----------
Each entry in CHECKS carries, beside its function, a note saying what would make it fail
-- because a check nothing can break is a check that measures nothing. The report prints
every check's invocation count and that note, and flags any check named in fixtures.json
with no implementation, or implemented and never invoked.

Run
---
  python3 tests/check_fixtures.py <run-dir>
  python3 tests/check_fixtures.py <run-dir> --verbose     every detail, not just failures
"""
import glob
import importlib.util
import json
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

FAIL, GAP, PASS = "FAIL", "GAP", "PASS"
RANK = {FAIL: 0, GAP: 1, PASS: 2}

NO_COMPLEX_STAGED = (
    "No complex reached this step, so the claim cannot be measured here. Phase 2 IS "
    "built; what this means is that this case did not ask for the complex path -- "
    "PARAM_SCORE_INTERFACE in its params -- or node 01b refused before node 02 could "
    "write one. Either way it is a gap, not a pass."
)


def load_probe():
    """probe_design.py as a module. The independent measurer; see its docstring (§8.4)."""
    spec = importlib.util.spec_from_file_location(
        "probe_design", os.path.join(HERE, "probe_design.py")
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PROBE = load_probe()


def natural_key(name):
    return [int(p) if p.isdigit() else p for p in re.split(r"(\d+)", name)]


def near(got, want, tolerance):
    if got is None:
        return False
    return abs(float(got) - float(want)) <= float(tolerance)


class View:
    """Read-only view over one case's staged tree."""

    def __init__(self, record):
        self.record = record
        self.dir = record["dir"]
        self.id = record["id"]
        self.phase = record.get("phase")
        self.tokens = record.get("tokens") or {}
        self._probe = None

    # -- raw ---------------------------------------------------------------------

    def path(self, *parts):
        return os.path.join(self.dir, *parts)

    def text(self, *parts):
        full = self.path(*parts)
        return open(full, errors="replace").read() if os.path.exists(full) else None

    def json(self, *parts):
        raw = self.text(*parts)
        if raw is None:
            return None
        try:
            return json.loads(raw)
        except ValueError:
            return None

    def files(self, *parts):
        where = self.path(*parts)
        if not os.path.isdir(where):
            return []
        return sorted(n for n in os.listdir(where)
                      if os.path.isfile(os.path.join(where, n)))

    def size(self, *parts):
        """Bytes, or None if absent. A zero-byte file and a missing one are different."""
        full = self.path(*parts)
        return os.path.getsize(full) if os.path.exists(full) else None

    def step(self, name):
        for entry in self.record.get("steps") or []:
            if entry["step"] == name:
                return entry
        return None

    def rc(self, name):
        step = self.step(name)
        return None if step is None else step["rc"]

    def log(self, name):
        step = self.step(name)
        if step is None:
            return ""
        return (step.get("stdout") or "") + "\n" + (step.get("stderr") or "")

    def all_logs(self):
        return "\n".join(self.log(s["step"]) for s in (self.record.get("steps") or []))

    def failures(self):
        return [s for s in (self.record.get("steps") or []) if s["rc"] != 0]

    # -- named artifacts ---------------------------------------------------------

    def wants_interface(self):
        """Did this case ASK for the complex path? Read from its params, like the harness."""
        value = (self.record.get("params") or {}).get("PARAM_SCORE_INTERFACE", "")
        return str(value).strip().lower() == "true"

    def inventory(self, step="n1"):
        return self.json(step, "outputs", "design_inventory.json")

    def inputs_decl(self, step="n2"):
        return self.json(step, "outputs", "prediction_inputs.json")

    def manifest(self):
        for where in (("n3", "outputs"), ("pred",)):
            got = self.json(*where, "prediction_manifest.json")
            if got is not None:
                return got
        return None

    def truth(self):
        return self.json("pred", "truth.json")

    def scores(self):
        return self.json("n4", "outputs", "design_scores.json")

    def report(self):
        return self.text("n5", "outputs", "report.html")

    def static_html(self):
        """report.html with every <script> block removed.

        What is left is what a reader sees with JavaScript off, on a machine with no
        network. 037's report put its provenance in a JS renderer and lost all of it
        offline (cb418341), so 'the report says X' has to mean the static page says X.
        """
        html = self.report()
        if html is None:
            return None
        return re.sub(r"<script\b.*?</script>", "", html, flags=re.S | re.I)

    def fasta(self):
        return self.text("n5", "outputs", "passing_sequences.fasta")

    def metrics_csv(self):
        return self.text("n5", "outputs", "validation_metrics.csv")

    def design(self, name=None, step="n4"):
        scored = self.scores()
        if not scored:
            return None
        rows = scored.get("designs") or []
        if name is None:
            return rows[0] if rows else None
        for row in rows:
            if row.get("design") == name:
                return row
        return None

    # -- the independent measurement ---------------------------------------------

    def probe(self, where="g"):
        """probe_design's reading of the generated structures. Shares no code with a node."""
        if self._probe is None:
            self._probe = {}
        if where not in self._probe:
            paths = sorted(glob.glob(self.path(where, "*.cif")))
            self._probe[where] = [PROBE.probe(p) for p in paths]
        return self._probe[where]

    def probe_chain(self, chain, where="g", which=0):
        for entry in self.probe(where)[which]["chains"]:
            if entry["chain"] == chain:
                return entry
        return None


# ===========================================================================
# ingest / node 01
# ===========================================================================

def c_inventory_designed_chain(c, spec):
    inv = c.inventory()
    if not inv:
        return FAIL, "no design_inventory.json was written"
    row = inv["designs"][0]
    want_chain, want_n = spec["chain"], spec["n_residues"]
    # The count is taken from probe_design, not from the inventory it is checking.
    measured = c.probe_chain(want_chain)
    if measured is None:
        return FAIL, f"probe_design found no chain {want_chain} in the generated file"
    detail = (
        f"inventory says binder_chain={row['binder_chain']} with "
        f"n_binder_residues={row['n_binder_residues']}; probe_design measures chain "
        f"{want_chain} as {measured['n_residues']} residues, sequence "
        f"{measured['sequence'][:16]}..., and the inventory's binder_sequence is "
        f"{row['binder_sequence'][:16]}..."
    )
    if row["binder_chain"] != want_chain:
        return FAIL, f"expected chain {want_chain}. " + detail
    if row["n_binder_residues"] != want_n or measured["n_residues"] != want_n:
        return FAIL, f"expected {want_n} residues. " + detail
    if row["binder_sequence"] != measured["sequence"]:
        return FAIL, "the inventory's sequence is not the one in the file. " + detail
    return PASS, detail


def c_inventory_chain_source(c, spec):
    inv = c.inventory()
    if not inv:
        return FAIL, "no design_inventory.json was written"
    got = inv["designs"][0]["designed_chain_source"]
    detail = (f"designed_chain_source={got!r}, workflow-level sources="
              f"{inv.get('designed_chain_sources')}, metrics_present="
              f"{inv.get('metrics_present')}")
    return (PASS if got == spec["source"] else FAIL), (
        detail if got == spec["source"] else f"expected {spec['source']!r}. " + detail
    )


def c_fallback_stated_in_log(c, spec):
    log = c.log("n1")
    inv = c.inventory()
    source = (inv or {}).get("designs", [{}])[0].get("designed_chain_source")
    says_source = source and source in log
    names_file = "design_0_generated.cif" in log
    warns = "WARNING" in log or "fallback" in log.lower() or "assum" in log.lower()
    detail = (f"node 01's log mentions the source {source!r}: {bool(says_source)}; names "
              f"the file: {names_file}; warns in words: {warns}")
    if says_source and names_file and warns:
        return PASS, detail
    return FAIL, "a stated fallback has to be stated. " + detail


def c_source_reaches_report(c, spec):
    html = c.static_html()
    if html is None:
        return FAIL, "no report.html was written"
    inv = c.inventory()
    source = (inv or {}).get("designs", [{}])[0].get("designed_chain_source")
    csv = c.metrics_csv() or ""
    in_html = bool(source) and source in html
    in_csv = bool(source) and source in csv
    detail = (f"source {source!r} appears in the STATIC report.html: {in_html}; in "
              f"validation_metrics.csv: {in_csv}")
    return (PASS if in_html and in_csv else FAIL), detail


def c_inventory_target_chains(c, spec):
    inv = c.inventory()
    if not inv:
        return FAIL, "no design_inventory.json was written"
    row = inv["designs"][0]
    got = [t["chain"] for t in row["target_chains"]]
    binder = row["binder_chain"]
    present = [ch["chain"] for ch in c.probe()[0]["chains"]]
    expected_from_file = [ch for ch in present if ch != binder]
    detail = (f"inventory target chains {got}; probe_design sees chains {present} in the "
              f"file with binder {binder}, so the targets in the file are "
              f"{expected_from_file}")
    if got != spec["chains"]:
        return FAIL, f"expected {spec['chains']}. " + detail
    if sorted(got) != sorted(expected_from_file):
        return FAIL, "the inventory's targets are not the file's. " + detail
    for recorded in row["target_chains"]:
        measured = c.probe_chain(recorded["chain"])
        if recorded["sequence"] != measured["sequence"]:
            return FAIL, (f"chain {recorded['chain']}'s recorded sequence is not the "
                          f"file's. " + detail)
    return PASS, detail + "; every recorded target sequence matches the file"


def c_inventory_chain_order_recorded(c, spec):
    decl = c.inputs_decl()
    if not decl:
        return FAIL, "no prediction_inputs.json was written"
    rows = decl.get("inputs") or []
    if not rows:
        return FAIL, "prediction_inputs.json declares no inputs"
    bad = [r["design"] for r in rows
           if not r.get("chain_order") or
           len(r["chain_order"]) != len(r.get("chain_lengths") or [])]
    manifest = c.manifest() or {}
    carried = [p.get("chain_order") for p in (manifest.get("predictions") or [])]
    detail = (f"node 02 recorded chain_order={rows[0]['chain_order']} with "
              f"chain_lengths={rows[0]['chain_lengths']} for {len(rows)} input(s); the "
              f"manifest carries {carried or 'nothing (no prediction step in this case)'}; "
              f"chain_order_source={manifest.get('chain_order_source')}")
    if bad:
        return FAIL, f"order missing or mismatched for {bad}. " + detail
    return PASS, detail


def c_inventory_records_gaps(c, spec):
    inv = c.inventory()
    if not inv:
        return FAIL, "no design_inventory.json was written"
    row = inv["designs"][0]
    recorded = [g for t in row["target_chains"] for g in (t.get("gaps") or [])]
    measured = [g for ch in c.probe()[0]["chains"]
                if ch["chain"] != row["binder_chain"] for g in ch["gaps"]]
    detail = (f"inventory records {len(recorded)} target gap(s) {recorded} and a "
              f"workflow-level n_target_numbering_gaps="
              f"{inv.get('n_target_numbering_gaps')}; probe_design independently measures "
              f"{len(measured)} gap(s) {measured}")
    if len(recorded) < spec["min_gaps"]:
        return FAIL, f"expected at least {spec['min_gaps']}. " + detail
    if len(recorded) != len(measured):
        return FAIL, "the inventory and the file disagree on how many. " + detail
    return PASS, detail


def c_fusion_disclosed(c, spec):
    """Every numbering break in a folded target must be disclosed as a fusion, by name.

    029 trimmed a gapped receptor to the segment carrying the epitope and measured ~42
    pLDDT against 81 for the whole chain (8e4772a7), so 040 folds the chain WHOLE and the
    junction it creates is a real defect kept knowingly. A defect kept knowingly has to be
    stated, and stated where the numbers are read.

    What this does NOT accept: the string "fus" appearing in prediction_inputs.json. The
    payload carries a `fusion` key whether or not anything was fused, so the key name
    alone satisfied the first version of this check -- which is §8.4's case, a check that
    passes for the wrong reason. So this counts the gaps node 01 measured and requires a
    disclosure naming each one's boundary residues.
    """
    decl = c.inputs_decl() or {}
    modes = sorted({r.get("mode") for r in decl.get("inputs") or []})
    if not modes or modes == ["mono"]:
        return GAP, (f"modes={modes or 'none declared'}; the target is only folded on the "
                     f"complex path, so no fusion can be disclosed. " + NO_COMPLEX_STAGED)

    inv = c.inventory()
    if not inv:
        return FAIL, "no design_inventory.json was written, so there are no gaps to match"
    gaps = [{"chain": chain["chain"], **gap}
            for row in inv["designs"] if not row.get("skipped")
            for chain in (row.get("target_chains") or [])
            for gap in (chain.get("gaps") or [])]
    if not gaps:
        return FAIL, ("the inventory records no target numbering gap, so this case is not "
                      "the shape it says it is and the disclosure cannot be checked "
                      "(§8.14). A gapped target is the whole point of the case.")

    target = (decl.get("target") or {})
    disclosed = [str(line) for line in (target.get("fusion") or [])]
    blob = " ".join(disclosed).lower()
    missing = [g for g in gaps
               if not (str(g["after"]) in blob and str(g["before"]) in blob
                       and str(g["chain"]).lower() in blob)]
    detail = (f"modes={modes}; the inventory measured {len(gaps)} numbering gap(s) "
              f"{[(g['chain'], g['after'], g['before']) for g in gaps]}; node 02 declared "
              f"{len(disclosed)} fusion line(s): {disclosed}")
    if not disclosed:
        return FAIL, ("the target was folded and nothing was disclosed. " + detail)
    if missing:
        return FAIL, (
            f"{len(missing)} gap(s) are folded over without being named: "
            f"{[(g['chain'], g['after'], g['before']) for g in missing]}. A disclosure "
            f"that does not say WHERE the junction is cannot be acted on. " + detail
        )
    return PASS, detail


def c_target_not_trimmed(c, spec):
    inv = c.inventory()
    if not inv:
        return FAIL, "no design_inventory.json was written"
    row = inv["designs"][0]
    lines = []
    for recorded in row["target_chains"]:
        measured = c.probe_chain(recorded["chain"])
        lines.append(
            f"chain {recorded['chain']}: inventory keeps {len(recorded['sequence'])} "
            f"residues spanning ids {recorded['residue_ids'][0]}..."
            f"{recorded['residue_ids'][-1]}; the file holds {measured['n_residues']} "
            f"spanning {measured['first_id']}...{measured['last_id']}"
        )
        if len(recorded["sequence"]) != measured["n_residues"]:
            return FAIL, "a target chain lost residues on the way in. " + "; ".join(lines)
        if recorded["sequence"] != measured["sequence"]:
            return FAIL, "a target sequence is not the file's. " + "; ".join(lines)
    return PASS, ("every target residue present in the file is still in the inventory, "
                  "gap and all: " + "; ".join(lines))


def c_inventory_residue_ids_in_order(c, spec):
    inv = c.inventory()
    if not inv:
        return FAIL, "no design_inventory.json was written"
    row = inv["designs"][0]
    target = row["target_chains"][0]
    ids = target["residue_ids"]
    measured = c.probe_chain(target["chain"])
    ascending = all(b > a for a, b in zip(ids, ids[1:]))
    detail = (f"chain {target['chain']} residue_ids start at {ids[0]} and run to "
              f"{ids[-1]}, strictly ascending: {ascending}; probe_design measures first_id="
              f"{measured['first_id']} last_id={measured['last_id']}; the binder's own ids "
              f"start at {row['binder_residue_ids'][0]}")
    if ids[0] != spec["first_id"] or measured["first_id"] != spec["first_id"]:
        return FAIL, f"expected first_id={spec['first_id']}. " + detail
    if not ascending or ids != measured["residue_ids"]:
        return FAIL, "the recorded order is not the file's order. " + detail
    return PASS, detail


def c_rmsd_mapping_asserted(c, spec):
    """The RMSD must be mapped by POSITION, and the mapping must be asserted, not assumed.

    Shown on an input where matching by residue number would be wrong: the design's target
    is numbered from 3 and the predictor renumbers from 1, so any correspondence by number
    is off by two somewhere (ef74646e).
    """
    row = c.design()
    inv = c.inventory()
    if not row or not inv:
        return FAIL, "node 04 produced no scored row"
    pdbs = [n for n in c.files("n4", "inputs") if n.endswith(".pdb")]
    first_pred_id = None
    if pdbs:
        for line in open(c.path("n4", "inputs", pdbs[0])):
            if line.startswith("ATOM"):
                first_pred_id = int(line[22:26])
                break
    designed = inv["designs"][0]
    detail = (
        f"sc_rmsd_A={row['sc_rmsd_A']} over sc_rmsd_n_residues={row['sc_rmsd_n_residues']}; "
        f"the design's binder ids start at {designed['binder_residue_ids'][0]} and its "
        f"target at {designed['target_chains'][0]['residue_ids'][0]}, while the predicted "
        f"PDB's first residue is numbered {first_pred_id} -- so the two numberings differ "
        f"and the mapping cannot have been by residue number"
    )
    if row["sc_rmsd_n_residues"] != designed["n_binder_residues"]:
        return FAIL, ("the RMSD was taken over a different number of residues than the "
                      "design has. " + detail)
    if row["sc_rmsd_A"] is None:
        return FAIL, "no RMSD was reported at all. " + detail
    return PASS, detail


def c_unknown_residues_counted(c, spec):
    inv = c.inventory()
    if not inv:
        return FAIL, "no design_inventory.json was written"
    row = inv["designs"][0]
    measured = c.probe_chain(row["binder_chain"])
    detail = (f"inventory n_unknown={row['n_unknown']} for a binder_sequence of "
              f"{row['binder_sequence']}; probe_design counts {measured['n_unknown']} "
              f"unknown in {measured['sequence']}")
    if row["n_unknown"] < spec["min_unknown"]:
        return FAIL, f"expected at least {spec['min_unknown']}. " + detail
    if row["n_unknown"] != measured["n_unknown"]:
        return FAIL, "the two counts disagree. " + detail
    decl = c.inputs_decl() or {}
    carried = [r.get("n_unknown") for r in decl.get("inputs") or []]
    return PASS, detail + f"; node 02 carries n_unknown={carried} forward"


def c_all_unknown_design_skipped(c, spec):
    """A sequence that is entirely unknown has nothing to fold, so it must be skipped.

    Its own run: `also=allnonstandard` in the harness field, staged under n1b/n2b.
    """
    inv = c.inventory("n1b")
    if not inv:
        return FAIL, ("the second run (also=allnonstandard) produced no inventory; rc="
                      f"{c.rc('n1b')}, log: {c.log('n1b').strip()[-300:]}")
    rows = inv["designs"]
    probe = c.probe("g2")
    detail = (f"the all-unknown run ingested {inv['n_designs']} design(s) and skipped "
              f"{inv.get('n_skipped')}; probe_design measures the generated binder as "
              f"{probe[0]['chains'][0]['n_unknown']}/"
              f"{probe[0]['chains'][0]['n_residues']} unknown")
    skipped = [r for r in rows if r.get("skipped")]
    if not skipped:
        return FAIL, ("nothing was skipped, so an all-X sequence would be sent to the "
                      "predictor. " + detail)
    reason = skipped[0].get("skip_reason") or ""
    decl = c.inputs_decl("n2b") or {}
    written = [n for n in c.files("n2b", "outputs") if n.endswith(".fasta")]
    if written:
        return FAIL, (f"node 02 still wrote {written} for a skipped design. " + detail)
    return PASS, (detail + f"; skip_reason={reason!r}; node 02 declared "
                  f"{decl.get('n_inputs')} input(s) and wrote no FASTA")


# ===========================================================================
# node 02
# ===========================================================================

def _mono_files(c, step="n2"):
    return [n for n in c.files(step, "outputs")
            if n.startswith("mono_") and n.endswith(".fasta")]


def c_mono_input_written(c, spec):
    """A monomer design still has a backbone to recover, so it must be folded, not skipped."""
    files = _mono_files(c)
    decl = c.inputs_decl()
    inv = c.inventory()
    if decl is None or inv is None:
        return FAIL, ("node 02 produced no prediction_inputs.json; rc="
                      f"{c.rc('n2')}, log: {c.log('n2').strip()[-300:]}")
    sizes = {n: len(c.text("n2", "outputs", n) or "") for n in files}
    detail = (f"node 01 ingested {inv['n_designs']} design(s) with sources "
              f"{inv.get('designed_chain_sources')} and skipped {inv.get('n_skipped')}; "
              f"node 02 declared n_inputs={decl.get('n_inputs')} in mode "
              f"{decl.get('mode')!r} and wrote {sizes}")
    if len(files) != spec["n_designs"] or decl.get("n_inputs") != spec["n_designs"]:
        return FAIL, f"expected {spec['n_designs']} mono input(s). " + detail
    if any(v == 0 for v in sizes.values()):
        return FAIL, "an input file is empty. " + detail
    return PASS, detail


def c_cplx_input_not_written(c, spec):
    """And it must NOT be given a complex input: there is no target to score against.

    An absence, asserted. 026's defect is a monomer prediction whose confidence is
    reported as an interface, and the quiet way to arrive there is to build a complex
    input out of a design that has no second chain.
    """
    decl = c.inputs_decl()
    if decl is None:
        return FAIL, f"node 02 produced no prediction_inputs.json; rc={c.rc('n2')}"
    cplx = [n for n in c.files("n2", "outputs") if n.startswith("cplx_")]
    modes = sorted({r.get("mode") for r in decl.get("inputs") or []})
    colons = {n: (c.text("n2", "outputs", n) or "").count(":") for n in _mono_files(c)}
    orders = [r.get("chain_order") for r in decl.get("inputs") or []]
    detail = (f"node 02 wrote {cplx or 'no'} cplx_ file(s); per-record modes {modes}; "
              f"chain_order per record {orders}; colons in the mono files {colons}; "
              f"score_interface={decl.get('score_interface')}, phase={decl.get('phase')}")
    if cplx:
        return FAIL, detail
    if modes != ["mono"] or any("target" in (o or []) for o in orders):
        return FAIL, "a record asks for a target that is not there. " + detail
    if any(colons.values()):
        return FAIL, "a mono file carries a chain break. " + detail
    return PASS, detail


def c_rerun_idempotent_or_refused(c, spec):
    """Running again over a populated directory has to be safe, or refused without damage.

    Two runs, because the two halves are different claims. The second run over node 02's
    own outputs may legitimately redo the work, so either exit 0 with byte-identical files
    or refuse. The third run adds a file from a design this run knows nothing about, which
    is the state 029's defect actually fired in -- a directory left by an EARLIER run
    (4d1a523e) -- and there a refusal is the only safe answer, with the existing file
    untouched.
    """
    again, stale = c.step("n2-again"), c.step("n2-stale")
    if again is None or stale is None:
        return FAIL, "the harness recorded no second or third node 02 run"
    again_note = json.loads(again["note"] or "{}")
    stale_note = json.loads(stale["note"] or "{}")
    detail = (
        f"second run over its own outputs: rc={again['rc']}, files byte-identical "
        f"afterwards: {again_note.get('identical')} ({again_note.get('before')} -> "
        f"{again_note.get('after')}). Third run with a stale mono_design_99.fasta from a "
        f"design this run does not know: rc={stale['rc']}, the stale file survived "
        f"unchanged: {stale_note.get('stale_survived')}. Its message: "
        + " ".join(((stale.get("stderr") or stale.get("stdout") or "")
                    .strip().splitlines() or [""])[-1].split())[:220]
    )
    if again["rc"] == 0 and not again_note.get("identical"):
        return FAIL, ("the second run exited 0 and changed the files, which is the "
                      "non-idempotent case 029 hit. " + detail)
    if again["rc"] != 0 and not again_note.get("identical"):
        return FAIL, "the refusal did not leave the existing files alone. " + detail
    if stale["rc"] == 0:
        return FAIL, ("a stale file from an unknown design was accepted, so a later node "
                      "would predict it. " + detail)
    if not stale_note.get("stale_survived"):
        return FAIL, "the refusal deleted or rewrote the stale file. " + detail
    return PASS, detail


def c_mono_record_count(c, spec):
    files = _mono_files(c)
    if not files:
        return FAIL, "node 02 wrote no mono_*.fasta"
    counts = {}
    for name in files:
        body = c.text("n2", "outputs", name)
        counts[name] = body.count(">")
    bad = {k: v for k, v in counts.items() if v != spec["n_records"]}
    decl = c.inputs_decl() or {}
    verified = [r.get("verified_records_in_file") for r in decl.get("inputs") or []]
    detail = (f"record counts read off disk: {counts}; node 02's own read-back recorded "
              f"verified_records_in_file={verified}")
    return (PASS if not bad else FAIL), (
        detail if not bad else f"expected {spec['n_records']} per file. " + detail
    )


def c_mono_has_no_colon(c, spec):
    files = _mono_files(c)
    if not files:
        return FAIL, "node 02 wrote no mono_*.fasta"
    found = {}
    for name in files:
        body = c.text("n2", "outputs", name)
        found[name] = body.count(":")
    decl = c.inputs_decl() or {}
    verified = [r.get("verified_colons_in_sequence") for r in decl.get("inputs") or []]
    detail = (f"colons counted in the files on disk: {found}; node 02's own read-back "
              f"recorded verified_colons_in_sequence={verified}; the alphabet node 02 "
              f"accepts is {decl.get('alphabet')!r}, which contains no colon")
    return (PASS if not any(found.values()) else FAIL), detail


def _cplx_gap(c, what):
    """No cplx_ file to read. A gap if none was asked for, a FAILURE if one was.

    Node 02 writing no complex when score_interface asked for one is the defect this
    whole section is about, and reporting it as a gap would be the §8.5 case: a failure
    that does not appear as a failure. So the two are told apart here, by what the case
    asked for rather than by what came out.
    """
    decl = c.inputs_decl() or {}
    files = [n for n in c.files("n2", "outputs") if n.startswith("cplx_")]
    detail = (f"node 02 declared mode={decl.get('mode')!r} and wrote {files or 'no'} "
              f"cplx_ file(s), so {what} cannot be measured")
    if c.wants_interface():
        return FAIL, (
            f"{detail}. This case set PARAM_SCORE_INTERFACE=true, so a complex WAS asked "
            f"for and node 02 did not write one. node 02 rc={c.rc('n2')}, node 01b "
            f"rc={c.rc('nalign')}."
        )
    return GAP, detail + ". " + NO_COMPLEX_STAGED


def c_cplx_colon_count(c, spec):
    files = [n for n in c.files("n2", "outputs") if n.startswith("cplx_")]
    if not files:
        return _cplx_gap(c, "the colon count")
    counts = {n: c.text("n2", "outputs", n).count(":") for n in files}
    bad = {k: v for k, v in counts.items() if v != spec["n_colons"]}
    return (PASS if not bad else FAIL), f"colons per cplx_ file: {counts}"


def c_cplx_halves_nonempty(c, spec):
    files = [n for n in c.files("n2", "outputs") if n.startswith("cplx_")]
    if not files:
        return _cplx_gap(c, "whether both halves are non-empty")
    sizes = {}
    for name in files:
        body = "".join(l for l in c.text("n2", "outputs", name).splitlines()
                       if not l.startswith(">"))
        sizes[name] = [len(h) for h in body.split(":")]
    bad = {k: v for k, v in sizes.items() if len(v) != 2 or min(v) == 0}
    return (PASS if not bad else FAIL), f"half lengths per cplx_ file: {sizes}"


def c_mono_file_per_design(c, spec):
    inv = c.inventory()
    files = _mono_files(c)
    if not inv:
        return FAIL, "no design_inventory.json was written"
    designs = [r["design"] for r in inv["designs"] if not r.get("skipped")]
    expected = sorted(f"mono_{d}.fasta" for d in designs)
    detail = (f"{len(designs)} design(s) {designs} produced {len(files)} file(s) {files}; "
              f"one per design would be {expected}")
    if len(designs) != spec["n_designs"]:
        return FAIL, f"expected {spec['n_designs']} designs. " + detail
    if sorted(files) != expected:
        return FAIL, "the files are not one per design, named after it. " + detail
    bodies = {n: c.text("n2", "outputs", n) for n in files}
    if len({v for v in bodies.values()}) != len(bodies):
        return FAIL, "two designs' files have identical contents. " + detail
    return PASS, detail


def c_no_batch_file_written(c, spec):
    outputs = c.files("n2", "outputs")
    multi = {}
    for name in outputs:
        if name.endswith(".fasta"):
            multi[name] = c.text("n2", "outputs", name).count(">")
    detail = (f"every FASTA node 02 wrote, with its record count: {multi}. A batch file "
              f"would hold more than one, and 026 conflates designs precisely because its "
              f"predictor output lands in one undistinguished directory (818954e9)")
    if not multi:
        return FAIL, "node 02 wrote no FASTA at all. " + detail
    if any(v != 1 for v in multi.values()):
        return FAIL, "a file holds more than one record. " + detail
    return PASS, detail


def c_record_ids_filesystem_safe(c, spec):
    decl = c.inputs_decl()
    if not decl:
        return FAIL, "no prediction_inputs.json was written"
    safe = re.compile(r"^[A-Za-z0-9._-]+$")
    ids = [r["record_id"] for r in decl["inputs"]]
    bad = [i for i in ids if not safe.match(i)]
    headers = []
    for row in decl["inputs"]:
        body = c.text("n2", "outputs", row["path"]) or ""
        headers.append(body.splitlines()[0][1:] if body.startswith(">") else None)
    detail = (f"record ids {ids}, all matching ^[A-Za-z0-9._-]+$: {not bad}; the headers "
              f"actually written are {headers}; ColabFold names its output files after "
              f"these, which is why 037 makes them safe (54bf21b6)")
    if bad:
        return FAIL, f"unsafe: {bad}. " + detail
    if headers != ids:
        return FAIL, "the header written is not the declared record id. " + detail
    return PASS, detail


def c_recorded_lengths_sum_asserted(c, spec):
    decl = c.inputs_decl()
    if not decl:
        return FAIL, "no prediction_inputs.json was written"
    lines, bad = [], []
    for row in decl["inputs"]:
        total = sum(row["chain_lengths"])
        body = "".join(l for l in (c.text("n2", "outputs", row["path"]) or "").splitlines()
                       if not l.startswith(">"))
        residues = len(body.replace(":", ""))
        lines.append(f"{row['design']}: chain_lengths={row['chain_lengths']} sum={total} "
                     f"expected_residues={row['expected_residues']} and the file on disk "
                     f"holds {residues}")
        if not (total == row["expected_residues"] == residues):
            bad.append(row["design"])
    detail = "; ".join(lines)
    return (PASS if not bad else FAIL), (
        detail if not bad else f"disagreement for {bad}. " + detail
    )


# ===========================================================================
# node 03
# ===========================================================================

def _argvs(c, mode=None):
    """{'<design>/<mode>': argv}. Keyed by BOTH, because the complex path runs two.

    The first version keyed on design alone. On the complex path node 03 runs each design
    twice -- mono through alphafold2_ptm and cplx through alphafold2_multimer_v3 -- so one
    of the two silently replaced the other in this dict, and whichever survived was the
    one every argv check then read. The model-type check would have compared the monomer's
    argv against the multimer value and failed a correct node, or passed a wrong one,
    depending only on manifest order (§8.16).
    """
    manifest = c.manifest()
    if not manifest:
        return None
    out = {}
    for prediction in manifest["predictions"]:
        this = prediction.get("mode", "mono")
        if mode is not None and this != mode:
            continue
        out[f"{prediction['design']}/{this}"] = prediction.get("argv") or []
    return out


def c_argv_excludes(c, spec):
    argvs = _argvs(c)
    if argvs is None:
        return FAIL, f"no prediction_manifest.json; node 03 rc={c.rc('n3')}"
    flag = spec["flag"]
    offenders = {d: v for d, v in argvs.items() if flag in v}
    manifest = c.manifest()
    detail = (f"{len(argvs)} resolved argv(s) over {sorted(argvs)}, none containing "
              f"{flag}: {not offenders}. First: {argvs[sorted(argvs)[0]]}. It records "
              f"host_url_passed={manifest.get('host_url_passed')} and "
              f"network_expected={manifest.get('network_expected')!r}")
    return (PASS if not offenders else FAIL), (
        detail if not offenders else f"{flag} present for {sorted(offenders)}. " + detail
    )


def c_argv_includes(c, spec):
    """The flag is present with the stated value, in every argv the spec is about.

    `mode` in the spec narrows it to one path. Without that, a flag whose correct value
    DIFFERS by path -- --model-type above all -- could only be asserted by picking one
    argv and hoping, which is what keying on design alone used to do for us.
    """
    mode = spec.get("mode")
    argvs = _argvs(c, mode)
    if argvs is None:
        return FAIL, f"no prediction_manifest.json; node 03 rc={c.rc('n3')}"
    flag, want = spec["flag"], spec["value"]

    if not argvs:
        # No argv of this mode was resolved. Whether that is a gap or a defect depends on
        # whether a complex was ever declared -- and that distinction is the check, so it
        # is read off node 02's declaration rather than assumed either way (§8.3).
        decl = c.inputs_decl() or {}
        modes = sorted({r.get("mode") for r in (decl.get("inputs") or [])})
        if mode and mode not in modes:
            return GAP, (f"node 03 resolved no argv with mode={mode!r}; node 02 declared "
                         f"modes={modes or 'none'}. " + NO_COMPLEX_STAGED)
        return FAIL, (f"node 02 declared modes={modes} including {mode!r}, but node 03 "
                      f"resolved no argv for it, so the prediction was never run")

    got = {key: (argv[argv.index(flag) + 1] if flag in argv else None)
           for key, argv in argvs.items()}
    wrong = {k: v for k, v in got.items() if v != want}
    detail = (f"{flag} per design/mode{'' if not mode else f' (mode={mode})'}: {got}; "
              f"wanted {want!r}")
    return (PASS if not wrong else FAIL), detail


def c_no_network_reached(c, spec):
    manifest = c.manifest()
    if not manifest:
        return FAIL, f"no prediction_manifest.json; node 03 rc={c.rc('n3')}"
    argvs = _argvs(c)
    urls = {d: [a for a in v if "://" in a or "api.colabfold" in a]
            for d, v in argvs.items()}
    offenders = {d: v for d, v in urls.items() if v}
    log = c.log("n3")
    mentions = [l for l in log.splitlines() if "api.colabfold" in l]
    detail = (f"host_url_passed={manifest.get('host_url_passed')}, "
              f"network_expected={manifest.get('network_expected')!r}, no URL-shaped "
              f"argument in any argv: {not offenders}, no api.colabfold.com line in the "
              f"log: {not mentions}. The suite can say this much about what node 03 "
              f"ASKS FOR; that nothing leaves the host is established by the container "
              f"running with --network none, which this harness passes and the run "
              f"records in run.json")
    bad = offenders or mentions or manifest.get("host_url_passed")
    return (PASS if not bad else FAIL), detail


def c_sigterm_ignoring_child_killed(c, spec):
    step = c.step("worker")
    if step is None:
        return FAIL, "the harness recorded no surviving-worker check"
    seen = json.loads(step["note"])
    manifest = c.manifest() or {}
    statuses = [(p["design"], p["status"], (p.get("reason") or "")[:60])
                for p in manifest.get("predictions") or []]
    detail = (f"the stand-in installed TERM-ignoring traps and left a worker at pid "
              f"{seen.get('pid')}; after node 03 returned, /proc/{seen.get('pid')} "
              f"existed: {seen.get('alive')}. node 03 rc={c.rc('n3')}, manifest statuses "
              f"{statuses}")
    if not seen.get("observable", True):
        # Reported as a gap rather than a pass, which is the only honest answer: an
        # unobservable claim that returns PASS is the check passing for the wrong reason,
        # and that is worse than no check at all (§8.4).
        return GAP, (f"not observable in this run: {seen.get('why')} The claim is "
                     f"established by the host-interpreter run of this same case, where "
                     f"the pid is this host's. " + detail)
    if seen.get("pid") is None:
        return FAIL, "the worker never recorded a pid. " + detail
    if seen.get("alive"):
        return FAIL, ("the worker outlived the node, which is 037's measured failure: "
                      "the client dies and the GPU stays held (755d4af2). " + detail)
    return PASS, detail


def c_timeout_message_carries_output(c, spec):
    manifest = c.manifest() or {}
    reasons = [(p["design"], p.get("reason") or "")
               for p in manifest.get("predictions") or []]
    log = c.log("n3")
    body = " ".join(r for _, r in reasons) + log
    # The stand-in's own words. A timeout message that carried none of them would be
    # telling whoever reads it that something timed out and nothing about what it was doing.
    marks = {
        "names the timeout": bool(re.search(r"timed out|timeout", body, re.I)),
        "quotes the predictor's stdout": "STAND-IN colabfold_batch" in body,
        "quotes the worker line": "worker started" in body,
        "names the budget": "5" in " ".join(r for _, r in reasons),
    }
    retained = [n for n in c.files("n3", "outputs") if n.endswith(".log")]
    detail = (f"{marks}; reasons recorded in the manifest: {reasons}; log files retained "
              f"in node 03's outputs: {retained}")
    return (PASS if all(marks.values()) else FAIL), detail


def c_timeout_tail_from_log_file(c, spec):
    """When stdout was buffered away, the reason must come from the predictor's own log.

    Measured in control C6: real colabfold_batch stdout is block-buffered, so a killed
    prediction leaves a ZERO-byte `<design>_predict.log` and the reason used to read
    "Last lines: " with nothing after it. The lines were never lost -- colabfold logs
    through the logging module to its own `<record>_log.txt`, which node 03 has already
    collected by then.

    So three things have to hold, and the middle one is the point: the reason has to carry
    the log file's content, it has to NAME that file as the source rather than implying
    stdout, and node 03's own stdout capture has to be genuinely empty -- otherwise this
    case is testing the ordinary path under a different name.
    """
    manifest = c.manifest() or {}
    rows = manifest.get("predictions") or []
    if not rows:
        return FAIL, f"no predictions in the manifest; node 03 rc={c.rc('n3')}"

    retained = [n for n in c.files("n3", "outputs") if n.endswith("_log.txt")]
    marks, detail = {}, []
    for row in rows:
        reason = row.get("reason") or ""
        source = row.get("tail_source") or ""
        marks[f"{row['design']}: names the timeout"] = bool(
            re.search(r"timed out", reason, re.I))
        # The log file's own words, which the stand-in writes nowhere else.
        marks[f"{row['design']}: carries the log file's lines"] = (
            "STAND-IN logfile" in reason)
        marks[f"{row['design']}: names the fallback source"] = (
            "_log.txt" in source and "buffered" in source)
        marks[f"{row['design']}: does not claim stdout"] = (
            "the predictor's stdout" != source)
        detail.append(f"{row['design']}: tail_source={source[:80]!r}")

    # A stand-in that printed to stdout would make all of the above pass through the
    # ordinary path. The whole case rests on this being empty.
    predict_logs = [n for n in c.files("n3", "outputs") if n.endswith("_predict.log")]
    sizes = {n: c.size("n3", "outputs", n) for n in predict_logs}
    marks["stdout capture really was empty"] = bool(sizes) and all(
        s == 0 for s in sizes.values())

    return (PASS if all(marks.values()) else FAIL), (
        f"{marks}; {'; '.join(detail)}; colabfold logs collected: {retained}; "
        f"stdout capture sizes: {sizes}")


def c_failed_design_in_manifest(c, spec):
    manifest = c.manifest()
    if not manifest:
        return FAIL, f"no prediction_manifest.json; node 03 rc={c.rc('n3')}"
    rows = {p["design"]: p["status"] for p in manifest["predictions"]}
    failed = [d for d, s in rows.items() if s != "ok"]
    reasons = {p["design"]: (p.get("reason") or "")[:90]
               for p in manifest["predictions"] if p["status"] != "ok"}
    detail = (f"manifest declares {manifest['n_declared']}, ok {manifest['n_ok']}, failed "
              f"{manifest['n_failed']}; statuses {rows}; reasons {reasons}")
    if not failed:
        return FAIL, "nothing is recorded as failed, so the failure vanished. " + detail
    if not all(reasons.values()):
        return FAIL, "a failed design carries no reason. " + detail
    return PASS, detail


def c_batch_continues_after_failure(c, spec):
    manifest = c.manifest()
    if not manifest:
        return FAIL, f"no prediction_manifest.json; node 03 rc={c.rc('n3')}"
    declared = len((c.inputs_decl() or {}).get("inputs") or [])
    detail = (f"node 02 declared {declared} design(s); the manifest holds "
              f"{manifest['n_declared']} with {manifest['n_ok']} ok and "
              f"{manifest['n_failed']} failed, and node 03 exited {c.rc('n3')}")
    if manifest["n_declared"] != declared:
        return FAIL, "a declared design is missing from the manifest. " + detail
    if manifest["n_ok"] < 1 or manifest["n_failed"] < 1:
        return FAIL, "this case needs at least one of each to say anything. " + detail
    return PASS, detail


def c_predictor_output_retained(c, spec):
    manifest = c.manifest()
    if not manifest:
        return FAIL, f"no prediction_manifest.json; node 03 rc={c.rc('n3')}"
    failed = [p for p in manifest["predictions"] if p["status"] != "ok"]
    outputs = c.files("n3", "outputs")
    kept = {}
    for prediction in failed:
        log = prediction.get("log")
        kept[prediction["design"]] = {
            "log_recorded": log,
            "log_on_disk": bool(log) and log in outputs,
            "reason_length": len(prediction.get("reason") or ""),
        }
    detail = (f"failed design(s) {[p['design'] for p in failed]}; what was kept for each: "
              f"{kept}; node 03's outputs directory holds {outputs}")
    if not failed:
        return FAIL, "no design failed, so nothing is being retained. " + detail
    if not all(v["log_on_disk"] for v in kept.values()):
        return FAIL, ("a failed design's predictor output was not kept, which is what "
                      "made 029 unable to tell a bad design from an unmodelled target "
                      "(2359ffa4). " + detail)
    return PASS, detail


# ===========================================================================
# node 04
# ===========================================================================

def _truth_row(c, design):
    truth = c.truth() or {}
    for row in truth.get("designs") or []:
        if row["design"] == design:
            return row
    return {}


def c_sc_rmsd(c, spec):
    row = c.design()
    if not row:
        return FAIL, f"node 04 produced no scored row; rc={c.rc('n4')}"
    truth = _truth_row(c, row["design"])
    got, want, tol = row["sc_rmsd_A"], spec["value"], spec["tolerance"]
    detail = (
        f"node 04 reports sc_rmsd_A={got} over {row['sc_rmsd_n_residues']} residues. "
        f"The fixture's constructed value is {want} +/-{tol}; make_predictions.py, which "
        f"runs its own Kabsch that node 04 never sees, recorded "
        f"sc_rmsd_A={truth.get('sc_rmsd_A')} by construction and "
        f"{truth.get('sc_rmsd_A_as_written')} measured off the written PDB. shape="
        f"{(c.truth() or {}).get('shape')!r}"
    )
    if not near(got, want, tol):
        return FAIL, detail
    as_written = truth.get("sc_rmsd_A_as_written")
    if as_written is not None and not near(got, as_written, 1e-3):
        return FAIL, ("node 04 and the independent measurement of the same file "
                      "disagree. " + detail)
    return PASS, detail


def _plddt(c, spec, key):
    row = c.design()
    if not row:
        return FAIL, f"node 04 produced no scored row; rc={c.rc('n4')}"
    truth = _truth_row(c, row["design"])
    manifest = c.manifest() or {}
    order = (manifest.get("predictions") or [{}])[0].get("chain_order")
    lengths = (manifest.get("predictions") or [{}])[0].get("chain_lengths")
    got, want, tol = row[key], spec["value"], spec["tolerance"]
    detail = (f"node 04 reports {key}={got} against a fill value of {truth.get(key)} "
              f"written into the {key.split('_')[-1]} block; the prediction's chain_order "
              f"is {order} with lengths {lengths}, so the block node 04 had to slice is "
              f"{'first' if order and order[0] in key else 'not first'}. Fixture wants "
              f"{want} +/-{tol}")
    if not near(got, want, tol):
        return FAIL, detail
    return PASS, detail


def c_plddt_binder_in_complex(c, spec):
    """The binder's confidence IN the complex is the BINDER's block, not the target's.

    Not a fill-value comparison. This re-slices the raw plddt array itself, out of the
    scores file node 04 read, using the chain_order and chain_lengths the manifest
    records -- arithmetic that shares nothing with the node under test (§8.4). Node 04's
    number has to equal the binder block's mean and not the target block's.

    And it refuses a case whose two blocks are indistinguishable. If the target and binder
    fill values are equal, a node that sliced the wrong block passes this check, and a
    check that cannot tell the two apart is not measuring the thing it is named after.
    That refusal is the reason the case carries two different fills.
    """
    row = c.design(spec.get("design"))
    if not row:
        return _no_row(c, spec)
    manifest = c.manifest()
    if not manifest:
        return FAIL, f"no prediction_manifest.json; node 03 rc={c.rc('n3')}"

    cplx = [p for p in manifest["predictions"] if p.get("mode") == "cplx"]
    if not cplx:
        return GAP, ("no prediction in this run has mode='cplx', so there is no complex "
                     "block to slice. " + NO_COMPLEX_STAGED)
    if len(cplx) > 1 and not spec.get("design"):
        return FAIL, (f"{len(cplx)} complex predictions and the fixture names no design, "
                      "so which row this is about is undecided")
    prediction = next((p for p in cplx if p["design"] == row["design"]), None)
    if prediction is None:
        return FAIL, (f"node 04 scored {row['design']!r} but no complex prediction was "
                      f"recorded for it; complexes are for "
                      f"{[p['design'] for p in cplx]}")

    name = prediction.get("scores_rank_001")
    scores = None
    for where in (("n4", "inputs"), ("n3", "outputs"), ("pred",)):
        if name:
            scores = c.json(*where, name)
        if scores is not None:
            break
    if scores is None or "plddt" not in scores:
        return FAIL, (f"the complex scores file {name!r} is not readable or carries no "
                      f"'plddt' key, so the slice cannot be checked independently")

    order, lengths = prediction["chain_order"], prediction["chain_lengths"]
    array = scores["plddt"]
    if sum(lengths) != len(array):
        return FAIL, (f"the recorded lengths {lengths} sum to {sum(lengths)} but the "
                      f"plddt array holds {len(array)}, so every slice is of the wrong "
                      "residues and each one still produces a plausible number")
    blocks, at = {}, 0
    for chain, length in zip(order, lengths):
        blocks[chain] = array[at:at + length]
        at += length
    means = {k: round(sum(v) / len(v), 2) for k, v in blocks.items() if v}

    got = row.get("plddt_binder_in_complex")
    detail = (f"node 04 reports plddt_binder_in_complex={got}; re-slicing the array by "
              f"chain_order={order} lengths={lengths} gives {means}")
    if "binder" not in means or "target" not in means:
        return FAIL, f"the complex has no binder and target block to tell apart. {detail}"
    if abs(means["binder"] - means["target"]) < 1.0:
        return FAIL, (
            "the binder and target blocks are within 1.0 pLDDT of each other, so a node "
            "that sliced the wrong one would pass this check. The case has to write two "
            f"clearly different fills for it to measure anything (§8.4). {detail}"
        )
    if got is None:
        return FAIL, ("the metric is absent on a complex that has a binder block. "
                      f"null_reasons={row.get('null_reasons')}. {detail}")
    if not near(got, means["binder"], spec.get("tolerance", 0.05)):
        if near(got, means["target"], spec.get("tolerance", 0.05)):
            return FAIL, ("node 04 reported the TARGET's block as the binder's "
                          "confidence, which is 026's defect by a different route. "
                          + detail)
        return FAIL, "the value matches neither block. " + detail
    return PASS, detail


def c_plddt_binder(c, spec):
    return _plddt(c, spec, "plddt_binder")


def c_plddt_target(c, spec):
    return _plddt(c, spec, "plddt_target")


def _row_for(c, spec):
    """The row this expectation is about.

    `design` in the spec where it matters. A case with one good design and one that lost
    its scores file has two rows, and reading the first would ask the question of the
    design that worked -- which is a check that passes while measuring the wrong thing.
    """
    return c.design(spec.get("design"))


def _no_row(c, spec):
    which = spec.get("design") or "the first design"
    rows = [r["design"] for r in (c.scores() or {}).get("designs") or []]
    return FAIL, (f"node 04 wrote no row for {which}; it wrote {rows}; "
                  f"rc={c.rc('n4')}")


def c_metric_is_null(c, spec):
    row = _row_for(c, spec)
    if not row:
        return _no_row(c, spec)
    metric = spec["metric"]
    if metric not in row:
        return FAIL, f"{metric} is not a column in design_scores.json at all"
    scored = c.scores()
    reasons = (scored.get("null_reasons") or {}).get(row["design"], {})
    detail = (f"{row['design']}: {metric}={row[metric]!r}, status={row['status']!r}; "
              f"null_reasons for it name {sorted(reasons)}")
    return (PASS if row[metric] is None else FAIL), detail


def c_metric_not_value(c, spec):
    row = _row_for(c, spec)
    if not row:
        return _no_row(c, spec)
    metric, forbidden = spec["metric"], spec["value"]
    got = row[metric]
    if forbidden == "inf":
        is_forbidden = got is not None and isinstance(got, float) and math.isinf(got)
        shown = "inf"
    else:
        is_forbidden = got == forbidden
        shown = repr(forbidden)
    csv = c.metrics_csv() or ""
    detail = (f"{metric}={got!r}, which is not {shown}. 029 writes 0.0 for a missing iptm "
              f"and inf for a missing PAE, and both read as measurements (§8.7). The CSV "
              f"renders it as an empty cell: {'yes' if ',,' in csv or csv.rstrip().endswith(',') else 'check the row'}")
    return (PASS if not is_forbidden else FAIL), detail


def c_null_reason_recorded(c, spec):
    row = _row_for(c, spec)
    scored = c.scores()
    if not row or not scored:
        return _no_row(c, spec)
    metric = spec["metric"]
    reasons = (scored.get("null_reasons") or {}).get(row["design"], {})
    reason = reasons.get(metric)
    detail = (f"null_reasons[{row['design']}][{metric}] = {reason!r}")
    if not reason:
        return FAIL, ("a null with no recorded reason is indistinguishable from a null "
                      "nobody noticed. " + detail)
    if len(reason) < 40:
        return FAIL, "the reason is too short to say anything. " + detail
    return PASS, detail


def c_pae_interaction(c, spec):
    row = c.design()
    if not row:
        return FAIL, f"node 04 produced no scored row; rc={c.rc('n4')}"
    truth = _truth_row(c, row["design"])
    options = (c.truth() or {}).get("options") or {}
    got, want, tol = row["pae_interaction"], spec["value"], spec["tolerance"]
    detail = (f"pae_interaction={got} from blocks set to "
              f"target->binder={options.get('pae_target_binder')} and "
              f"binder->target={options.get('pae_binder_target')}; the mean over BOTH is "
              f"{truth.get('pae_interaction')}, and averaging one block alone would give "
              f"{options.get('pae_target_binder')} or "
              f"{options.get('pae_binder_target')}. Fixture wants {want} +/-{tol}")
    if got is None:
        return FAIL, detail
    return (PASS if near(got, want, tol) else FAIL), detail


def c_metric_does_not_gate(c, spec):
    row = c.design()
    scored = c.scores()
    if not row or not scored:
        return FAIL, f"node 04 produced no scored row; rc={c.rc('n4')}"
    metric = spec["metric"]
    gating = scored.get("gating_metrics") or []
    diagnostic = scored.get("diagnostic_metrics") or []
    failed_on = row.get("failed_thresholds") or ""
    detail = (f"{metric}={row[metric]} with thresholds {scored.get('thresholds')}; "
              f"gating_metrics={gating}; diagnostic_metrics={diagnostic}; this design "
              f"passes={row['passes']} and failed_thresholds={failed_on!r}; there is no "
              f"{metric} threshold in the output: "
              f"{all(metric not in k for k in (scored.get('thresholds') or {}))}")
    if metric in gating:
        return FAIL, f"{metric} is listed as gating. " + detail
    if metric in failed_on:
        return FAIL, f"{metric} appears in failed_thresholds. " + detail
    if metric not in diagnostic:
        return FAIL, f"{metric} is in neither list, so nothing says what it does. " + detail
    if not row["passes"]:
        return FAIL, ("the design does not pass, so this case cannot show that a low "
                      f"{metric} did not cost it anything. " + detail)
    return PASS, detail


def c_low_target_plddt_warned(c, spec):
    row = c.design()
    if not row:
        return FAIL, f"node 04 produced no scored row; rc={c.rc('n4')}"
    log = c.log("n4")
    value = row.get("plddt_target")
    warned = [l.strip() for l in log.splitlines()
              if "WARNING" in l and ("plddt" in l.lower() or "target" in l.lower())]
    names_value = any(str(value).split(".")[0] in l for l in warned)
    detail = (f"plddt_target={value}; node 04's warning lines: {warned}; the warning names "
              f"the value: {names_value}")
    if value is None:
        return GAP, ("there is no target in this prediction, so there is no target "
                     "confidence to warn about. " + NO_COMPLEX_STAGED)
    if not warned:
        return FAIL, ("029's wrong diagnosis came from this number being absent, not from "
                      "it being low (2359ffa4), so a low value has to be said out "
                      "loud. " + detail)
    return PASS, detail


def c_metrics_merged_by_name(c, spec):
    scored = c.scores()
    truth = c.truth()
    if not scored or not truth:
        return FAIL, f"node 04 produced no design_scores.json; rc={c.rc('n4')}"
    rows = {r["design"]: r for r in scored["designs"]}
    expected = {r["design"]: r for r in truth["designs"]}
    manifest = c.manifest() or {}
    order_in_manifest = [p["design"] for p in manifest.get("predictions") or []]
    mutation = next((s["note"] for s in c.record["steps"]
                     if s["step"].startswith("mutate:")), None)
    lines, bad = [], []
    for design in sorted(rows, key=natural_key):
        got = rows[design].get("plddt_binder")
        want = expected.get(design, {}).get("plddt_binder")
        lines.append(f"{design}: plddt_binder={got} vs {want} written for it")
        if got is None or want is None or abs(got - want) > 1e-6:
            bad.append(design)
    detail = (f"{len(rows)} design(s); the manifest's records are in the order "
              f"{order_in_manifest} after the harness mutation, while node 04 reports "
              + "; ".join(lines)
              + f". Mutation: {mutation}")
    if len(rows) != spec["n_designs"]:
        return FAIL, f"expected {spec['n_designs']} designs. " + detail
    if len({expected[d].get("plddt_binder") for d in expected}) < len(expected):
        return FAIL, ("two designs were given the same confidence, so a positional merge "
                      "would pass this check by accident. " + detail)
    if bad:
        return FAIL, f"wrong numbers for {bad}. " + detail
    for design in rows:
        inv_row = next((r for r in (c.inventory() or {})["designs"]
                        if r["design"] == design), {})
        if rows[design].get("binder_sequence") != inv_row.get("binder_sequence"):
            return FAIL, f"{design}'s sequence came from another design. " + detail
    return PASS, detail


def c_design_recorded_failed(c, spec):
    scored = c.scores()
    if not scored:
        return FAIL, f"node 04 produced no design_scores.json; rc={c.rc('n4')}"
    failed = [r for r in scored["designs"] if r["status"] != "scored"]
    wanted = spec["reason"]
    reasons = {r["design"]: (r.get("reason") or "")[:120] for r in failed}
    detail = (f"n_scored={scored['n_scored']} n_failed={scored['n_failed']}; failed rows "
              f"{reasons}; the fixture expects a reason naming {wanted!r}")
    if not failed:
        return FAIL, "nothing is recorded as failed. " + detail
    hit = [d for d, r in reasons.items()
           if wanted in r or wanted.replace("_", " ") in r]
    if not hit:
        return FAIL, f"no failed row's reason names {wanted!r}. " + detail
    return PASS, detail


def c_design_not_passing(c, spec):
    scored = c.scores()
    if not scored:
        return FAIL, f"node 04 produced no design_scores.json; rc={c.rc('n4')}"
    rows = {r["design"]: (r["passes"], r.get("failed_thresholds"), r.get("reason"))
            for r in scored["designs"]}
    not_passing = [d for d, (p, _, _) in rows.items() if not p]
    fasta = c.fasta() or ""
    exported = [l[1:].split()[0] for l in fasta.splitlines() if l.startswith(">")]
    detail = (f"per design (passes, failed_thresholds, reason): {rows}; n_passing="
              f"{scored['n_passing']}; node 05 exported {exported}")
    if not not_passing:
        return FAIL, ("every design passes, so a missing measurement cost nothing. "
                      + detail)
    leaked = [d for d in not_passing if d in exported]
    if leaked:
        return FAIL, f"a non-passing design reached the FASTA: {leaked}. " + detail
    return PASS, detail


# ===========================================================================
# node 05
# ===========================================================================

def _row_order(html):
    """Design names in the order the static table lists them."""
    body = re.search(r"<table id=\"metrics\".*?</table>", html or "", re.S)
    scope = body.group(0) if body else (html or "")
    return re.findall(r"<td[^>]*>(design_\d+)</td>", scope)


def c_report_row_present(c, spec):
    html = c.static_html()
    if html is None:
        return FAIL, "no report.html was written"
    design = spec["design"]
    order = _row_order(html)
    scored = c.scores() or {}
    row = c.design(design)
    detail = (f"the static table lists {order}; {design} is in it: {design in order}. "
              f"node 04 recorded it with status={(row or {}).get('status')!r}; the report "
              f"covers {len(order)} of {len(scored.get('designs') or [])} scored designs")
    if design not in order:
        return FAIL, ("a design with no prediction must reach the report as a failed row "
                      "rather than vanishing from it (§8.5). " + detail)
    if len(order) != len(scored.get("designs") or []):
        return FAIL, "the report drops designs. " + detail
    return PASS, detail


def c_report_row_states_reason(c, spec):
    html = c.static_html()
    if html is None:
        return FAIL, "no report.html was written"
    design = spec["design"]
    row = c.design(design) or {}
    reason = row.get("reason") or ""
    cell = re.search(
        r"<tr[^>]*>(?:(?!</tr>).)*?" + re.escape(design) + r"(?:(?!</tr>).)*?</tr>",
        html, re.S,
    )
    text = re.sub(r"<[^>]+>", " ", cell.group(0)) if cell else ""
    words = [w for w in re.findall(r"[a-z]{5,}", reason.lower())][:6]
    present = [w for w in words if w in text.lower()]
    detail = (f"node 04's reason for {design} is {reason[:140]!r}; its row in the static "
              f"table reads {' '.join(text.split())[:200]!r}; distinctive words from the "
              f"reason found in the row: {present} of {words}")
    if not reason:
        return FAIL, f"node 04 recorded no reason for {design}. " + detail
    if not cell:
        return FAIL, f"no table row for {design}. " + detail
    if len(present) < max(1, len(words) // 2):
        return FAIL, "the row does not carry the reason. " + detail
    return PASS, detail


def c_report_columns_from_data(c, spec):
    html = c.static_html()
    csv = c.metrics_csv()
    if html is None or csv is None:
        return FAIL, "node 05 wrote no report.html or no validation_metrics.csv"
    row = c.design() or {}
    header = csv.splitlines()[0].split(",")
    found = {}
    for key in spec["columns"]:
        value = row.get(key)
        rendered = value is not None and str(value) in html
        found[key] = {
            "in design_scores": key in row,
            "value": value,
            "value_rendered_in_static_html": rendered,
        }
    labels_present = [h for h in header if h]
    detail = (f"validation_metrics.csv header is {labels_present}; per requested metric: "
              f"{found}")
    missing = [k for k, v in found.items() if not v["in design_scores"]]
    if missing:
        return FAIL, f"{missing} is not a column at all. " + detail
    unrendered = [k for k, v in found.items()
                  if v["value"] is not None and not v["value_rendered_in_static_html"]]
    if unrendered:
        return FAIL, (f"{unrendered} has a value that the static page never shows, so the "
                      "number is not a column a reader can see. " + detail)
    nulls = [k for k, v in found.items() if v["value"] is None]
    if len(nulls) == len(spec["columns"]):
        return GAP, ("every requested metric is null on this run, so nothing was "
                     "rendered to check. " + detail)
    return PASS, detail


def c_report_shows_chain_source(c, spec):
    html = c.static_html()
    if html is None:
        return FAIL, "no report.html was written"
    row = c.design() or {}
    got = row.get("designed_chain_source")
    want = spec["source"]
    detail = (f"node 04 recorded designed_chain_source={got!r}; the fixture expects "
              f"{want!r}; it appears in the static report: {want in html}; the report also "
              f"shows the column label 'Chain resolved by': "
              f"{'Chain resolved by' in html}")
    if got != want:
        return FAIL, "the run did not produce the source this case is about. " + detail
    if want not in html:
        return FAIL, detail
    return PASS, detail


def c_design_rows_in_numeric_order(c, spec):
    html = c.static_html()
    csv = c.metrics_csv()
    if html is None or csv is None:
        return FAIL, "node 05 wrote no report.html or no validation_metrics.csv"
    order = _row_order(html)
    csv_order = [l.split(",")[1] for l in csv.splitlines()[1:] if "," in l]
    natural = sorted(order, key=natural_key)
    lexical = sorted(order)
    detail = (f"the static table lists {order}; numeric order is {natural}; lexical order "
              f"would be {lexical}; validation_metrics.csv lists {csv_order}; node 04 "
              f"declares its order as {(c.scores() or {}).get('order')!r}")
    if len(order) < 11:
        return FAIL, ("fewer than 11 designs, so design_10 never exists and lexical and "
                      "numeric order cannot differ. " + detail)
    if natural == lexical:
        return FAIL, "the two orders coincide here, so this proves nothing. " + detail
    if order != natural:
        return FAIL, detail
    if csv_order != natural:
        return FAIL, "the CSV is in a different order from the report. " + detail
    return PASS, detail


def c_report_renders_provenance_offline(c, spec):
    html = c.static_html()
    full = c.report()
    if html is None:
        return FAIL, "no report.html was written"
    needed = {
        "the uncalibrated-threshold claim": "uncalibrated",
        "the no-binding-evidence claim": "NOTHING about binding",
        "the unscreened notice": "UNSCREENED",
        "the chain-resolution column": "Chain resolved by",
        "the predictor's model type": (c.scores() or {}).get("predictor", {})
                                      .get("model_type", "alphafold2_ptm"),
        "the msa mode": "single_sequence",
    }
    missing = {k: v for k, v in needed.items() if v not in html}
    only_in_script = {k: v for k, v in missing.items() if v in (full or "")}
    detail = (f"with every <script> block removed, report.html is "
              f"{len(html)} of {len(full)} bytes and still carries "
              f"{len(needed) - len(missing)} of {len(needed)} required statements. "
              f"Missing: {sorted(missing)}. Present only inside <script>: "
              f"{sorted(only_in_script)}. External references: "
              f"{sorted(set(re.findall(r'https?://[a-z0-9./-]+', full or '')))}")
    return (PASS if not missing else FAIL), detail


def c_report_renders_table_offline(c, spec):
    html = c.static_html()
    if html is None:
        return FAIL, "no report.html was written"
    order = _row_order(html)
    scored = c.scores() or {}
    designs = [r["design"] for r in scored.get("designs") or []]
    values = []
    for row in scored.get("designs") or []:
        if row.get("sc_rmsd_A") is not None:
            values.append((row["design"], str(row["sc_rmsd_A"]) in html))
    detail = (f"the static table holds rows for {order}, against {designs} scored; each "
              f"design's sc_rmsd_A rendered server-side: {values}; the page has "
              f"{html.count('<tr')} <tr> elements outside any script")
    if sorted(order) != sorted(designs):
        return FAIL, "the static table is not the full set of designs. " + detail
    if any(not rendered for _, rendered in values):
        return FAIL, "a number is missing from the static page. " + detail
    return PASS, detail


def c_passing_fasta_tagged(c, spec):
    fasta = c.fasta()
    if fasta is None:
        return FAIL, "no passing_sequences.fasta was written"
    headers = [l for l in fasta.splitlines() if l.startswith(">")]
    tagged = [h for h in headers if "UNSCREENED" in h]
    scored = c.scores() or {}
    detail = (f"{len(headers)} record(s), {len(tagged)} carrying the UNSCREENED tag; node "
              f"04 counted n_passing={scored.get('n_passing')}; the first header is "
              f"{headers[0][:200] if headers else None!r}")
    if len(headers) != spec["n_records"]:
        return FAIL, f"expected {spec['n_records']} record(s). " + detail
    if len(tagged) != len(headers):
        return FAIL, "a record is untagged. " + detail
    return PASS, detail


def c_tag_has_no_colon(c, spec):
    fasta = c.fasta()
    if fasta is None:
        return FAIL, "no passing_sequences.fasta was written"
    headers = [l for l in fasta.splitlines() if l.startswith(">")]
    tags = [h[h.index("UNSCREENED"):] for h in headers if "UNSCREENED" in h]
    offenders = [t for t in tags if ":" in t]
    detail = (f"{len(tags)} tag(s); none contains a colon: {not offenders}. The tag as "
              f"written is {tags[0][:160]!r}. 026 reads a colon in a record id or "
              f"description as a chain-break marker, so a colon here would turn one "
              f"sequence into two chains at the next predictor (7709809b)")
    if not tags:
        return FAIL, "no tag to check. " + detail
    return (PASS if not offenders else FAIL), detail


def c_fasta_has_no_comment_banner(c, spec):
    fasta = c.fasta()
    if fasta is None:
        return FAIL, "no passing_sequences.fasta was written"
    lines = fasta.splitlines()
    comments = [l for l in lines if l.startswith(";")]
    detail = (f"{len(lines)} line(s), of which {len(comments)} start with ';'. Biopython "
              f"1.85 -- the version in the image -- warns rather than raises on a FASTA "
              f"comment line, so a banner here would travel silently into whatever reads "
              f"this next. The first two lines are {lines[:2]}")
    return (PASS if not comments else FAIL), detail


def c_thresholds_in_output(c, spec):
    scored = c.scores()
    fasta = c.fasta()
    csv = c.metrics_csv()
    if not scored:
        return FAIL, f"node 04 produced no design_scores.json; rc={c.rc('n4')}"
    thresholds = scored.get("thresholds") or {}
    row = c.design() or {}
    html = c.static_html() or ""
    lines, bad = [], []
    for metric in spec["metrics"]:
        key = {"sc_rmsd_A": "max_sc_rmsd_a", "plddt_binder": "min_plddt_binder"}[metric]
        limit = thresholds.get(key)
        beside_in_row = key in row
        in_fasta = bool(fasta) and str(limit) in (fasta or "")
        in_html = str(limit) in html
        lines.append(f"{metric}={row.get(metric)} with {key}={limit}: carried on the row "
                     f"itself: {beside_in_row}, in the FASTA header: {in_fasta}, on the "
                     f"static page: {in_html}")
        if limit is None or not beside_in_row or not in_html:
            bad.append(metric)
    detail = ("; ".join(lines) + f"; thresholds_are_calibrated="
              f"{scored.get('thresholds_are_calibrated')}, which is the claim the number "
              f"is allowed to support")
    if scored.get("thresholds_are_calibrated") is not False:
        return FAIL, ("the output does not say the thresholds are uncalibrated. " + detail)
    return (PASS if not bad else FAIL), (
        detail if not bad else f"{bad} is reported without its threshold. " + detail
    )


# ===========================================================================
# the table
# ===========================================================================

CHECKS = {
    # name: (function, what would make it fail)
    "inventory_designed_chain": (
        c_inventory_designed_chain,
        "probe_design measures the chain independently, so naming the wrong chain, "
        "reporting the wrong length, or copying a sequence that is not the file's all "
        "fail. The swap case makes it fail for a node that takes the first chain.",
    ),
    "inventory_chain_source": (
        c_inventory_chain_source,
        "four cases ask for four different values -- manifest-target-chains, "
        "assumed-first-chain, single-chain -- so a node that hardcoded any one of them "
        "fails the others.",
    ),
    "fallback_stated_in_log": (
        c_fallback_stated_in_log,
        "requires the source, the filename and a word of warning in node 01's own log; "
        "037's fallback was silent for weeks and right every time, which is how nobody "
        "noticed it.",
    ),
    "source_reaches_report": (
        c_source_reaches_report,
        "looks in the STATIC report and in the CSV, so a value that only reaches the "
        "manifest fails.",
    ),
    "inventory_target_chains": (
        c_inventory_target_chains,
        "compares against the chains probe_design finds in the file, so silently "
        "reducing a two-chain target to its first chain fails.",
    ),
    "inventory_chain_order_recorded": (
        c_inventory_chain_order_recorded,
        "requires chain_order and chain_lengths to be the same length for every input; a "
        "node that re-derived the order downstream instead of recording it leaves this "
        "empty.",
    ),
    "inventory_records_gaps": (
        c_inventory_records_gaps,
        "the count must match probe_design's independent reading of the numbering, so "
        "both missing the gap and inventing one fail.",
    ),
    "fusion_disclosed": (
        c_fusion_disclosed,
        "Phase 2. Every numbering gap the inventory measured must be named, with its "
        "boundary residues, in the fusion the declaration discloses. GAP only when no "
        "complex was declared; a `fusion` key carrying an empty list fails.",
    ),
    "target_not_trimmed": (
        c_target_not_trimmed,
        "compares recorded residues against the ones in the file per chain, so trimming "
        "around the gap -- which was argued for and then measured to be much worse "
        "(8e4772a7) -- fails.",
    ),
    "inventory_residue_ids_in_order": (
        c_inventory_residue_ids_in_order,
        "the input is numbered from 3, so a node that renumbered from 1 fails, and so "
        "does one that sorts the ids into a different order from the file's.",
    ),
    "rmsd_mapping_asserted": (
        c_rmsd_mapping_asserted,
        "the design's target starts at residue 3 and the prediction's at 1, so a node "
        "that mapped by residue number produces the wrong count and fails.",
    ),
    "unknown_residues_counted": (
        c_unknown_residues_counted,
        "must equal probe_design's own count of X in the sequence.",
    ),
    "all_unknown_design_skipped": (
        c_all_unknown_design_skipped,
        "runs a second, all-unknown shape and requires both a skip and no FASTA for it; "
        "a node that folded the X-only sequence fails on the written file.",
    ),
    "mono_input_written": (
        c_mono_input_written,
        "a monomer that was skipped rather than folded leaves no file and fails; so does "
        "an empty one, and so does node 02 declaring a different count from what it wrote.",
    ),
    "cplx_input_not_written": (
        c_cplx_input_not_written,
        "an asserted absence: any cplx_ file, any record asking for a target that is not "
        "in the design, or any colon in a mono file fails. This is the quiet route into "
        "026's defect, where a monomer's confidence is reported as an interface.",
    ),
    "rerun_idempotent_or_refused": (
        c_rerun_idempotent_or_refused,
        "three runs of node 02 over one directory. Exiting 0 while changing the files "
        "fails; accepting a stale file from a design this run does not know fails; and so "
        "does a refusal that deletes or rewrites what was already there.",
    ),
    "mono_record_count": (
        c_mono_record_count,
        "counts '>' in the bytes on disk, not in what node 02 says it wrote.",
    ),
    "mono_has_no_colon": (
        c_mono_has_no_colon,
        "counts colons in the bytes on disk; one colon would make the predictor treat one "
        "design as two chains.",
    ),
    "cplx_colon_count": (
        c_cplx_colon_count,
        "Phase 2. Exactly one colon per cplx_ file: zero folds a monomer and reports its "
        "confidence as an interface, two predicts three chains. A case that asked for the "
        "complex path and got no cplx_ file FAILS; one that never asked for it is a gap.",
    ),
    "cplx_halves_nonempty": (
        c_cplx_halves_nonempty,
        "Phase 2. Both halves non-empty, because 'TARGET:' has exactly one colon and "
        "folds one chain, so counting separators is not enough on its own.",
    ),
    "mono_file_per_design": (
        c_mono_file_per_design,
        "requires the file set to be exactly one per design, named after it, with no two "
        "files identical -- the conflation bug writes one file and reports it as several "
        "(818954e9).",
    ),
    "no_batch_file_written": (
        c_no_batch_file_written,
        "any FASTA holding more than one record fails.",
    ),
    "record_ids_filesystem_safe": (
        c_record_ids_filesystem_safe,
        "also requires the header actually written to BE the declared record id, so a "
        "node that sanitised the id but wrote the raw name fails.",
    ),
    "recorded_lengths_sum_asserted": (
        c_recorded_lengths_sum_asserted,
        "three numbers have to agree -- the sum of chain_lengths, expected_residues, and "
        "the residues in the file -- so any one of them being stale fails.",
    ),
    "argv_excludes": (
        c_argv_excludes,
        "reads every resolved argv out of the manifest; 026 passes --host-url "
        "unconditionally, which is exactly what this would catch.",
    ),
    "argv_includes": (
        c_argv_includes,
        "the mono case wants alphafold2_ptm and the complex case "
        "alphafold2_multimer_v3, so a hardcoded model type fails one of them.",
    ),
    "no_network_reached": (
        c_no_network_reached,
        "fails on any URL-shaped argument, any api.colabfold.com line in the log, or "
        "host_url_passed being true. It does NOT by itself prove no packet left the host; "
        "the detail says so and names what does.",
    ),
    "sigterm_ignoring_child_killed": (
        c_sigterm_ignoring_child_killed,
        "the stand-in and the worker it leaves behind both install `trap '' TERM` in one\n"
        "process group, so a node relying on subprocess's own timeout -- which signals the\n"
        "direct child only -- leaves the worker running and fails. A worker that called\n"
        "setsid would survive any process-group kill; that is a stated limit of the fix,\n"
        "not something this check pretends to cover.",
    ),
    "timeout_message_carries_output": (
        c_timeout_message_carries_output,
        "requires the predictor's own stdout lines in the recorded reason, so a bare "
        "'timed out' fails.",
    ),
    "timeout_tail_from_log_file": (
        c_timeout_tail_from_log_file,
        "requires the reason to fall back to the predictor's own log file when stdout "
        "was buffered away, and to name that file as the source. A node that reported an "
        "empty tail, or that claimed stdout, fails.",
    ),
    "failed_design_in_manifest": (
        c_failed_design_in_manifest,
        "a failure that aborted the batch or vanished from the manifest fails, and so "
        "does one recorded without a reason.",
    ),
    "batch_continues_after_failure": (
        c_batch_continues_after_failure,
        "needs at least one ok and one failed in the same manifest, and the manifest to "
        "hold every design node 02 declared.",
    ),
    "predictor_output_retained": (
        c_predictor_output_retained,
        "the failed design's log has to be named in the manifest AND present on disk; 029 "
        "discarded it and then could not diagnose a run (2359ffa4).",
    ),
    "sc_rmsd": (
        c_sc_rmsd,
        "three cases share this check with constructed answers 0, 0 and 2.0: identity, a "
        "rigid rotation plus a 260 A translation, and a rigid-motion-orthogonal "
        "displacement. A node that skipped superposition passes the first and fails the "
        "second by hundreds of angstroms. It also has to agree with make_predictions' own "
        "Kabsch, which it never sees.",
    ),
    "plddt_binder": (
        c_plddt_binder,
        "the binder block is written first in one case and last in another with different "
        "fill values, so a node that hardcoded plddt[:binder_length] passes one and fails "
        "the other.",
    ),
    "plddt_target": (
        c_plddt_target,
        "same slice from the other end, and a third case sets it to 38.4 to show it is "
        "read rather than assumed.",
    ),
    "plddt_binder_in_complex": (
        c_plddt_binder_in_complex,
        "Phase 2. Re-slices the raw plddt array by the recorded chain_order and "
        "chain_lengths and requires node 04's value to be the BINDER block's mean, not "
        "the target's. Refuses a case whose two blocks are within 1.0 of each other, "
        "because a swapped slice would pass that one.",
    ),
    "metric_is_null": (
        c_metric_is_null,
        "the key is genuinely absent from the scores file, so a node with a .get(key, 0.0) "
        "default reports a number and fails.",
    ),
    "metric_not_value": (
        c_metric_not_value,
        "names the specific wrong value each defect produces -- 0.0 for iptm, inf for PAE "
        "-- so reproducing 029's defaults fails.",
    ),
    "null_reason_recorded": (
        c_null_reason_recorded,
        "requires a recorded reason of at least 40 characters for that exact metric; a "
        "null nobody explained fails.",
    ),
    "pae_interaction": (
        c_pae_interaction,
        "the two off-diagonal blocks are 6.0 and 14.0, so averaging one block alone gives "
        "6.0 or 14.0 and fails against 10.0.",
    ),
    "metric_does_not_gate": (
        c_metric_does_not_gate,
        "the value is 38.4, far below any plausible threshold, and the design still has "
        "to pass; a node that gated on it fails.",
    ),
    "low_target_plddt_warned": (
        c_low_target_plddt_warned,
        "requires a WARNING line naming the value. GAP on a monomer run, where there is "
        "no target confidence at all.",
    ),
    "metrics_merged_by_name": (
        c_metrics_merged_by_name,
        "the harness reverses the manifest's records and --plddt-ramp gives every design "
        "a different confidence, so a positional merge gives every design another "
        "design's numbers. The check also refuses to pass if the ramp was left off, "
        "because then it could not fail.",
    ),
    "design_recorded_failed": (
        c_design_recorded_failed,
        "requires a failed row whose reason names the cause; a node that dropped the "
        "design or marked it scored with a null fails.",
    ),
    "design_not_passing": (
        c_design_not_passing,
        "requires at least one design not passing AND that none of those reached the "
        "exported FASTA.",
    ),
    "report_row_present": (
        c_report_row_present,
        "reads the static table, and also requires the report to cover every scored "
        "design, so dropping the failed one fails.",
    ),
    "report_row_states_reason": (
        c_report_row_states_reason,
        "looks for distinctive words from node 04's own reason inside that design's row, "
        "so a blank cell or a generic 'failed' fails.",
    ),
    "report_columns_from_data": (
        c_report_columns_from_data,
        "each metric's actual value has to appear in the static page; a report that "
        "described the metrics in prose instead of tabulating them fails (179003b9).",
    ),
    "report_shows_chain_source": (
        c_report_shows_chain_source,
        "the case withholds design_metrics.json so the source really is the fallback, and "
        "the fallback string has to be on the static page.",
    ),
    "design_rows_in_numeric_order": (
        c_design_rows_in_numeric_order,
        "12 designs, so design_10 exists and lexical order differs from numeric; the "
        "check refuses to pass when the two coincide, because then it proves nothing.",
    ),
    "report_renders_provenance_offline": (
        c_report_renders_provenance_offline,
        "every <script> block is stripped before looking, so provenance that lived in a "
        "JS renderer fails -- which is what 037's report did (cb418341).",
    ),
    "report_renders_table_offline": (
        c_report_renders_table_offline,
        "the table and every design's sc_rmsd_A have to survive stripping the scripts.",
    ),
    "passing_fasta_tagged": (
        c_passing_fasta_tagged,
        "the record count has to be right and every header tagged; an untagged export "
        "fails.",
    ),
    "tag_has_no_colon": (
        c_tag_has_no_colon,
        "a colon anywhere in the tag fails, because the next predictor reads it as a "
        "chain break (7709809b).",
    ),
    "fasta_has_no_comment_banner": (
        c_fasta_has_no_comment_banner,
        "any line starting with ';' fails; biopython 1.85 only warns, so nothing "
        "downstream would.",
    ),
    "thresholds_in_output": (
        c_thresholds_in_output,
        "each gating metric has to carry its own threshold on its own row and on the "
        "static page, and the output has to state that the thresholds are uncalibrated.",
    ),
}


# ===========================================================================
# falsification -- Gate A item 3, done rather than asserted
# ===========================================================================
#
# "Every check demonstrated non-vacuous" means: suppress the behaviour the check asserts
# and confirm the check fails. 037 shipped four checks that passed for the wrong reason,
# and a check that passes vacuously is worse than no check because it is also evidence
# (§8.4). So this is mechanical, not a note in a PR body.
#
# How: the checker reads a tree on disk and nothing else, so a corruption of that tree is
# a suppression of the behaviour. Each corruption below edits one artifact in a copy of a
# case's tree, in a way a plausibly-wrong node would have produced. Every check that
# passed on the real tree then has to FAIL under at least one corruption. One that no
# corruption can break is reported as VACUOUS and fails the gate.
#
# A corruption returns a note, or None when it does not apply to that case -- which is
# ordinary, since a case with no report has nothing for a report corruption to do.

def _edit_json(path, change):
    if not os.path.exists(path):
        return None
    payload = json.load(open(path))
    note = change(payload)
    if note is None:
        return None
    json.dump(payload, open(path, "w"), indent=2)
    return note


def _edit_text(path, change):
    if not os.path.exists(path):
        return None
    before = open(path).read()
    after = change(before)
    if after is None or after == before:
        return None
    open(path, "w").write(after)
    return True


def x_drop_inventory(d):
    path = os.path.join(d, "n1", "outputs", "design_inventory.json")
    if not os.path.exists(path):
        return None
    os.remove(path)
    return "deleted design_inventory.json"


def x_drop_scores(d):
    path = os.path.join(d, "n4", "outputs", "design_scores.json")
    if not os.path.exists(path):
        return None
    os.remove(path)
    return "deleted design_scores.json"


def x_drop_report(d):
    path = os.path.join(d, "n5", "outputs", "report.html")
    if not os.path.exists(path):
        return None
    os.remove(path)
    return "deleted report.html"


def x_report_js_only(d):
    """What 037's report was: everything in a renderer, nothing server-side (cb418341)."""
    path = os.path.join(d, "n5", "outputs", "report.html")
    kept = []

    def change(html):
        for block in re.findall(r"<script\b.*?</script>", html, re.S | re.I):
            kept.append(block)
        return "<html><body><div id=\"app\"></div>" + "".join(kept) + "</body></html>"

    return ("moved the whole page into its <script> blocks"
            if _edit_text(path, change) else None)


def x_wrong_binder_chain(d):
    def change(payload):
        row = payload["designs"][0]
        others = [t["chain"] for t in row.get("target_chains") or []]
        if not others:
            return None
        row["binder_chain"] = others[0]
        row["designed_chain_source"] = "assumed-first-chain"
        return f"binder_chain changed to {others[0]} and the source to assumed-first-chain"

    return _edit_json(os.path.join(d, "n1", "outputs", "design_inventory.json"), change)


def x_truncate_target(d):
    def change(payload):
        row = payload["designs"][0]
        targets = row.get("target_chains") or []
        if not targets or len(targets[0]["sequence"]) < 2:
            return None
        before = len(targets[0]["sequence"])
        targets[0]["sequence"] = targets[0]["sequence"][:-1]
        targets[0]["residue_ids"] = targets[0]["residue_ids"][:-1]
        targets[0]["gaps"] = []
        return (f"target chain {targets[0]['chain']} trimmed from {before} to "
                f"{before - 1} residues and its gaps list emptied")

    return _edit_json(os.path.join(d, "n1", "outputs", "design_inventory.json"), change)


def x_renumber_from_one(d):
    def change(payload):
        row = payload["designs"][0]
        targets = row.get("target_chains") or []
        if not targets or targets[0]["residue_ids"][0] == 1:
            return None
        n = len(targets[0]["residue_ids"])
        targets[0]["residue_ids"] = list(range(1, n + 1))
        return "target residue ids renumbered from 1, as a node matching on number would"

    return _edit_json(os.path.join(d, "n1", "outputs", "design_inventory.json"), change)


def x_zero_unknown(d):
    def change(payload):
        row = payload["designs"][0]
        if not row.get("n_unknown"):
            return None
        was = row["n_unknown"]
        row["n_unknown"] = 0
        return f"n_unknown reset from {was} to 0"

    return _edit_json(os.path.join(d, "n1", "outputs", "design_inventory.json"), change)


def x_unskip(d):
    changed = []
    for step in ("n1", "n1b"):
        def change(payload):
            rows = [r for r in payload["designs"] if r.get("skipped")]
            if not rows:
                return None
            for row in rows:
                row["skipped"] = False
                row["skip_reason"] = None
            payload["n_skipped"] = 0
            return True

        if _edit_json(os.path.join(d, step, "outputs", "design_inventory.json"), change):
            changed.append(step)
    return f"un-skipped every skipped design in {changed}" if changed else None


def x_batch_fasta(d):
    """One file holding every record: 026's undistinguished directory (818954e9)."""
    where = os.path.join(d, "n2", "outputs")
    if not os.path.isdir(where):
        return None
    files = sorted(n for n in os.listdir(where)
                   if n.startswith("mono_") and n.endswith(".fasta"))
    if len(files) < 2:
        return None
    body = "".join(open(os.path.join(where, n)).read() for n in files)
    for name in files:
        os.remove(os.path.join(where, name))
    open(os.path.join(where, "batch.fasta"), "w").write(body)
    return f"collapsed {len(files)} per-design FASTAs into one batch.fasta"


def x_second_record_in_mono(d):
    """Two records in one per-design file.

    Not the same corruption as batch-fasta, and the difference matters: batch-fasta needs
    at least two designs to merge, so on a one-design case nothing changed and
    mono_record_count came out vacuous. This one reaches it. The predictor would fold both
    records and name its output after the first, which is 026's conflation arriving inside
    a single file (818954e9).
    """
    where = os.path.join(d, "n2", "outputs")
    if not os.path.isdir(where):
        return None
    for name in sorted(os.listdir(where)):
        if name.startswith("mono_") and name.endswith(".fasta"):
            with open(os.path.join(where, name), "a") as handle:
                handle.write(">second_record\nAAAAAAAAAA\n")
            return f"appended a second record to {name}"
    return None


def x_blank_fusion_disclosure(d):
    """Fold the gapped target and say nothing about it. 8e4772a7's defect, undisclosed.

    The gap stays in the inventory, the complex is still declared, and only the disclosure
    goes -- which is the one thing fusion_disclosed is about. Keeping the rest intact is
    the point: a corruption that also removes the shape would make the check fail for a
    reason that is not its own.
    """
    path = os.path.join(d, "n2", "outputs", "prediction_inputs.json")
    if not os.path.exists(path):
        return None
    payload = json.load(open(path))
    target = payload.get("target")
    if not isinstance(target, dict):
        return None
    before = list(target.get("fusion") or [])
    if not before:
        return None
    target["fusion"] = []
    json.dump(payload, open(path, "w"), indent=2)
    return (f"emptied target.fusion in prediction_inputs.json, dropping {len(before)} "
            f"disclosure(s): {before}")


def _cplx_fastas(d):
    where = os.path.join(d, "n2", "outputs")
    if not os.path.isdir(where):
        return []
    return [os.path.join(where, n) for n in sorted(os.listdir(where))
            if n.startswith("cplx_") and n.endswith(".fasta")]


def x_drop_colon_from_cplx(d):
    """Remove the chain break. The monomer prediction that reports itself as an interface.

    This is 026's bug in one byte: zero colons folds one chain, and every interface
    quantity computed off it is a number about a monomer. cplx_colon_count exists to catch
    exactly this, so this is what has to make it fail.
    """
    hit = {}
    for path in _cplx_fastas(d):
        body = open(path).read()
        open(path, "w").write(body.replace(":", ""))
        hit[os.path.basename(path)] = body.count(":")
    return (f"removed the chain break from {len(hit)} cplx_ fasta(s): {hit} colon(s) each"
            if hit else None)


def x_empty_cplx_half(d):
    """Keep the colon and empty the binder half, which is the other way to fold one chain.

    Counting colons is not enough on its own: 'TARGET:' has exactly one. 029 asserts both
    halves non-empty for this reason (04_fold_and_filter/script.py:161), and a check that
    only counted separators would pass this file while the prediction it produced was a
    monomer again.
    """
    hit = []
    for path in _cplx_fastas(d):
        lines = open(path).read().splitlines()
        for index, line in enumerate(lines):
            if ":" in line:
                lines[index] = line.split(":")[0] + ":"
                break
        open(path, "w").write("\n".join(lines) + "\n")
        hit.append(os.path.basename(path))
    return f"emptied the binder half of {hit}, leaving the colon in place" if hit else None


def x_colon_in_mono(d):
    where = os.path.join(d, "n2", "outputs")
    if not os.path.isdir(where):
        return None
    for name in sorted(os.listdir(where)):
        if name.startswith("mono_") and name.endswith(".fasta"):
            lines = open(os.path.join(where, name)).read().splitlines()
            lines[1] = lines[1][:4] + ":" + lines[1][4:]
            open(os.path.join(where, name), "w").write("\n".join(lines) + "\n")
            return f"inserted a chain-break colon into {name}'s sequence"
    return None


def x_break_declared_lengths(d):
    def change(payload):
        rows = payload.get("inputs") or []
        if not rows:
            return None
        rows[0]["expected_residues"] = rows[0]["expected_residues"] + 1
        rows[0]["record_id"] = rows[0]["record_id"] + " two"
        rows[0]["chain_order"] = []
        rows[0]["chain_lengths"] = []
        return ("expected_residues bumped by one, chain_order and chain_lengths emptied, "
                "and a space put in the record id")

    return _edit_json(os.path.join(d, "n2", "outputs", "prediction_inputs.json"), change)


def x_add_host_url(d):
    """026's defect, reinstated: --host-url, defaulting to the public API (02_predict:58)."""
    def change(payload):
        rows = payload.get("predictions") or []
        if not rows:
            return None
        for prediction in rows:
            prediction["argv"] = (prediction.get("argv") or []) + [
                "--host-url", "https://api.colabfold.com"
            ]
        payload["host_url_passed"] = True
        return "added --host-url https://api.colabfold.com to every argv"

    return _edit_json(_manifest_path(d), change)


def x_wrong_model_type(d):
    """Kept separate from add-host-url so each argv check names its own falsifier."""
    def change(payload):
        touched = 0
        for prediction in payload.get("predictions") or []:
            argv = prediction.get("argv") or []
            for flag, wrong in (("--model-type", "alphafold2_ptm_WRONG"),
                                ("--msa-mode", "mmseqs2_uniref_env")):
                if flag in argv:
                    argv[argv.index(flag) + 1] = wrong
                    touched += 1
        if not touched:
            return None
        return (f"changed {touched} resolved flag value(s): --model-type to a name no "
                "predictor has and --msa-mode to the one that needs the network")

    return _edit_json(_manifest_path(d), change)


def _manifest_path(d):
    for where in (("n3", "outputs"), ("pred",)):
        path = os.path.join(d, *where, "prediction_manifest.json")
        if os.path.exists(path):
            return path
    return os.path.join(d, "n3", "outputs", "prediction_manifest.json")


def x_mark_everything_ok(d):
    def change(payload):
        failed = [p for p in payload.get("predictions") or [] if p["status"] != "ok"]
        if not failed:
            return None
        for prediction in failed:
            prediction["status"] = "ok"
            prediction["reason"] = None
            prediction.pop("log", None)
        payload["n_failed"] = 0
        payload["n_ok"] = len(payload.get("predictions") or [])
        return (f"marked {len(failed)} failed prediction(s) ok, dropped their reasons and "
                "their retained log")

    return _edit_json(_manifest_path(d), change)


def x_claim_stdout_for_buffered_tail(d):
    """Node 03 as it was before C6: an empty tail, reported as if it came from stdout.

    The precise defect, rather than the blunt one. `mark-everything-ok` also breaks the
    buffered-tail check, but it breaks it by destroying the whole failure record, which
    would equally break a check that only looked at the status. This corruption leaves
    everything else intact and changes only where the lines are said to have come from,
    so the check has to be reading the thing it claims to read.
    """
    def change(payload):
        touched = 0
        for prediction in payload.get("predictions") or []:
            if prediction.get("status") == "ok":
                continue
            prediction["tail_source"] = "the predictor's stdout"
            reason = prediction.get("reason") or ""
            if "Last lines, from" in reason:
                prediction["reason"] = (
                    reason.split("Last lines, from")[0] + "Last lines: "
                )
            touched += 1
        if not touched:
            return None
        return (f"rewrote {touched} failed prediction(s) to claim the tail came from "
                "stdout, and emptied the tail -- node 03 as it behaved before C6")

    return _edit_json(_manifest_path(d), change)


def x_drop_failed_prediction(d):
    def change(payload):
        rows = payload.get("predictions") or []
        failed = [p for p in rows if p["status"] != "ok"]
        if not failed or len(rows) < 2:
            return None
        payload["predictions"] = [p for p in rows if p["status"] == "ok"]
        payload["n_declared"] = len(payload["predictions"])
        payload["n_failed"] = 0
        return f"dropped {len(failed)} failed prediction(s) from the manifest entirely"

    return _edit_json(_manifest_path(d), change)


def x_shift_rmsd(d):
    def change(payload):
        rows = [r for r in payload["designs"] if r.get("sc_rmsd_A") is not None]
        if not rows:
            return None
        for row in rows:
            row["sc_rmsd_A"] = round(row["sc_rmsd_A"] + 5.0, 4)
            row["sc_rmsd_n_residues"] = (row["sc_rmsd_n_residues"] or 0) + 2
        return "added 5.0 A to every sc_rmsd_A and 2 to the residue count it was taken over"

    return _edit_json(os.path.join(d, "n4", "outputs", "design_scores.json"), change)


def x_claim_target_as_binder_plddt(d):
    """Report the TARGET block's confidence as the binder's. 026's defect by another route.

    The one corruption that makes plddt_binder_in_complex fail for its OWN reason. Without
    it that check is falsified only by deleting design_scores.json, which proves it reads
    the file and nothing about whether it reads the right block -- the same weakness
    fusion_disclosed had (§8.4).
    """
    def change(payload):
        hit = 0
        for row in payload["designs"]:
            if (row.get("plddt_binder_in_complex") is not None
                    and row.get("plddt_target") is not None):
                row["plddt_binder_in_complex"] = row["plddt_target"]
                hit += 1
        if not hit:
            return None
        return (f"plddt_binder_in_complex set to plddt_target on {hit} row(s), which is "
                "what slicing the complex array from the wrong end produces")

    return _edit_json(os.path.join(d, "n4", "outputs", "design_scores.json"), change)


def x_swap_plddt(d):
    def change(payload):
        swapped = 0
        for row in payload["designs"]:
            if row.get("plddt_target") is not None:
                row["plddt_binder"], row["plddt_target"] = (
                    row["plddt_target"], row["plddt_binder"]
                )
                swapped += 1
        if not swapped:
            return None
        return (f"swapped plddt_binder and plddt_target on {swapped} row(s), which is what "
                "slicing the array from the wrong end produces")

    return _edit_json(os.path.join(d, "n4", "outputs", "design_scores.json"), change)


def x_shift_pae(d):
    def change(payload):
        rows = [r for r in payload["designs"] if r.get("pae_interaction") is not None]
        if not rows:
            return None
        for row in rows:
            row["pae_interaction"] = 6.0
        return ("pae_interaction set to 6.0, which is what averaging the target->binder "
                "block alone gives")

    return _edit_json(os.path.join(d, "n4", "outputs", "design_scores.json"), change)


def x_sentinel_nulls(d):
    """029's defaults, reinstated: 0.0 for a missing iptm, inf for a missing PAE (§8.7)."""
    def change(payload):
        touched = []
        for row in payload["designs"]:
            for metric, sentinel in (("iptm", 0.0), ("pae_interaction", 1e400),
                                     ("plddt_binder", 0.0), ("plddt_target", 0.0),
                                     ("sc_rmsd_A", 0.0)):
                if row.get(metric) is None:
                    row[metric] = sentinel
                    touched.append(metric)
        if not touched:
            return None
        payload["null_reasons"] = {}
        return (f"replaced every null with 029's sentinel ({sorted(set(touched))}) and "
                "emptied null_reasons")

    return _edit_json(os.path.join(d, "n4", "outputs", "design_scores.json"), change)


def x_drop_null_reasons(d):
    def change(payload):
        if not payload.get("null_reasons"):
            return None
        payload["null_reasons"] = {}
        return "emptied null_reasons, leaving every null unexplained"

    return _edit_json(os.path.join(d, "n4", "outputs", "design_scores.json"), change)


def x_gate_on_target(d):
    def change(payload):
        rows = [r for r in payload["designs"] if r.get("plddt_target") is not None]
        if not rows:
            return None
        payload["gating_metrics"] = (payload.get("gating_metrics") or []) + [
            "plddt_target"
        ]
        payload["diagnostic_metrics"] = [
            m for m in payload.get("diagnostic_metrics") or [] if m != "plddt_target"
        ]
        for row in rows:
            row["passes"] = False
            row["failed_thresholds"] = "plddt_target"
        return "made plddt_target gating and failed every design on it"

    return _edit_json(os.path.join(d, "n4", "outputs", "design_scores.json"), change)


def x_drop_thresholds(d):
    def change(payload):
        if not payload.get("thresholds"):
            return None
        payload["thresholds"] = {}
        payload["thresholds_are_calibrated"] = True
        for row in payload["designs"]:
            row.pop("max_sc_rmsd_a", None)
            row.pop("min_plddt_binder", None)
        return ("removed the thresholds, the per-row threshold columns, and the "
                "uncalibrated disclaimer")

    return _edit_json(os.path.join(d, "n4", "outputs", "design_scores.json"), change)


def x_scramble_merge(d):
    """Rotate each row's generation metrics and sequence onto the next design."""
    def change(payload):
        rows = payload["designs"]
        if len(rows) < 2:
            return None
        keys = [k for k in rows[0] if k.startswith("gen_")] + [
            "binder_sequence", "plddt_binder"
        ]
        values = [{k: r.get(k) for k in keys} for r in rows]
        for row, taken in zip(rows, values[1:] + values[:1]):
            row.update(taken)
        return f"rotated {keys} onto the next design, as a positional merge would"

    return _edit_json(os.path.join(d, "n4", "outputs", "design_scores.json"), change)


def x_clear_reasons(d):
    def change(payload):
        rows = [r for r in payload["designs"] if r.get("reason")]
        if not rows:
            return None
        for row in rows:
            row["reason"] = None
            row["status"] = "scored"
        return f"cleared the reason on {len(rows)} failed row(s) and called them scored"

    return _edit_json(os.path.join(d, "n4", "outputs", "design_scores.json"), change)


def x_pass_everything(d):
    def change(payload):
        rows = [r for r in payload["designs"] if not r.get("passes")]
        if not rows:
            return None
        for row in rows:
            row["passes"] = True
        payload["n_passing"] = len(payload["designs"])
        return f"marked {len(rows)} non-passing design(s) as passing"

    return _edit_json(os.path.join(d, "n4", "outputs", "design_scores.json"), change)


def x_lexical_order(d):
    """Reorder the report's rows lexically: design_1, design_10, design_11, design_2 (#51)."""
    path = os.path.join(d, "n5", "outputs", "report.html")
    if not os.path.exists(path):
        return None
    html = open(path).read()
    table = re.search(r"(<table id=\"metrics\".*?<tbody>)(.*?)(</tbody>)", html, re.S)
    if not table:
        return None
    rows = re.findall(r"<tr\b.*?</tr>", table.group(2), re.S)
    if len(rows) < 3:
        return None

    def key(row):
        found = re.search(r">(design_\d+)<", row)
        return found.group(1) if found else row

    reordered = sorted(rows, key=key)
    if reordered == rows:
        return None
    html = html.replace(table.group(2), "\n".join(reordered))
    open(path, "w").write(html)
    csv_path = os.path.join(d, "n5", "outputs", "validation_metrics.csv")
    if os.path.exists(csv_path):
        lines = open(csv_path).read().splitlines()
        body = sorted(lines[1:], key=lambda l: l.split(",")[1] if "," in l else l)
        open(csv_path, "w").write("\n".join([lines[0]] + body) + "\n")
    return "reordered the report's rows and the CSV's lexically"


def x_drop_row_from_report(d):
    path = os.path.join(d, "n5", "outputs", "report.html")
    if not os.path.exists(path):
        return None
    html = open(path).read()
    names = re.findall(r">(design_\d+)<", html)
    if len(set(names)) < 2:
        return None
    victim = sorted(set(names), key=natural_key)[-1]
    rows = re.findall(r"<tr\b(?:(?!</tr>).)*?" + victim + r"(?:(?!</tr>).)*?</tr>",
                      html, re.S)
    if not rows:
        return None
    for row in rows:
        html = html.replace(row, "")
    open(path, "w").write(html)
    return f"removed {victim}'s row from the report"


def x_strip_source_from_report(d):
    """Keep the value out of everything a reader sees, leaving it in the manifest only."""
    scores = os.path.join(d, "n4", "outputs", "design_scores.json")
    if not os.path.exists(scores):
        return None
    payload = json.load(open(scores))
    sources = {r.get("designed_chain_source") for r in payload["designs"]}
    sources.discard(None)
    if not sources:
        return None
    hit = False
    for name in ("report.html", "validation_metrics.csv"):
        path = os.path.join(d, "n5", "outputs", name)
        if not os.path.exists(path):
            continue
        text = open(path).read()
        for source in sources:
            text = text.replace(source, "")
        open(path, "w").write(text)
        hit = True
    return f"removed {sorted(sources)} from the report and the CSV" if hit else None


def x_fasta_colon_and_banner(d):
    path = os.path.join(d, "n5", "outputs", "passing_sequences.fasta")

    def change(text):
        lines = text.splitlines()
        if not any(l.startswith(">") for l in lines):
            return None
        out = ["; workflow-040 passing sequences"]
        for line in lines:
            out.append(line.replace("UNSCREENED -", "UNSCREENED:")
                       if line.startswith(">") else line)
        return "\n".join(out) + "\n"

    return ("put a ';' banner at the top and turned the tag's dash into a colon"
            if _edit_text(path, change) else None)


def x_untag_fasta(d):
    path = os.path.join(d, "n5", "outputs", "passing_sequences.fasta")

    def change(text):
        if "UNSCREENED" not in text:
            return None
        lines = [l.split(" | ")[0] if l.startswith(">") else l
                 for l in text.splitlines()]
        return "\n".join(lines) + "\n"

    return ("removed the UNSCREENED tag from every header"
            if _edit_text(path, change) else None)


def x_export_a_failure(d):
    """Put a non-passing design into the exported FASTA."""
    scores = os.path.join(d, "n4", "outputs", "design_scores.json")
    path = os.path.join(d, "n5", "outputs", "passing_sequences.fasta")
    if not (os.path.exists(scores) and os.path.exists(path)):
        return None
    payload = json.load(open(scores))
    losers = [r for r in payload["designs"] if not r.get("passes")]
    if not losers:
        return None
    row = losers[0]
    with open(path, "a") as handle:
        handle.write(f">{row['design']} | UNSCREENED - leaked\n"
                     f"{row.get('binder_sequence') or 'AAAA'}\n")
    return f"appended the non-passing {row['design']} to passing_sequences.fasta"


def x_remove_one_mono(d):
    """Lose one per-design input while the declaration still claims it."""
    where = os.path.join(d, "n2", "outputs")
    if not os.path.isdir(where):
        return None
    for name in sorted(os.listdir(where)):
        if name.startswith("mono_") and name.endswith(".fasta"):
            os.remove(os.path.join(where, name))
            return f"deleted {name} while prediction_inputs.json still declares it"
    return None


def x_blank_node_log(d, record):
    """Take the words out of the nodes' output, leaving every artifact intact.

    Takes the record as well as the directory, because a node's stdout and stderr are kept
    in run.json rather than only in the .log files beside the tree -- the first version of
    this corruption emptied the files alone and changed nothing, which is why two
    log-reading checks were credited to deleting an artifact instead.
    """
    hit = []
    for name in sorted(os.listdir(d)):
        if name.endswith(".log") and os.path.isfile(os.path.join(d, name)):
            open(os.path.join(d, name), "w").write("")
            hit.append(name)
    for step in record.get("steps") or []:
        if step.get("stdout") or step.get("stderr"):
            step["stdout"] = ""
            step["stderr"] = ""
            hit.append(step["step"])
    return f"silenced {sorted(set(hit))}" if hit else None


# ORDER MATTERS. falsify() stops at the first corruption that breaks a check, so whichever
# comes first is the one recorded against it -- and "the check noticed the file was gone"
# is a much weaker demonstration than "the check noticed the number was wrong". The
# specific, value-level corruptions therefore come first and the blunt deletions last, so
# a check credited to a deletion is one that nothing more precise could reach, which is
# itself worth seeing.
CORRUPTIONS = [
    # wrong values, right files
    #
    # blank-fusion-disclosure comes FIRST, ahead of truncate-target. Both make
    # fusion_disclosed fail, but only one makes it fail for its own reason: truncate-target
    # empties the gaps list, so the check fails because the case is no longer a gapped
    # target at all. That proves the check reads the inventory and says nothing about
    # whether it reads the disclosure. falsify() stops at the first corruption that flips a
    # check, so the order here decides which reason gets recorded -- and a check falsified
    # for the wrong reason is only a little better than one that is never falsified (§8.4).
    ("blank-fusion-disclosure", x_blank_fusion_disclosure),
    ("wrong-binder-chain", x_wrong_binder_chain),
    ("truncate-target", x_truncate_target),
    ("renumber-from-one", x_renumber_from_one),
    ("zero-unknown", x_zero_unknown),
    ("unskip", x_unskip),
    ("second-record-in-mono", x_second_record_in_mono),
    ("remove-one-mono", x_remove_one_mono),
    ("colon-in-mono", x_colon_in_mono),
    # Value-level before deletion, per the ordering note on this table: falsify() stops at
    # the first corruption that flips a check, and a deleted file would shadow both of
    # these with a weaker reason.
    ("empty-cplx-half", x_empty_cplx_half),
    ("drop-colon-from-cplx", x_drop_colon_from_cplx),
    ("batch-fasta", x_batch_fasta),
    ("break-declared-lengths", x_break_declared_lengths),
    ("wrong-model-type", x_wrong_model_type),
    ("add-host-url", x_add_host_url),
    # Before mark-everything-ok, which also breaks the buffered-tail check but does it by
    # destroying the whole failure record. falsify() stops at the first corruption that
    # flips a check, so the precise one has to come first or the check gets credited to
    # the blunt one and nobody learns anything.
    ("claim-stdout-for-buffered-tail", x_claim_stdout_for_buffered_tail),
    ("mark-everything-ok", x_mark_everything_ok),
    ("drop-failed-prediction", x_drop_failed_prediction),
    ("shift-rmsd", x_shift_rmsd),
    ("claim-target-as-binder-plddt", x_claim_target_as_binder_plddt),
    ("swap-plddt", x_swap_plddt),
    ("shift-pae", x_shift_pae),
    ("sentinel-nulls", x_sentinel_nulls),
    ("drop-null-reasons", x_drop_null_reasons),
    ("gate-on-target", x_gate_on_target),
    ("drop-thresholds", x_drop_thresholds),
    ("scramble-merge", x_scramble_merge),
    ("clear-reasons", x_clear_reasons),
    ("pass-everything", x_pass_everything),
    ("lexical-order", x_lexical_order),
    ("drop-row-from-report", x_drop_row_from_report),
    ("strip-source-from-report", x_strip_source_from_report),
    ("report-js-only", x_report_js_only),
    ("fasta-colon-and-banner", x_fasta_colon_and_banner),
    ("untag-fasta", x_untag_fasta),
    ("export-a-failure", x_export_a_failure),
    ("blank-node-log", x_blank_node_log),
    # and only then, the artifact missing entirely
    ("drop-inventory", x_drop_inventory),
    ("drop-scores", x_drop_scores),
    ("drop-report", x_drop_report),
]

# A few claims do not live in an artifact at all, so no corruption of the tree can reach
# them. Each is named here with what DOES falsify it, because "no corruption broke it" has
# to mean one of two very different things and the gate must not confuse them.
NOT_IN_THE_TREE = {
    "sigterm_ignoring_child_killed":
        "the claim is about a process that no longer exists, recorded live by the harness "
        "as /proc/<pid> existing or not. It was falsified for real during the build: the "
        "first stand-in put its worker behind setsid, the worker survived, and the check "
        "failed. See validation_pipeline.md §6. It is also the one check the IMAGE run "
        "cannot make -- a container pid means nothing on the host -- so there it reports a "
        "gap, and the host run is what establishes it.",
    "no_network_reached":
        "partly in the tree -- add-host-url breaks the argv half -- and partly in the "
        "container's --network none, which is a property of the run rather than of any "
        "file. The check's own detail says which half it is establishing.",
    "rerun_idempotent_or_refused":
        "the claim is about three exit codes and a byte comparison the harness performed "
        "between runs, not about the final state of any file. Falsified by making node 02 "
        "accept the stale file, which is what 029 did (4d1a523e).",
    "expect_error":
        "the claim is a non-zero exit and a phrase in the message, so what falsifies it is "
        "a node that accepts the input -- which is a change to a node, not a corruption of "
        "the tree. cplx-input-missing-colon used to be the live demonstration: it ran "
        "every step to 0 and was reported as a gap for exactly that reason. It no longer "
        "is, because the guard it was aimed at now refuses, so that demonstration is gone "
        "and is not replaced by a weaker one. What stops 'refuse everything' from passing "
        "all 14 of these is structural rather than a corruption: boundary-value-accepted "
        "is an input sitting exactly on every limit that MUST be accepted, so a guard made "
        "strict enough to pass these fails that one.",
    "expect_ok":
        "the counter-entry to every refusal: a guard that refused this input would fail "
        "it. It is the one case in the suite whose whole job is to fail if the guards are "
        "made too strict.",
}


def falsify(run, by_id, root):
    """For each check that passed, find a corruption of the tree that makes it fail."""
    import copy as copy_module
    import shutil
    import tempfile

    broken_by = {}
    attempted = {}
    for record in run["cases"]:
        section, entry = by_id[record["id"]]
        specs = list(entry.get("expect") or [])
        if not specs:
            continue
        live = View(record)
        wanted = []
        for spec in specs:
            name = spec["check"]
            if name not in CHECKS:
                continue
            status, _ = CHECKS[name][0](live, spec)
            if status == PASS:
                wanted.append((name, spec))
                attempted.setdefault(name, set()).add(record["id"])
        if not wanted:
            continue
        work = tempfile.mkdtemp(prefix="falsify-", dir=root)
        try:
            for label, corrupt in CORRUPTIONS:
                if not any(name not in broken_by for name, _ in wanted):
                    break
                copy = os.path.join(work, label)
                shutil.copytree(record["dir"], copy)
                # Deep, because a corruption may edit the recorded stdout and stderr and
                # must not reach the record the real verdict was taken from.
                shadow = copy_module.deepcopy(record)
                shadow["dir"] = copy
                try:
                    note = (corrupt(copy, shadow)
                            if corrupt.__code__.co_argcount == 2 else corrupt(copy))
                except Exception as exc:
                    note = f"the corruption itself raised {type(exc).__name__}: {exc}"
                if not note:
                    shutil.rmtree(copy, ignore_errors=True)
                    continue
                view = View(shadow)
                for name, spec in wanted:
                    if name in broken_by:
                        continue
                    try:
                        status, detail = CHECKS[name][0](view, spec)
                    except Exception as exc:
                        status = FAIL
                        detail = f"raised {type(exc).__name__}: {exc}"
                    if status == FAIL:
                        broken_by[name] = (label, record["id"], note,
                                           " ".join(str(detail).split())[:150])
                shutil.rmtree(copy, ignore_errors=True)
        finally:
            shutil.rmtree(work, ignore_errors=True)
    return broken_by, attempted


# ===========================================================================
# running the table
# ===========================================================================

def judge_error_case(case, entry):
    """A case whose whole claim is a refusal."""
    want = entry["expect_error"]
    failures = case.failures()
    if not failures:
        return FAIL, (f"every step exited 0, so nothing refused anything. Expected a "
                      f"message containing {want!r}. Steps: "
                      + ", ".join(f"{s['step']}={s['rc']}"
                                  for s in case.record["steps"]))
    first = failures[0]
    body = (first.get("stdout") or "") + (first.get("stderr") or "")
    line = next((l.strip() for l in body.splitlines() if want in l), None)
    left = case.files(first["step"], "outputs")
    detail = (f"{first['step']} exited {first['rc']}; its ./outputs/ holds {left or 'no'} "
              f"file(s); the message says: "
              + " ".join((body.strip().splitlines() or [""])[-1].split())[:300])
    if line is None:
        return FAIL, f"expected a message containing {want!r}. " + detail
    return PASS, f"refused at {first['step']} with {want!r} in: {line[:240]}"


def judge_ok_case(case, entry):
    failures = case.failures()
    scored = case.scores()
    detail = ("every step exited 0: "
              + ", ".join(f"{s['step']}={s['rc']}" for s in case.record["steps"])
              + (f"; node 04 scored {scored.get('n_scored')} design(s) with "
                 f"n_passing={scored.get('n_passing')} at thresholds "
                 f"{scored.get('thresholds')}" if scored else ""))
    if failures:
        first = failures[0]
        body = ((first.get("stderr") or "") + (first.get("stdout") or "")).strip()
        return FAIL, (f"{first['step']} exited {first['rc']} on an input that is exactly "
                      f"on its limit and must be accepted: "
                      + " ".join(body.splitlines()[-1:])[:300])
    if not scored:
        return FAIL, "nothing was scored, so acceptance was not shown. " + detail
    return PASS, detail


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    root = os.path.abspath(sys.argv[1])
    verbose = "--verbose" in sys.argv
    run = json.load(open(os.path.join(root, "run.json")))

    width = 78
    print("=" * width)
    print("workflow-040 fixture suite")
    if run.get("is_gate"):
        print(f"  ran in {run['image']} -- THIS IS THE GATE")
    else:
        print("  NOT THE GATE: ran on " + run["runner"])
        print(f"  because: {run.get('runner_reason')}")
        print("  Tier 1 is only satisfied by a run in the image.")
    print(f"  staged: {run['root']}")
    print(f"  cases:  {run['n_cases']}")
    print("=" * width)

    invoked = {}
    # The fixture entries live beside the run, not in it, so the two are matched by id.
    # Deliberately re-read here rather than copied into run.json: a changed expectation
    # should be able to re-judge a run that is already on disk without re-running it.
    payload = json.load(open(os.path.join(HERE, "fixtures.json")))
    by_id = {}
    for section in ("ingest", "predict", "scoring", "rejections", "controls"):
        for entry in payload[section]:
            by_id[entry["id"]] = (section, entry)

    verdicts = []
    for record in run["cases"]:
        case = View(record)
        section, entry = by_id[record["id"]]
        rows = []
        if record.get("harness_error"):
            rows.append(("harness", FAIL,
                         "the harness raised before the nodes could be judged: "
                         + record["harness_error"]))
        elif "expect_error" in entry:
            status, detail = judge_error_case(case, entry)
            rows.append(("expect_error", status, detail))
            invoked.setdefault("expect_error", []).append(record["id"])
        elif entry.get("expect_ok"):
            status, detail = judge_ok_case(case, entry)
            rows.append(("expect_ok", status, detail))
            invoked.setdefault("expect_ok", []).append(record["id"])
        for spec in entry.get("expect") or []:
            name = spec["check"]
            invoked.setdefault(name, []).append(record["id"])
            if name not in CHECKS:
                rows.append((name, FAIL, "no implementation for this check"))
                continue
            try:
                status, detail = CHECKS[name][0](case, spec)
            except Exception as exc:
                import traceback
                status = FAIL
                detail = (f"the check raised {type(exc).__name__}: {exc}\n"
                          + traceback.format_exc().strip().splitlines()[-3])
            rows.append((name, status, detail))
        # Phase 2 is specced and not built. A phase-2 claim that does not hold is a gap in
        # what exists, not a defect in it -- but the observation is kept verbatim, because
        # several of these already hold and the ones that do not say why.
        if entry.get("phase") == 2:
            rows = [(n, GAP if s == FAIL else s,
                     ("phase 2, not built -- recorded, not counted against the gate: "
                      + d) if s == FAIL else d)
                    for n, s, d in rows]
        verdicts.append((section, entry, case, rows))

    counts = {FAIL: 0, GAP: 0, PASS: 0}
    for section, entry, case, rows in verdicts:
        worst = min((RANK[s] for _, s, _ in rows), default=RANK[PASS])
        headline = [k for k, v in RANK.items() if v == worst][0]
        counts[headline] += 1
        show = verbose or headline != PASS
        mark = {PASS: "ok  ", FAIL: "FAIL", GAP: "gap "}[headline]
        print(f"\n[{mark}] {section}/{entry['id']}"
              + (f"  (phase {entry['phase']})" if entry.get("phase") == 2 else ""))
        for name, status, detail in rows:
            if status == PASS and not show:
                continue
            print(f"    {status:4}  {name}")
            for line in wrap(detail, width - 10):
                print(f"          {line}")

    print("\n" + "=" * width)
    print(f"cases: {counts[PASS]} pass, {counts[FAIL]} FAIL, {counts[GAP]} gap "
          f"(of {len(verdicts)})")

    named = set()
    for section in ("ingest", "predict", "scoring", "rejections"):
        for entry in payload[section]:
            for spec in entry.get("expect") or []:
                named.add(spec["check"])
    unimplemented = sorted(named - set(CHECKS))
    unused = sorted(set(CHECKS) - set(invoked))
    print(f"checks: {len(CHECKS)} implemented, {len(named)} named in fixtures.json, "
          f"{len([k for k in invoked if k in CHECKS])} invoked")
    if unimplemented:
        print(f"  NAMED BUT NOT IMPLEMENTED: {unimplemented}")
    if unused:
        print(f"  IMPLEMENTED BUT NEVER INVOKED: {unused}")

    # Gate A item 3, done rather than asserted. Skipped only with --no-falsify, and the
    # gate then says so.
    broken_by, attempted = ({}, {})
    if "--no-falsify" not in sys.argv:
        print("\nfalsifying every passing check by corrupting the tree it read...")
        broken_by, attempted = falsify(run, by_id, root)

    vacuous = sorted(
        name for name in attempted
        if name not in broken_by and name not in NOT_IN_THE_TREE
    )

    print("\nnon-vacuity -- what would make each check fail, and what did:")
    for name in sorted(CHECKS):
        cases = invoked.get(name, [])
        print(f"  {name} ({len(cases)} case(s): {', '.join(cases) or 'NONE'})")
        for line in wrap(CHECKS[name][1], width - 6):
            print(f"      {line}")
        if name in broken_by:
            label, where, note, detail = broken_by[name]
            for line in wrap(f"FALSIFIED by '{label}' on {where}: {note} -> the check "
                             f"then failed with: {detail}", width - 6):
                print(f"      {line}")
        elif name in NOT_IN_THE_TREE:
            for line in wrap("NOT IN THE TREE: " + NOT_IN_THE_TREE[name], width - 6):
                print(f"      {line}")
        elif name in attempted:
            for line in wrap("VACUOUS: no corruption of the tree made this check fail, "
                             "and it is not listed in NOT_IN_THE_TREE. Either the check "
                             "is measuring nothing or the corruption it needs is "
                             "missing.", width - 6):
                print(f"      {line}")
        else:
            for line in wrap("not falsified: it had no passing invocation to corrupt "
                             "(every case it appears in is a gap).", width - 6):
                print(f"      {line}")
    for pseudo in ("expect_error", "expect_ok"):
        if pseudo in invoked:
            print(f"  {pseudo} ({len(invoked[pseudo])} case(s))")
            for line in wrap("NOT IN THE TREE: " + NOT_IN_THE_TREE[pseudo], width - 6):
                print(f"      {line}")

    print("=" * width)
    falsified = len(broken_by)
    print(f"falsification: {falsified} check(s) broken by a corruption, "
          f"{len([n for n in attempted if n in NOT_IN_THE_TREE])} not in the tree and "
          f"accounted for, {len(vacuous)} VACUOUS")
    if vacuous:
        print(f"  VACUOUS: {vacuous}")
    gate = (counts[FAIL] == 0 and not unimplemented and not unused and not vacuous
            and "--no-falsify" not in sys.argv)
    print("Gate A: zero FAIL, every check implemented, invoked and falsifiable: "
          + ("yes" if gate else "NO"))
    if not run.get("is_gate"):
        print("Gate A NOT satisfied regardless: this run was not in "
              f"{run.get('image') or 'the image'}, so it does not count. Everything above "
              "is a dry run.")
    return 0 if counts[FAIL] == 0 and not vacuous else 1


def wrap(text, width):
    words, line, out = str(text).split(), "", []
    for word in words:
        if len(line) + len(word) + 1 > width:
            out.append(line)
            line = word
        else:
            line = (line + " " + word).strip()
    if line:
        out.append(line)
    return out


if __name__ == "__main__":
    sys.exit(main())
