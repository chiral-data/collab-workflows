#!/usr/bin/env python3
"""Assert what the nodes REPORTED against what the input actually CONTAINS.

The reason this file exists: every fixture in the manifest exits 0 on nodes 01 and 02.
Exit codes detect none of the input-shape bugs, because the failure mode is silent --
a residue that cannot be addressed, an ensemble passed through whole, a cofactor
dropped without mention. So the suite compares node output against a measurement of
the input taken by a separate parser (tests/probe_entry.py).

Expectations live in tests/fixtures.json, not here. An entry marked "status": "open" is
a known gap: it is reported as GAP and does not fail the run. Deleting that field makes
the suite enforce it.

Run:  python3 -I tests/check_fixtures.py <results-dir> <downloads-dir>

Both directories are produced by tests/run_fixtures.sh.
"""
import glob
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from probe_entry import atom_site, probe  # noqa: E402

PASS, FAIL, GAP, SKIP = "PASS", "FAIL", "GAP", "SKIP"


class Fixture:
    """Everything measured or read for one entry, loaded once."""

    def __init__(self, results, downloads, spec):
        self.id = spec["id"]
        self.spec = spec
        base = os.path.join(results, self.id)
        self.inv_path = os.path.join(base, "n1", "outputs", "target_inventory.json")
        self.out_cif = os.path.join(base, "n1", "outputs", "target.cif")
        self.src_cif = os.path.join(downloads, spec.get("source", self.id + ".cif"))
        self.n1_log = self._read(os.path.join(base, "n1.log"))
        self.n2_log = self._read(os.path.join(base, "n2.log"))
        self.spec_yaml = self._read(
            os.path.join(base, "n2", "outputs", "design_spec.yaml"))
        self.inventory = None
        self.chains = {}
        if os.path.exists(self.inv_path):
            with open(self.inv_path) as handle:
                self.inventory = json.load(handle)
            self.chains = {c["label_asym_id"]: c for c in self.inventory["chains"]}
        self.truth = probe(self.src_cif) if os.path.exists(self.src_cif) else {}

    @staticmethod
    def _read(path):
        if not os.path.exists(path):
            return ""
        with open(path, errors="replace") as handle:
            return handle.read()


# --- checks: each returns (status, detail) ---------------------------------

def check_label_to_auth_complete(fx, _):
    """Every residue of every chain must have an entry. label_seq_id is unique by
    construction, so a short map means residues that node 02 will refuse to resolve."""
    short = []
    for cid, chain in fx.chains.items():
        n_map, n_res = len(chain.get("label_to_auth", {})), chain["n_residues"]
        if n_map != n_res:
            short.append(f"chain {cid}: {n_map}/{n_res} ({n_res - n_map} unaddressable)")
    if short:
        return FAIL, "; ".join(short)
    total = sum(len(c.get("label_to_auth", {})) for c in fx.chains.values())
    return PASS, f"{total} entries across {len(fx.chains)} chain(s)"


def check_inserted_residues_addressable(fx, _):
    """Residues carrying an insertion code must be keys in label_to_auth, and the value
    must carry the code -- otherwise two residues report the same author number."""
    found, bad = 0, []
    for cid, chain in fx.chains.items():
        for residue in chain.get("residues", []):
            if not residue.get("ins_code"):
                continue
            found += 1
            key = str(residue["label_seq_id"])
            value = chain.get("label_to_auth", {}).get(key)
            expected = f"{residue['auth_seq_id']}{residue['ins_code']}"
            if value != expected:
                bad.append(f"{cid}:{key} -> {value!r} (want {expected!r})")
    if not found:
        return FAIL, "no residue carries an insertion code; fixture no longer fits"
    if bad:
        return FAIL, f"{len(bad)} wrong: " + "; ".join(bad[:4])
    return PASS, f"{found} inserted residue(s) addressable and uniquely labelled"


def check_insertion_codes_warned_at_use(fx, _):
    """Node 02 must name the residues an author-numbered hotspot cannot reach.

    The spec itself is not wrong: author numbering means the first occurrence, and that
    is what gets conditioned. What is checked here is that the loss is announced where
    the residue is named -- node 01's warning is about the structure and arrives in a
    different node's log -- and that the shadowed residues really are absent from the
    spec, so the message and the spec agree.

    Which residues are shadowed is measured from the inventory, not listed here, so the
    check still fits if the fixture's hotspot span changes.
    """
    tokens = fx.spec.get("params", {}).get("PARAM_HOTSPOT_RESIDUES", "")
    wanted = set()
    for token in tokens.split(","):
        match = re.match(r"^([A-Za-z]*)(\d+)(?:-(\d+))?$", token.strip())
        if match:
            start = int(match.group(2))
            wanted.update(range(start, int(match.group(3) or start) + 1))
    if not wanted:
        return FAIL, "fixture names no hotspot residues, so nothing can be shadowed"

    # Only the target chain: hotspots may name no other, and a residue numbered 100A on
    # a different chain of the same Fab would otherwise stand in for the one at issue.
    wanted_chain = fx.spec.get("params", {}).get("PARAM_TARGET_CHAIN", "")
    target = fx.chains.get(wanted_chain) or next(
        (c for c in fx.chains.values() if c.get("auth_asym_id") == wanted_chain), None)
    if target is None:
        return FAIL, f"target chain {wanted_chain!r} is not in the inventory"

    # The residue each author number actually resolves to is not shadowed by anything --
    # it is the one being conditioned. On 12E8 that is always the bare residue, so the
    # distinction is invisible there. On 1IGT author 100 has no bare form and resolves to
    # 100H, and demanding that node 02 call 100H unreachable demands a message that
    # contradicts the spec written beside it.
    resolved = {}
    for auth_num in wanted:
        label_seq = target.get("auth_to_label", {}).get(str(auth_num))
        if label_seq is not None:
            resolved[auth_num] = int(label_seq)

    shadowed, conditioned = [], []
    for residue in target.get("residues", []):
        if not residue.get("ins_code"):
            continue
        auth_num = int(residue["auth_seq_id"])
        if auth_num not in wanted:
            continue
        name = f"{residue['auth_seq_id']}{residue['ins_code']}"
        if resolved.get(auth_num) == residue["label_seq_id"]:
            conditioned.append((name, residue["label_seq_id"]))
        else:
            shadowed.append((name, residue["label_seq_id"]))
    if not shadowed:
        return FAIL, "no requested hotspot carries an insertion code; fixture no longer fits"

    warnings = [l for l in fx.n2_log.splitlines() if l.startswith("WARNING")]
    missing = [name for name, _ in shadowed
               if not any(name in line for line in warnings)]
    if missing:
        return FAIL, f"node 02 never named {', '.join(missing)} as unreachable"
    leaked = [str(seq) for _, seq in shadowed if f"/{seq}/" in fx.spec_yaml]
    if leaked:
        return FAIL, f"spec conditions label {', '.join(leaked)}, which auth cannot name"
    # The other half, and the half that was missing: an insertion-coded residue the
    # lookup landed on must be conditioned, and must not also be called unreachable.
    for name, label_seq in conditioned:
        if any(name in line for line in warnings):
            return FAIL, (f"node 02 calls {name} unreachable, but author "
                          f"{name[:-1]} resolves to it at label {label_seq}")
        if f"/{label_seq}/" not in fx.spec_yaml:
            return FAIL, (f"author {name[:-1]} resolves to {name} at label {label_seq}, "
                          "which the spec does not condition")
    detail = (f"{len(shadowed)} shadowed residue(s) named in node 02's log "
              "and absent from the spec")
    if conditioned:
        detail += (f"; {', '.join(n for n, _ in conditioned)} conditioned as the "
                   "residue(s) author numbering resolves to")
    return PASS, detail


def check_log_mentions_models(fx, _):
    models = fx.truth.get("models", 1)
    if models <= 1:
        return FAIL, "input is single-model; fixture no longer fits"
    if "model" not in fx.n1_log.lower():
        return FAIL, f"input has {models} models, node 01 never says so"
    return PASS, f"{models} models in input, node 01 reports it"


def check_single_model_output(fx, _):
    if not os.path.exists(fx.out_cif):
        return FAIL, "no target.cif"
    out = probe(fx.out_cif)
    got, src = out.get("models", 1), fx.truth.get("models", 1)
    if got != 1:
        return FAIL, f"target.cif still holds {got} models"
    return PASS, (f"{src} model(s) in, 1 out; "
                  f"{fx.truth.get('n_atom_rows')} -> {out.get('n_atom_rows')} atom rows")


def check_nonpolymer_named(fx, expect):
    """A non-polymer in the file must reach both the inventory and the printed output.

    Checked against the inventory rather than by grepping the log, because a two-letter
    component like NA matches inside NAG and inside 'nan' -- a log scan would have
    reported this as passing while the entity was still invisible.
    """
    want = expect.get("components", [])
    present = [c for c in want if c in fx.truth.get("nonpolymer_list", [])]
    if not present:
        return FAIL, f"{want} not in the input at all; fixture no longer fits"
    entities = fx.inventory.get("non_polymer") or []
    reported = {e["comp_id"] for e in entities}
    missing = [c for c in present if c not in reported]
    if missing:
        return FAIL, f"{present} in the input, missing from the inventory: {missing}"
    shown = [c for c in present if f" {c} " in fx.n1_log]
    if len(shown) != len(present):
        return FAIL, f"{present} inventoried but not printed"
    nearest = [f"{e['comp_id']}->{e['nearest_chain']}@{e['min_distance_A']}A"
               for e in entities if e["comp_id"] in present][:3]
    return PASS, f"{len(entities)} entit(y/ies) inventoried and printed; {nearest}"


def check_dropped_chains_named(fx, _):
    if len(fx.chains) <= 2:
        return FAIL, "fewer than 3 chains; fixture no longer fits"
    contig = ""
    for line in fx.spec_yaml.splitlines():
        if line.strip().startswith("contig:"):
            contig = line.strip()
            break
    dropped = [c for c in fx.chains if c not in contig]
    if not dropped:
        return PASS, "contig preserves every chain"
    # Each dropped chain must be named individually. A phrase match would pass on the
    # unrelated "derived contig ... discarded" note that a contig override prints.
    warnings = "\n".join(
        line for line in fx.n2_log.splitlines() if line.startswith("WARNING"))
    unnamed = [c for c in dropped if f"{c} (auth " not in warnings]
    detail = f"{len(dropped)} dropped {dropped}"
    if unnamed:
        return FAIL, detail + f"; not named in any warning: {unnamed}"
    return PASS, detail + "; each one named in a warning"


def check_chain_length(fx, expect):
    cid, want = expect["chain"], expect["n_residues"]
    chain = fx.chains.get(cid)
    if chain is None:
        return FAIL, f"chain {cid} absent (have {sorted(fx.chains)})"
    got = chain["n_residues"]
    return (PASS, f"chain {cid}: {got} residues") if got == want else (
        FAIL, f"chain {cid}: {got} residues, expected {want}")


def check_chain_nonstandard(fx, expect):
    """Non-standard residues must be counted. If they are dropped from the chain instead,
    the count reads 0 and the heavily-modified-chain warning can never fire."""
    cid, want = expect["chain"], expect["n_nonstandard"]
    chain = fx.chains.get(cid)
    if chain is None:
        return FAIL, f"chain {cid} absent (have {sorted(fx.chains)})"
    got = chain["n_nonstandard"]
    detail = f"chain {cid}: {got} non-standard of {chain['n_residues']}"
    return (PASS, detail) if got == want else (FAIL, detail + f", expected {want}")


def check_no_mixed_asym_id(fx, _):
    """No label_asym_id may hold both polymer and non-polymer atoms.

    RFD4 refuses such a file outright ("mixes defined and missing label_seq_id values"),
    so this is not a style point -- it decides whether the structure can be used at all.
    """
    if not os.path.exists(fx.out_cif):
        return FAIL, "no target.cif"
    idx, rows = atom_site(fx.out_cif)
    i_asym = idx.get("_atom_site.label_asym_id")
    i_seq = idx.get("_atom_site.label_seq_id")
    i_comp = idx.get("_atom_site.label_comp_id")
    if i_asym is None or i_seq is None:
        return FAIL, "target.cif has no label_asym_id/label_seq_id columns"
    tally = {}
    for row in rows:
        has_seq = row[i_seq] not in (".", "?", "")
        entry = tally.setdefault(row[i_asym], [0, 0, set()])
        entry[0 if has_seq else 1] += 1
        if not has_seq and i_comp is not None:
            entry[2].add(row[i_comp])
    mixed = [(a, v) for a, v in tally.items() if v[0] and v[1]]
    if mixed:
        asym, (withseq, without, comps) = mixed[0]
        return FAIL, (f"{len(mixed)} mixed; label_asym_id {asym} has {withseq} polymer "
                      f"and {without} non-polymer atoms {sorted(comps)[:4]}")
    return PASS, f"{len(tally)} asym id(s), none mixing polymer and non-polymer"


def check_no_warnings(fx, _):
    warnings = [l for l in fx.n1_log.splitlines() if l.startswith("WARNING")]
    if warnings:
        return FAIL, f"{len(warnings)} warning(s): {warnings[0][:70]}"
    return PASS, "clean run"


# --- node 04 checks: these read a node 04 run rather than a node 01 inventory ------

def load_metrics(case_dir):
    path = os.path.join(case_dir, "n4", "outputs", "design_metrics.json")
    if not os.path.exists(path):
        return None
    with open(path) as handle:
        return json.load(handle)


def check_node04_designed_chain(case_dir, expect):
    """Node 04 names the chain that was not in the target as the design. These inputs
    are built with a known designed chain, so the answer is checkable rather than
    plausible -- including when the chains are written in the other order."""
    payload = load_metrics(case_dir)
    if payload is None:
        return FAIL, "node 04 produced no design_metrics.json"
    designs = payload.get("designs") or []
    if not designs:
        return FAIL, "no designs in the metrics payload"
    row = designs[0]
    got_chain, got_len = row.get("designed_chain"), row.get("designed_length")
    want_chain, want_len = expect["chain"], expect["n_residues"]
    detail = (f"designed_chain={got_chain} designed_length={got_len} "
              f"(chain_lengths {row.get('chain_lengths')})")
    if got_chain == want_chain and got_len == want_len:
        return PASS, detail
    return FAIL, detail + f"; expected {want_chain}/{want_len}"


def check_design_rows_in_numeric_order(case_dir, _):
    """design_10 must not sort between design_1 and design_2."""
    import csv
    path = os.path.join(case_dir, "n4", "outputs", "design_metrics.csv")
    if not os.path.exists(path):
        return FAIL, "node 04 produced no design_metrics.csv"
    with open(path) as handle:
        names = [row["design"] for row in csv.DictReader(handle)]

    def index(name):
        digits = "".join(c for c in name.split("_")[1] if c.isdigit())
        return int(digits) if digits else -1

    got = [index(n) for n in names]
    if got == sorted(got):
        return PASS, f"{len(got)} rows in numeric order"
    shown = " ".join(n.replace("_generated.cif", "") for n in names[:6])
    return FAIL, f"row order is lexicographic: {shown} ..."


def check_node04_geometry(case_dir, expect):
    """Node 04 must measure the two geometric flags, and fire them where they belong.

    Nothing else in the pipeline would notice a broken model: contacts, net charge,
    composition and radius of gyration are all computable on an unfolded chain, and all
    look ordinary. This is the only check that would fail if the model itself were broken.

    The two flags are not equally strong, and the fixtures hold both readings. A chain
    break means the output is not a polypeptide. Not-compact means the ratio is past 1.8
    and nothing more: 'unfolded-design' is a straight chain, which deserves the flag, and
    'helical-design' is a textbook alpha-helix, which gets it too and is a real fold. Both
    are expected to be flagged, which is the point -- the flag is what is asserted here,
    not a verdict about what the design is.
    """
    payload = load_metrics(case_dir)
    if payload is None:
        return FAIL, "node 04 produced no design_metrics.json"
    designs = payload.get("designs") or []
    if not designs:
        return FAIL, "no designs in the metrics payload"
    row = designs[0]
    breaks, ratio = row.get("chain_breaks"), row.get("compactness_ratio")
    detail = (f"chain_breaks={breaks} Rg={row.get('radius_of_gyration_A')} "
              f"expected={row.get('expected_radius_of_gyration_A')} ratio={ratio}")
    log = ""
    log_path = os.path.join(case_dir, "n4.log")
    if os.path.exists(log_path):
        with open(log_path, errors="replace") as handle:
            log = handle.read()

    if "min_chain_breaks" in expect:
        if (breaks or 0) < expect["min_chain_breaks"]:
            return FAIL, detail + f"; expected at least {expect['min_chain_breaks']}"
        if "peptide bond spans" not in log:
            return FAIL, detail + "; breaks found but node 04 printed no warning"
    if expect.get("max_chain_breaks") is not None:
        if (breaks or 0) > expect["max_chain_breaks"]:
            return FAIL, detail + f"; expected at most {expect['max_chain_breaks']}"
    # Match the compactness warning on the two numbers it reports, not on its wording.
    # This used to look for the phrase "more extended than a folded protein", so
    # rewording the warning -- which is what narrowing the claim required -- broke the
    # check while the behaviour was unchanged. "<Rg> A vs <expected> A expected" is that
    # warning's own format and nothing else here prints it, so the assertion is pinned to
    # what node 04 measured rather than to how it phrased it.
    signature = (f"{row.get('radius_of_gyration_A')} A vs "
                 f"{row.get('expected_radius_of_gyration_A')} A expected")
    warned = any(line.startswith("WARNING") and signature in line
                 for line in log.splitlines())
    if expect.get("overextended") is True:
        if ratio is None or ratio <= 1.8:
            return FAIL, detail + "; expected a ratio above the 1.8 limit"
        if not warned:
            return FAIL, detail + "; overextended but node 04 printed no warning"
    if expect.get("overextended") is False:
        if ratio is not None and ratio > 1.8:
            return FAIL, detail + "; flagged as overextended but should not be"
        if warned:
            return FAIL, detail + "; node 04 warned about a compact design"
    return PASS, detail


def check_node05_clash_ranked_last(case_dir, expect):
    """A clashing design must not be presented first."""
    data, problem = report_payload(case_dir)
    if data is None:
        return FAIL, problem
    designs = (data.get("metrics") or {}).get("designs") or []
    order = [row.get("design", "?") for row in designs]
    if not order:
        return FAIL, "no design names found in the report"
    clashing = expect.get("clashing", "design_0_generated.cif")
    if clashing not in order:
        return FAIL, "%s absent from the report; got %s" % (clashing, order)
    if not designs[order.index(clashing)].get("steric_clash"):
        return FAIL, "%s is not flagged clashing, so its rank proves nothing" % clashing
    shown = [name.replace("_generated.cif", "") for name in order]
    if order[-1] != clashing:
        return FAIL, "clashing design is not last: %s" % shown
    return PASS, "order %s; clashing design last" % shown
def check_node05_reports_clash(case_dir, _):
    """The Clash column is rendered, and the clashing design carries the flag.

    Read from the payload, not the page text: node 05 names its columns in prose, so
    "Clash" appears in report.html whether or not the column is built.
    """
    data, problem = report_payload(case_dir)
    if data is None:
        return FAIL, problem
    rendered = [key for key, _label in data.get("columns", [])]
    if "steric_clash" not in rendered:
        return FAIL, "report renders no Clash column"
    designs = (data.get("metrics") or {}).get("designs") or []
    if not any(row.get("steric_clash") for row in designs):
        return FAIL, "report does not carry a steric_clash flag"
    return PASS, "Clash column rendered and a design is flagged"
def report_payload(case_dir):
    """The DATA object node 05 embeds, parsed.

    Checking the report by searching its text does not work: the page explains its own
    columns in prose, so every column label appears in the file whether or not a column
    is rendered. The table is built by inline JS from this object, so this is what a
    reader actually sees.
    """
    path = os.path.join(case_dir, "n5", "outputs", "report.html")
    if not os.path.exists(path):
        return None, "node 05 produced no report.html"
    with open(path, errors="replace") as handle:
        html = handle.read()
    marker = "const DATA = "
    if marker not in html:
        return None, "report embeds no DATA object"
    blob = html[html.index(marker) + len(marker):]
    cut = blob.find(chr(10) + "const TOP_N")
    if cut == -1:
        return None, "report's DATA object is not delimited as expected"
    try:
        return json.loads(blob[:cut].rstrip().rstrip(";")), ""
    except ValueError as exc:
        return None, f"report's DATA object is not valid JSON: {exc}"


def check_node05_reports_geometry(case_dir, expect):
    """The report must render node 04's flags, not merely carry them.

    chain_breaks and compactness_ratio were in the payload and in the CSV long before
    they were columns, and a design node 04 had flagged still arrived in the report
    as an ordinary row. So this reads the rendered column list rather than the page text,
    and reads the flag out of node 05's own copy of the metrics rather than node 04's
    file -- node04_geometry already tests node 04 on these same fixtures, and a check
    that re-read that file would pass on a report that showed none of it.
    """
    data, problem = report_payload(case_dir)
    if data is None:
        return FAIL, problem
    rendered = [key for key, _ in data.get("columns", [])]
    missing = [c for c in ("compactness_ratio", "chain_breaks", "overextended")
               if c not in rendered]
    if missing:
        return FAIL, f"report renders no column for {', '.join(missing)}"

    designs = (data.get("metrics") or {}).get("designs") or []
    if not designs:
        return FAIL, "report carries no design rows"
    row = designs[0] if len(designs) == 1 else None
    if row is None:
        flagged = [d for d in designs
                   if d.get("overextended") or d.get("chain_breaks")]
        if not flagged:
            return FAIL, "no design in the report carries a geometry flag"
        row = flagged[0]
    if expect.get("overextended") and not row.get("overextended"):
        return FAIL, "the report's own row does not carry node 04's overextended flag"
    if expect.get("min_chain_breaks") and not row.get("chain_breaks"):
        return FAIL, "the report's own row records no chain break"
    return PASS, (f"columns rendered; overextended={row.get('overextended')} "
                  f"breaks={row.get('chain_breaks')}")


def check_node05_unfolded_ranked_last(case_dir, expect):
    """A design node 04 flagged on geometry must sort below every unflagged one.

    The clash fixture cannot stand in for this. Clashes sort into a bucket of their own,
    so a report that sank clashes and nothing else would still pass every other ordering
    check in the suite -- which is exactly what the sink rule was added to change.

    Sinking is an ordering rule, not a verdict. The fixture is a straight chain, which has
    genuinely not folded, but the same sink catches an elongated real fold -- see the
    'helical-design' fixture. The designs are moved, not dropped, and the report says
    'Not compact' rather than naming them unfolded.
    """
    data, problem = report_payload(case_dir)
    if data is None:
        return FAIL, problem
    designs = (data.get("metrics") or {}).get("designs") or []
    if len(designs) < 2:
        return FAIL, "ordering needs two designs or more; report has %d" % len(designs)
    order = [row.get("design", "?") for row in designs]
    shown = [name.replace("_generated.cif", "") for name in order]
    named = expect.get("unfolded", "design_0_generated.cif")
    if named not in order:
        return FAIL, "%s absent from the report; got %s" % (named, order)

    def flagged(row):
        return bool(row.get("overextended")) or bool(row.get("chain_breaks"))

    if not flagged(designs[order.index(named)]):
        return FAIL, "%s carries no geometry flag, so its rank proves nothing" % named
    flags = [flagged(row) for row in designs]
    if all(flags):
        return FAIL, "every design is flagged, so the ordering proves nothing: %s" % shown
    # The property, rather than just "the named one is last": once the list reaches a
    # flagged design no ordinary one may follow it.
    first_flagged = flags.index(True)
    if False in flags[first_flagged:]:
        return FAIL, "an unflagged design is ranked below a flagged one: %s" % shown
    return PASS, "order %s; flagged design(s) sunk below the unflagged ones" % shown


def check_node06_per_design_fasta(case_dir, expect):
    """One FASTA per design, each holding exactly one record.

    The batch file going to workflow-026 as a single multi-record run is what conflates
    designs there, so the per-design files are the whole point of the export.
    """
    out = os.path.join(case_dir, "n6", "outputs")
    per = sorted(glob.glob(os.path.join(out, "refold_*.fasta")))
    want_n = expect["n_designs"]
    if len(per) != want_n:
        return FAIL, f"{len(per)} refold_*.fasta for {want_n} design(s)"
    if not os.path.exists(os.path.join(out, "sequences.fasta")):
        return FAIL, "no sequences.fasta batch file"
    want_seq = expect.get("sequence")
    for path in per:
        with open(path, errors="replace") as handle:
            lines = [l.strip() for l in handle if l.strip()]
        headers = [l for l in lines if l.startswith(">")]
        sequence = "".join(l for l in lines if not l.startswith(">"))
        if len(headers) != 1:
            return FAIL, f"{os.path.basename(path)} holds {len(headers)} records"
        if want_seq and sequence != want_seq:
            return FAIL, (f"{os.path.basename(path)} sequence {sequence!r} "
                          f"!= designed chain {want_seq!r}")
    return PASS, f"{len(per)} single-record file(s), each matching the designed chain"


NODE04_CHECKS = {
    "node04_designed_chain": check_node04_designed_chain,
    "design_rows_in_numeric_order": check_design_rows_in_numeric_order,
    "node04_geometry": check_node04_geometry,
    "node05_clash_ranked_last": check_node05_clash_ranked_last,
    "node05_reports_clash": check_node05_reports_clash,
    "node05_reports_geometry": check_node05_reports_geometry,
    "node05_unfolded_ranked_last": check_node05_unfolded_ranked_last,
    "node06_per_design_fasta": check_node06_per_design_fasta,
}


CHECKS = {
    "label_to_auth_complete": check_label_to_auth_complete,
    "inserted_residues_addressable": check_inserted_residues_addressable,
    "insertion_codes_warned_at_use": check_insertion_codes_warned_at_use,
    "log_mentions_models": check_log_mentions_models,
    "single_model_output": check_single_model_output,
    "nonpolymer_named": check_nonpolymer_named,
    "dropped_chains_named": check_dropped_chains_named,
    "chain_length": check_chain_length,
    "chain_nonstandard": check_chain_nonstandard,
    "no_mixed_asym_id": check_no_mixed_asym_id,
    "no_warnings": check_no_warnings,
}


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        return 2
    results, downloads = sys.argv[1], sys.argv[2]
    with open(os.path.join(HERE, "fixtures.json")) as handle:
        manifest = json.load(handle)

    rows, tally = [], {PASS: 0, FAIL: 0, GAP: 0, SKIP: 0}
    for spec in manifest["fixtures"]:
        fx = Fixture(results, downloads, spec)
        if fx.inventory is None:
            rows.append((fx.id, "inventory written", FAIL, "node 01 produced none"))
            tally[FAIL] += 1
            continue
        for expect in spec.get("expect", []):
            name = expect["check"]
            runner = CHECKS.get(name)
            if runner is None:
                status, detail = SKIP, "no such check"
            else:
                status, detail = runner(fx, expect)
            if status == FAIL and expect.get("status") == "open":
                status = GAP
            tag = expect.get("finding", "")
            rows.append((fx.id, f"{name} [{tag}]" if tag else name, status, detail))
            tally[status] += 1

    # --- node 02 must refuse these, and must not refuse the last one ---
    for case in manifest.get("rejections", []):
        slug = case["name"].replace(" ", "-")
        case_dir = os.path.join(results, "_rejections", slug)
        rc_path = os.path.join(case_dir, "rc")
        if not os.path.exists(rc_path):
            rows.append(("node02", case["name"], FAIL, "case was not run"))
            tally[FAIL] += 1
            continue
        with open(rc_path) as handle:
            rc = handle.read().strip()
        log = ""
        log_path = os.path.join(case_dir, "log")
        if os.path.exists(log_path):
            with open(log_path, errors="replace") as handle:
                log = handle.read()
        if case.get("expect_ok"):
            status = PASS if rc == "0" else FAIL
            detail = f"rc={rc}, accepted as it should be" if rc == "0" else f"rc={rc}, refused"
        else:
            fragment = case.get("expect_error", "")
            if rc == "0":
                status, detail = FAIL, "accepted; node 02 raised nothing"
            elif fragment and fragment not in log:
                status, detail = FAIL, f"refused (rc={rc}) but not for the stated reason"
            else:
                status, detail = PASS, f"refused (rc={rc}) with the expected reason"
        rows.append(("node02", case["name"], status, detail))
        tally[status] += 1

    # --- nodes 04, 05 and 06, on structures built by make_designs.py ---
    for case in manifest.get("downstream", []):
        case_dir = os.path.join(results, "_downstream", case["name"])
        for expect in case.get("expect", []):
            name = expect["check"]
            runner = NODE04_CHECKS.get(name)
            if runner is None:
                status, detail = SKIP, "no such check"
            elif not os.path.isdir(case_dir):
                status, detail = FAIL, "case was not run"
            else:
                status, detail = runner(case_dir, expect)
            if status == FAIL and expect.get("status") == "open":
                status = GAP
            tag = expect.get("finding", "")
            label = name.split("_")[0] if name.startswith("node0") else "node04"
            rows.append((label, f"{name} [{tag}]" if tag else name, status, detail))
            tally[status] += 1

    width = max(len(r[1]) for r in rows) + 2
    print(f"{'fixture':8s} {'check':{width}s} {'':5s} detail")
    print("-" * (20 + width + 60))
    for fixture, name, status, detail in rows:
        print(f"{fixture:8s} {name:{width}s} {status:5s} {detail}")
    print("-" * (20 + width + 60))
    print(f"{len(rows)} checks: {tally[PASS]} pass, {tally[FAIL]} fail, "
          f"{tally[GAP]} known gap, {tally[SKIP]} skipped")
    if tally[GAP]:
        print("Known gaps are tracked in fixtures.json; remove \"status\": \"open\" to enforce.")
    return 1 if tally[FAIL] else 0


if __name__ == "__main__":
    sys.exit(main())
