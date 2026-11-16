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
import json
import os
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
    said = any(word in fx.n2_log for word in ("not preserved", "dropped", "discard"))
    detail = f"{len(dropped)} dropped {dropped}; node 02 names them: {said}"
    return (PASS, detail) if said else (FAIL, detail)


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
    path = os.path.join(case_dir, "outputs", "design_metrics.json")
    if not os.path.exists(path):
        return None
    with open(path) as handle:
        return json.load(handle)


def check_node04_designed_chain(case_dir, expect):
    """Node 04 names the first chain of the file as the design. These inputs are built
    with a known designed chain, so the answer is checkable rather than plausible."""
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
    path = os.path.join(case_dir, "outputs", "design_metrics.csv")
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


NODE04_CHECKS = {
    "node04_designed_chain": check_node04_designed_chain,
    "design_rows_in_numeric_order": check_design_rows_in_numeric_order,
}


CHECKS = {
    "label_to_auth_complete": check_label_to_auth_complete,
    "inserted_residues_addressable": check_inserted_residues_addressable,
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

    # --- node 04, on structures built by make_designs.py ---
    for case in manifest.get("node04", []):
        case_dir = os.path.join(results, "_node04", case["name"])
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
            rows.append(("node04", f"{name} [{tag}]" if tag else name, status, detail))
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
