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
from probe_entry import probe  # noqa: E402

PASS, FAIL, GAP, SKIP = "PASS", "FAIL", "GAP", "SKIP"


class Fixture:
    """Everything measured or read for one entry, loaded once."""

    def __init__(self, results, downloads, spec):
        self.id = spec["id"]
        self.spec = spec
        base = os.path.join(results, self.id)
        self.inv_path = os.path.join(base, "n1", "outputs", "target_inventory.json")
        self.out_cif = os.path.join(base, "n1", "outputs", "target.cif")
        self.src_cif = os.path.join(downloads, self.id + ".cif")
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
    want = expect.get("components", [])
    present = [c for c in want if c in fx.truth.get("nonpolymer_list", [])]
    if not present:
        return FAIL, f"{want} not in the input at all; fixture no longer fits"
    named = [c for c in present if c in fx.n1_log]
    if not named:
        return FAIL, f"{present} present in input, never named in node 01 output"
    return PASS, f"named: {named}"


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


def check_no_warnings(fx, _):
    warnings = [l for l in fx.n1_log.splitlines() if l.startswith("WARNING")]
    if warnings:
        return FAIL, f"{len(warnings)} warning(s): {warnings[0][:70]}"
    return PASS, "clean run"


CHECKS = {
    "label_to_auth_complete": check_label_to_auth_complete,
    "inserted_residues_addressable": check_inserted_residues_addressable,
    "log_mentions_models": check_log_mentions_models,
    "single_model_output": check_single_model_output,
    "nonpolymer_named": check_nonpolymer_named,
    "dropped_chains_named": check_dropped_chains_named,
    "chain_length": check_chain_length,
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
