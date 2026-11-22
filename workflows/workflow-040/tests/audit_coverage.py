#!/usr/bin/env python3
"""Which of node 04's columns does no check in the suite ever assert?

Why this exists
---------------
plddt_binder_in_complex was rendered by node 05, carried in design_scores.json and
design_scores.csv, and asserted by nothing. It appeared in no check and nowhere in
fixtures.json. A new metric with no check is the thing this suite exists to prevent, and
it got in anyway during the Phase 2 build -- which means noticing it cannot be left to
whoever happens to read the diff.

So this asks the question mechanically: for every name in node 04's COLUMNS, does any
check function or any fixture expectation mention it?

What a finding here does and does not mean
------------------------------------------
A column with no mention is not automatically a defect. Several columns are provenance
rather than measurement -- the thresholds a row was judged against, the design's own name
-- and are asserted indirectly by checks that compare whole rows or whole tables. So the
output is a LIST TO JUSTIFY, not a verdict, and each name either gets a check or gets an
entry in ASSERTED_INDIRECTLY saying what covers it. The point is that the answer is
written down somewhere, not that it is always "add a check".

Usage:  python3 tests/audit_coverage.py [--verbose]
"""
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
WF = os.path.dirname(HERE)

# Columns no check names directly, with what covers them instead. Anything not here and
# not mentioned by a check is a finding -- so this table cannot quietly absorb the next
# unchecked metric.
ASSERTED_INDIRECTLY = {
    "design": "the key every check looks a row up by; a wrong one means no row at all",
    "status": "read by design_not_passing, design_recorded_failed and the report checks",
    "reason": "report_row_states_reason requires it to reach the rendered row",
    "binder_sequence": "passing_fasta_tagged and tag_has_no_colon read the FASTA it writes",
    "sc_rmsd_n_residues": "rmsd_mapping_asserted fails when it disagrees with the mapping",
    "failed_thresholds": "design_not_passing requires the failing metric to be named here",
    "max_sc_rmsd_a": "thresholds_in_output requires every gating threshold on the row",
    "min_plddt_binder": "thresholds_in_output, same",
    "min_iptm": "thresholds_in_output, and metric_does_not_gate for the opt-in case",
    "max_pae_interaction": "thresholds_in_output, same",
    "interface_scored": "metric_does_not_gate reads it to tell Phase 1 from Phase 2",
    "interface_gated": "metric_does_not_gate: this is the flag it asserts is false",
    "interface_reason": "null_reason_recorded covers the reason strings as a family",
    "cplx_msa_mode": "argv_includes asserts the mode on the argv, which is the live value",
    "cplx_binder_aligned": "no_network_reached and the 01b outputs carry the same claim",
    "ptm": "metric_is_null and metric_not_value take it by name from the fixture",
    "max_pae": "same, by name from the fixture",
    "n_unknown": "unknown_residues_counted, which probes the sequence independently",
    "gen_compactness_ratio": "metrics_merged_by_name asserts the 037 block merges by name",
    "gen_overextended": "metrics_merged_by_name",
    "gen_chain_breaks": "metrics_merged_by_name",
    "gen_steric_clash": "metrics_merged_by_name",
    "gen_n_interface_contacts": "metrics_merged_by_name",
}


def columns():
    source = open(os.path.join(WF, "04-score", "script.py")).read()
    block = re.search(r"^COLUMNS\s*=\s*\[(.*?)^\]", source, re.S | re.M)
    if not block:
        raise SystemExit("audit_coverage: could not find COLUMNS in 04-score/script.py")
    return re.findall(r'"([^"]+)"', block.group(1))


def mentioned():
    """Every column name any check function or any fixture expectation refers to."""
    checks = open(os.path.join(HERE, "check_fixtures.py")).read()
    fixtures = json.load(open(os.path.join(HERE, "fixtures.json")))
    text = checks + json.dumps(fixtures)
    return text


def main(argv):
    verbose = "--verbose" in argv
    text = mentioned()
    unmentioned, justified = [], []
    for column in columns():
        # Word-boundary, so plddt_binder does not count as covering
        # plddt_binder_in_complex -- which is the exact pair that went unchecked.
        if re.search(r'"%s"' % re.escape(column), text):
            if verbose:
                print(f"  checked    {column}")
            continue
        if column in ASSERTED_INDIRECTLY:
            justified.append(column)
            if verbose:
                print(f"  indirect   {column}: {ASSERTED_INDIRECTLY[column]}")
            continue
        unmentioned.append(column)

    print(f"\n{len(columns())} column(s) in node 04's COLUMNS")
    print(f"  named by a check or a fixture: "
          f"{len(columns()) - len(justified) - len(unmentioned)}")
    print(f"  asserted indirectly, justified in this file: {len(justified)}")
    print("=" * 78)
    if unmentioned:
        print(f"{len(unmentioned)} column(s) that NO check and NO fixture mentions:")
        for column in unmentioned:
            print(f"  - {column}")
        print("\nEach one needs a check, or an entry in ASSERTED_INDIRECTLY saying what")
        print("covers it. plddt_binder_in_complex sat in this list unnoticed through a")
        print("whole phase of work, which is why the answer has to be written down.")
        print("=" * 78)
        return 1
    print("Every column is either named by a check or justified as indirect.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
