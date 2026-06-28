"""Node 05 — Report.

Applies filters, ranks survivors, generates scrambled negative controls,
computes success rate, writes HTML report + summary JSON + ranked CSV.
"""
from __future__ import annotations
import json
import pathlib

PARAMS   = json.loads(pathlib.Path("inputs/global_params.json").read_text())
SCORES   = json.loads(pathlib.Path("inputs/scores.json").read_text())
MANIFEST = json.loads(pathlib.Path("inputs/manifest.json").read_text())

FILTER_IPTM    = PARAMS["filter_iptm_min"]
FILTER_PLDDT   = PARAMS["filter_binder_plddt_min"]
FILTER_RMSD    = PARAMS["filter_self_consistency_rmsd_max"]

OUT = pathlib.Path("outputs")

# TODO: implement apply_filters, rank, make_scrambled_controls,
#       compute_success_rate, write_html_report, write_summary_json, write_csv

def main() -> None:
    raise NotImplementedError("generate_report.py is not yet implemented")

if __name__ == "__main__":
    main()
