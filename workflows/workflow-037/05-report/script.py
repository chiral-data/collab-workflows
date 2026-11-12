#!/usr/bin/env python3
"""Node 05 — Report.

Build a single self-contained report.html: Bootstrap 5 + Plotly + Mol* from CDN, with the
metrics table and the top_n structures inlined as JavaScript so the file works standalone.

The report must carry the biosafety notice prominently: the RFD4 CLI and Python API do not
run SecureDNA screening (only the HF Spaces apps do), so designs leaving this workflow are
unscreened and require review before synthesis.
"""
import glob
import os
import sys

# TODO(scaffold): implement HTML rendering with inlined data and Mol* viewer.


def main():
    os.makedirs("./outputs", exist_ok=True)

    top_n = int(os.environ.get("PARAM_TOP_N", "10"))

    metrics_path = "./inputs/design_metrics.json"
    if not os.path.exists(metrics_path):
        print(f"ERROR: Metrics not found: {metrics_path}", flush=True)
        sys.exit(1)

    designs = sorted(glob.glob("./inputs/*_generated.cif"))
    if not designs:
        print("ERROR: No generated structures found in ./inputs/ (*_generated.cif)", flush=True)
        sys.exit(1)

    print(f"Reporting {min(top_n, len(designs))} of {len(designs)} design(s)", flush=True)

    print("ERROR: Node 05 is not implemented yet (scaffold only).", flush=True)
    sys.exit(1)


if __name__ == "__main__":
    main()
