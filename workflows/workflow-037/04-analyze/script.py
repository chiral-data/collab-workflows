#!/usr/bin/env python3
"""Node 04 — Analyze Designs.

Parse each *_generated.cif into a per-design metrics row (chain lengths, designed-segment
length, composition, target contacts) and write design_metrics.csv / design_metrics.json.

Generation runs with `~metrics`, so RFD4 emits no scores of its own — every number here is
computed from the geometry of the output structures. There is no refolding oracle in this
workflow; self-consistency scoring is a follow-up that chains into workflow-026.
"""
import glob
import os
import sys

# TODO(scaffold): implement CIF parsing and metrics computation.


def main():
    os.makedirs("./outputs", exist_ok=True)

    designs = sorted(glob.glob("./inputs/*_generated.cif"))
    if not designs:
        print("ERROR: No generated structures found in ./inputs/ (*_generated.cif)", flush=True)
        sys.exit(1)

    print(f"Found {len(designs)} generated design(s)", flush=True)

    print("ERROR: Node 04 is not implemented yet (scaffold only).", flush=True)
    sys.exit(1)


if __name__ == "__main__":
    main()
