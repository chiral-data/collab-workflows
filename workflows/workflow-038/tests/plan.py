#!/usr/bin/env python3
"""Turn fixtures.json into the run plan, so the driver script holds no second copy of it.

    python3 -I tests/plan.py entries          # the structures to download
    python3 -I tests/plan.py fixtures         # slug  id
    python3 -I tests/plan.py rejections       # slug  id  stage  setup  env...
    python3 -I tests/plan.py downstream       # slug  source  n_frames  shape  physicality  setup  env...

Tab separated, one case per line. A case added to fixtures.json is run by the driver with
no change to the driver -- which is the point: the expectations and the run cannot drift
apart, because they are the same file.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))


def slug(name):
    return "".join(c if c.isalnum() or c in "-_." else "-" for c in name).strip("-")


def env(entry):
    return [f"{k}={v}" for k, v in sorted((entry.get("params") or {}).items())]


def main():
    section = sys.argv[1] if len(sys.argv) > 1 else ""
    with open(os.path.join(HERE, "fixtures.json")) as handle:
        spec = json.load(handle)

    if section == "entries":
        ids = {e["id"] for e in spec["fixtures"]}
        ids |= {e["id"] for e in spec["rejections"] if e.get("id")}
        ids |= {e["source"] for e in spec["downstream"] if e.get("source")}
        print("\n".join(sorted(ids)))
        return 0

    if section == "fixtures":
        for entry in spec["fixtures"]:
            print("\t".join([slug(entry["id"]), entry["id"]]))
        return 0

    if section == "rejections":
        for entry in spec["rejections"]:
            fields = [
                slug(entry["name"]),
                entry.get("id", ""),
                entry.get("stage", ""),
                entry.get("setup") or "none",
            ]
            print("\t".join(fields + env(entry)))
        return 0

    if section == "downstream":
        for entry in spec["downstream"]:
            fields = [
                slug(entry["name"]),
                entry.get("source", ""),
                str(entry.get("n_frames", 6)),
                entry.get("ensemble", "spread"),
                entry.get("physicality", "passed"),
                entry.get("setup") or "none",
            ]
            print("\t".join(fields + env(entry)))
        return 0

    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main())
