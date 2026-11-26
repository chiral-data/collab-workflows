#!/usr/bin/env python3
"""Build a stand-in target alignment for node 01b, from an inventory node 01 really wrote.

Why this exists
---------------
Node 01b aligns the TARGET chain, from a file or from the MSA server. The suite does not
touch the network, so the file arm is the one it can exercise -- and the file has to be a
real alignment of the real target, because 01b validates the query row against the target
sequence it derived from the inventory (`01b-align-target/script.py` `from_file`). A
hard-coded a3m would mismatch and 01b would correctly refuse it, which would test the
refusal instead of the path.

So the query comes OUT of `design_inventory.json` and the homologs are built from it, one
point substitution each. They are synthetic on purpose: nothing is sent anywhere to build
a test fixture, and a fixture whose provenance is a web request is not a fixture.

What it deliberately does NOT do
--------------------------------
It does not refuse anything. An inventory with no target chain (a monomer), or one whose
designs disagree about their target, is exactly what `monomer-on-complex-path` and §4's
one-target-per-run rule are about -- and those refusals belong to 01b. If this tool
refused them first, the suite would be asserting its own harness and the node's guard
would never be shown refusing (§8.13). In those cases it writes nothing, says why on
stdout, and exits 0 so the node is reached.

Usage
-----
  make_target_msa.py <design_inventory.json> <out.a3m> [n_homologs]

Default 4 homologs, which clears PARAM_MIN_MSA_SEQUENCES (2) with room to spare so a
change to that default does not silently turn these cases into the refusal case.
"""
import json
import sys

# A substitution has to land inside the alphabet node 01b checks the query against, and
# has to be a different residue from the one it replaces, or the "homolog" is the query.
SWAP = {"A": "G"}
DEFAULT_SWAP = "A"


def target_of(inventory):
    """(sequence, note). The fused target every live design agrees on, or None."""
    live = [row for row in inventory.get("designs") or [] if not row.get("skipped")]
    if not live:
        return None, "every design in the inventory is skipped"

    seen = {}
    for row in live:
        chains = row.get("target_chains") or []
        key = tuple((c["chain"], c["sequence"]) for c in chains)
        seen.setdefault(key, []).append(row.get("design"))

    if len(seen) > 1:
        return None, (
            f"the {len(live)} live design(s) declare {len(seen)} different targets, which "
            "is 01b's refusal to make, not this tool's"
        )

    key = next(iter(seen))
    if not key:
        return None, (
            "no design declares a target chain, so there is no target to align. On the "
            "complex path that is 01b's refusal to make, not this tool's"
        )
    # Concatenated in the inventory's own chain order, which is what 01b fuses (§8.12).
    return "".join(sequence for _, sequence in key), None


def homologs(sequence, count):
    """`count` rows, each the query with one residue substituted, spread along the chain."""
    rows = []
    for index in range(count):
        at = (index * 7) % len(sequence)
        here = sequence[at]
        rows.append((f"homolog{index + 1}",
                     sequence[:at] + SWAP.get(here, DEFAULT_SWAP) + sequence[at + 1:]))
    return rows


def main(argv):
    if len(argv) < 3:
        print(__doc__.strip().splitlines()[-4], file=sys.stderr)
        return 2
    inventory_path, out_path = argv[1], argv[2]
    count = int(argv[3]) if len(argv) > 3 else 4

    with open(inventory_path) as handle:
        inventory = json.load(handle)

    sequence, why = target_of(inventory)
    if sequence is None:
        print(f"make_target_msa: no alignment written -- {why}.")
        return 0

    rows = [("query", sequence)] + homologs(sequence, count)
    # No '#' header: 01b writes the paired header itself downstream, and a single-chain
    # a3m carrying one would be refused by node 02's load_target. Match columns only, no
    # lowercase insertions, so every row lines up with the query column for column.
    with open(out_path, "w") as handle:
        for name, row in rows:
            handle.write(f">{name}\n{row}\n")

    print(f"make_target_msa: {out_path} -- {len(rows)} row(s) "
          f"({len(rows) - 1} synthetic homolog(s)) for a {len(sequence)}-residue target.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
