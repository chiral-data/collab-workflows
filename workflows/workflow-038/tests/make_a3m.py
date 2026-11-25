#!/usr/bin/env python3
"""Write one A3M per protein chain, so node 02 can run without the ColabFold server.

    python3 -I tests/make_a3m.py <entry.cif> <system_inventory.json> <out-dir> [rows]

Node 02's upload path takes `<chain>.a3m` from its inputs and contacts nothing. That is
what makes the offline suite possible: nodes 02 and 03 are exercised over real chain sets,
with real sequences, and the only thing not exercised is the server call -- which
`run_gpu_campaign.sh` does for real.

The query row is the chain's sequence **as the system states it**, which is what the model is
conditioned on and what node 03 now requires an uploaded alignment to match. That is not the
same string as the sequence in the coordinates: the entity sequence includes residues the
structure does not resolve (1CLL: 148 against 144) and carries `X` where a component is
non-standard (2ZTA's N-terminal acetyl). The measured sequence from probe_entry.py is printed
alongside it whenever the two differ, so the difference stays visible instead of being
discovered later as a mismatch.

The remaining rows are that sequence with scattered substitutions: enough depth (default 40
rows) to clear node 02's shallow-alignment warning, which is itself the point for the control
entry.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from probe_entry import probe  # noqa: E402

SUBSTITUTIONS = {"A": "G", "G": "A", "L": "V", "V": "L", "K": "R", "R": "K",
                 "D": "E", "E": "D", "S": "T", "T": "S", "F": "Y", "Y": "F"}


def main():
    if len(sys.argv) < 4:
        print(__doc__)
        return 2
    cif, inventory_path, out_dir = sys.argv[1:4]
    rows = int(sys.argv[4]) if len(sys.argv) > 4 else 40

    with open(inventory_path) as handle:
        inventory = json.load(handle)
    chains = list(inventory.get("protein_chains") or [])
    if not chains:
        print("no protein chain in the inventory; writing nothing")
        return 0

    measured = probe(cif)
    resolved = measured.get("sequences") or {}
    stated = {
        str(chain.get("asym_id")): chain.get("sequence")
        for chain in (inventory.get("chains") or [])
        if chain.get("sequence")
    }
    os.makedirs(out_dir, exist_ok=True)

    for chain in chains:
        sequence = stated.get(chain)
        if not sequence:
            # The inventory names a protein chain but states no sequence for it. That is a
            # finding, not something to paper over with a placeholder alignment.
            print(f"ERROR: no sequence stated for protein chain {chain}", file=sys.stderr)
            return 1
        in_file = resolved.get(chain)
        if in_file and in_file != sequence:
            print(
                f"  chain {chain}: system states {len(sequence)} residues, the coordinates "
                f"resolve {len(in_file)}"
                + (" (and the system carries X for a non-standard component)"
                   if "X" in sequence else "")
            )
        lines = [">query", sequence]
        for index in range(1, rows):
            mutated = list(sequence)
            for position in range(index % 7, len(mutated), 11):
                mutated[position] = SUBSTITUTIONS.get(mutated[position], mutated[position])
            lines += [f">synthetic_{index:03d}", "".join(mutated)]
        path = os.path.join(out_dir, f"{chain}.a3m")
        with open(path, "w") as handle:
            handle.write("\n".join(lines) + "\n")
        print(f"  wrote {path}: {rows} rows, {len(sequence)} residues")
    return 0


if __name__ == "__main__":
    sys.exit(main())
