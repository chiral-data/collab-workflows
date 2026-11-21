#!/usr/bin/env python3
"""Write C12's provenance record, measured off the generated file rather than asserted.

C12 is the one control not built by make_controls.py -- it is workflow-037's own
make_designs.py output, an idealised helix carrying a designed sequence -- so its
provenance has to be written separately. It is written here, with the same reader and the
same field names the other four get, because the first version of this was a hand-written
dict with four keys and report_controls.py raised a KeyError on the fields it did not
have. That landed after the table and before the verdict, so the tee'd VERDICT.txt of
that run is a file that reads like a complete report and contains no verdict at all
(spec §8.3: the only sign was an exit code, and the artifact looked fine).

Two things follow, and both are implemented:
  - a record is MEASURED off the artifact, never asserted about it (§8.14);
  - and report_controls.py reports an absent field instead of dying on it, and refuses
    the gate when a gating control has no provenance.

Run:  python3 -I tests/make_c12_provenance.py <C12_generated.cif> <out.json>
"""
import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_controls  # noqa: E402  the same reader the other four controls were built with


def described(entry):
    """The field set make_controls.py records, so one printer reads every record."""
    return {
        "auth_asym_id": entry["auth"],
        "label_asym_id_written": entry["auth"],
        "n_residues": entry["n_residues"],
        "first_id": entry["first_id"],
        "last_id": entry["last_id"],
        "gaps": entry["gaps"],
        "sequence": entry["sequence"],
        "n_unknown": entry["n_unknown"],
    }


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    cif, out = sys.argv[1], sys.argv[2]

    chains, _ = make_controls.chains_in(cif)
    by_auth = {entry["auth"]: entry for entry in chains}
    if make_controls.BINDER_LABEL not in by_auth:
        print(f"make_c12_provenance: no chain {make_controls.BINDER_LABEL} in {cif}; "
              f"it holds {sorted(by_auth)}", file=sys.stderr)
        return 1
    binder = by_auth[make_controls.BINDER_LABEL]

    record = {
        "name": "C12",
        "source_file": "workflow-037/tests/make_designs.py shape=helix length=60",
        "source_sha256": hashlib.sha256(open(cif, "rb").read()).hexdigest(),
        "source_note": "generated, not deposited: the digest is of the generated file",
        "sequence_transform": "designed",
        "transform_note": (
            "a textbook alpha helix carrying make_designs.py's ALEKIAELKA design cycle. "
            "Synthetic, and that is the point: it is the class of input this workflow is "
            "for, where single_sequence prediction is the correct mode (spec §4). It "
            "is also the control that tells a mis-specified positive control apart from a "
            "broken scorer, which C1 to C4 could not do between them."
        ),
        "binder": described(binder),
        "targets": [described(by_auth[label])
                    for label in make_controls.TARGET_LABELS if label in by_auth],
    }
    with open(out, "w") as handle:
        json.dump(record, handle, indent=2)

    binder_record = record["binder"]
    print(f"C12 provenance: binder {binder_record['auth_asym_id']}, "
          f"{binder_record['n_residues']} residues, ids {binder_record['first_id']}.."
          f"{binder_record['last_id']}, {len(binder_record['gaps'])} gap(s), "
          f"{len(record['targets'])} target(s) -> {out}")
    print(f"  sequence {binder_record['sequence']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
