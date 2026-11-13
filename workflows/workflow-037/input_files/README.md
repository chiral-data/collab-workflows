# Sample input files

## 5o45_cropped.cif.zst

A cropped copy of PDB [5O45](https://www.rcsb.org/structure/5O45) (PD-L1), the target used by
the RFD4-Proteina binder quickstart and by this workflow's correctness gate.

- **Source:** `docs/examples/inputs/5o45_cropped.cif.zst` in
  [RosettaCommons/RFD4-Proteina](https://github.com/RosettaCommons/RFD4-Proteina)
- **Licence:** BSD-3-Clause, same as RFD4-Proteina itself
- **Retrieved:** 2026-11-13, 27,393 bytes

Upstream's own binder spec for this target uses `contig: "50-120, /0, A17-131"`, i.e. a 50 to 120
residue binder against chain A. This workflow derives the target segment from the structure itself
rather than hardcoding it, so node 02 should produce the same span.

To use it, place it in the workflow's input directory and leave `pdb_id` empty; an uploaded
structure always takes precedence over the accession. Any `.cif`, `.pdb` or `.cif.zst` file works.

Note that fetching `pdb_id=5o45` from RCSB instead gives the **full** entry, not this crop: chain A
runs to label 129 (author 17 to 145) and a second chain B is present. Both are valid targets; the
crop is simply what upstream's documented numbers refer to.
