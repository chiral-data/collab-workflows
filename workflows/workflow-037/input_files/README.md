# Sample input files

## 5o45_cropped.cif.zst

**Not yet vendored.** The RFD4-Proteina documentation quickstart uses a cropped copy of PDB
[5O45](https://www.rcsb.org/structure/5O45) (PD-L1) as its binder-design target. Place that file
here as `5o45_cropped.cif.zst` before running the correctness gate described in
[rfd4_pipeline.md](../rfd4_pipeline.md).

- **Source:** RFD4-Proteina documentation examples (RosettaCommons / NVIDIA GenAIR)
- **Licence:** BSD-3-Clause, same as RFD4-Proteina itself
- **Attribution:** retained here as the upstream sample; not redistributed under any other terms

Any `.cif`, `.pdb` or `.cif.zst` file dropped into the workflow's `inputs/` directory at run time
takes precedence over the `pdb_id` parameter.
