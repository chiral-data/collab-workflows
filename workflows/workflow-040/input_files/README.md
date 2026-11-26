# workflow-040 input files

This workflow's input is another workflow's output: the `*_generated.cif` files from
workflow-037 node 03, plus `design_metrics.json` from its node 04. There is no structure
to upload and no accession to fetch, so nothing is vendored here yet.

Two sets of files will land in this directory as the build proceeds, and both are listed
as blocked in [../tests/README.md](../tests/README.md) rather than left to memory:

**A real captured prediction.** One monomer prediction from `colabfold:2026_10_08` — the
`rank_001` scores JSON and the `rank_001` PDB — committed so node 04's parser is written
against a measured file rather than against a guess. `workflow-026/sample_outputs/confidence_summary.json`
is the reference for the shape today, from a real ColabFold 1.6.1 run, but it is 026's
post-processed summary rather than the raw file node 04 will read. The complex equivalent
has never been captured in this repo at all.

Workflow-037's `54dc7def` is the warning worth repeating here: it read the real package
before building and found that several of the documented guesses the workflow had been
built on were wrong.

**A few archived workflow-037 designs**, for control `C5` — the first honest ranking 037's
output has ever had. Either commit a handful of CIFs with their `design_metrics.json`, or
regenerate them on the GPU host.

A file earns its place in this directory by measurement. Probe it, and record what the
probe reported in `tests/fixtures.json` beside the entry that uses it.
