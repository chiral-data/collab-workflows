# Input files for workflow-039

## example_panel.csv

A two-clone smoke-test panel. Both entries are the paired VH/VL examples from the upstream
tools' own documentation, so they are known to be numberable by ANARCI and modellable by
ABodyBuilder2:

| clone_id | VH length | VL length | Source |
|----------|-----------|-----------|--------|
| `ib_example` | 121 | 107 | [ImmuneBuilder README](https://github.com/brennanaba/ImmuneBuilder) usage example |
| `ablang2_example` | 121 | 112 | [AbLang2 README](https://github.com/oxpig/AbLang2) usage example |

These are fixtures for exercising the pipeline, **not** a benchmark panel: no clinical-stage
reference values are attached to them, so they tell you the nodes run, not that the metrics are
right. See the pipeline doc's verification plan for what would establish that.

## Panel format

Either of these, selected by extension:

**CSV** — a `clone_id` column plus `vh` and `vl` columns:

```csv
clone_id,vh,vl
clone_01,EVQL...VTVSS,DIQM...KVEIK
```

For nanobodies, give `vh` (or `vhh`) and leave `vl` empty. Unknown extra columns are carried
through to the ranked output untouched, so target name, assay EC50 or library round survive into
the report.

**FASTA** — one record per chain, with the chain suffixed on the record id:

```text
>clone_01_VH
EVQL...VTVSS
>clone_01_VL
DIQM...KVEIK
```

Accepted suffixes are `_VH`/`_H`/`_HEAVY` and `_VL`/`_L`/`_LIGHT`/`_K`/`_KAPPA`/`_LAMBDA`,
case-insensitive. A record with no recognised suffix is an error rather than a guess.

## Sequence requirements

- **Variable domains only.** Give the whole VH or VL domain; ImmuneBuilder asserts the sequence
  is longer than 70 residues and ANARCI must recognise it as the chain type claimed. Fabs, full
  IgGs, scFv linkers and CDR-only fragments all fail node 01 — trim to the variable domain first.
- **One-letter amino acid codes.** `X` and other ambiguity codes are rejected by node 01 rather
  than silently modelled.
- **An scFv** must be split into its VH and VL at the linker before it is submitted. Node 01 does
  not detect or cut linkers.

## Verifying node 04 on its own

promb's documented example gives a fixed expected value, which makes it a useful independent
check of the humanness node:

```python
from promb import init_db
db = init_db('human-oas')
db.compute_peptide_content('QVQLVQSGVEVKKPGASVKVSCKASGYTFTNYYMYWVRQAPGQGLEWMGGINPSNGGTNFNEKFKNRVTLTTDSSTTTAYMELKSLQFDDTAVYYCARRDYRFDMGFDYWGQGTTVTVSS')
# 0.7142857142857143
```

If node 04 reports something other than ≈0.714 for that VH, the node is wrong, not the database.

## No bundled nanobody

The `nanobody` path (NanoBodyBuilder2) has no bundled example, because no VHH sequence here has
been verified against upstream documentation. Supply your own VHH to exercise it, and treat that
path as untested until you do.
