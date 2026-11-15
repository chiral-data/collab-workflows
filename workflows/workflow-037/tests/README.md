# workflow-037 fixture suite

Ten real RCSB entries, run through nodes 01 and 02, with the output checked against a
measurement of the input.

## Why this exists

Every one of these ten structures exits 0 on both nodes. Before this suite, that was the
whole test: the pipeline ran, so the pipeline worked.

It did not. The first run of these checks found twelve things wrong, none of which moved
an exit code, because the failure mode for a bad input shape is silence — a residue that
cannot be addressed, an NMR ensemble handed to the model whole, a cofactor deleted
without mention. Exit codes cannot see any of that. Comparing what a node *reported*
against what the file *contains* can.

That comparison is the entire design. `probe_entry.py` parses `_atom_site` by hand rather
than through Biopython, so the measurement does not share a parser — or a bug — with the
code under test.

## Running it

```sh
./tests/run_fixtures.sh              # writes to ./.fixture-run
./tests/run_fixtures.sh /some/dir    # or somewhere else
```

Needs `docker`, `curl`, `python3`, and the image from `../Dockerfile`
(override with `RFD4_IMAGE=...`). Structures are fetched once and cached.
Exit status is 0 when nothing in the `FAIL` column; known gaps do not fail the run.

Nodes 03–06 are out of scope: 03 needs a GPU and minutes per design, and every bug this
suite is built for bites in 01 or 02.

The pieces run on their own too:

```sh
python3 -I tests/probe_entry.py .fixture-run/downloads       # measure the inputs
python3 -I tests/check_fixtures.py .fixture-run/results .fixture-run/downloads
```

Downloaded structures are untrusted input. They stay in their own directory, nothing is
executed from it, and every `python3` here runs with `-I` so a file dropped next to a
script cannot be imported in place of a stdlib module.

## Choosing a fixture

A fixture earns its place by **measurement**, never by reputation. Each entry in
`fixtures.json` carries an `evidence` string recording what `probe_entry.py` reported
from the downloaded file. Re-run the probe before trusting it; RCSB re-releases entries.

This rule has already paid once. The suite was going to carry a NaN-coordinate fixture,
on the strength of a finding that said node 02 should guard against non-finite
coordinates. Probing twelve candidate entries found **zero** NaN coordinates — real mmCIF
files do not write them. The shape that actually occurs is *occupancy 0.0 with
real-looking coordinates*, present in four of the twelve (12E8: 80 atoms, 5O45: 44,
1A3N: 37, 1IGT: 18). The fixture set followed the measurement.

## What each fixture is for

| Entry | Shape | Catches |
|---|---|---|
| 12E8 | Fab, Kabat insertion codes 52A / 82A-C / 100A-C | residues in CDR-H2 and CDR-H3 that could not be addressed |
| 1IGT | intact IgG, insertion codes + N-glycans | same, plus branched non-polymer |
| 1G03 | NMR, 20 models | ensemble passed through whole |
| 2KOX | NMR, 640 models | same, at 787,840 atom rows for a 76-residue protein |
| 1ZNF | NMR + structural Zn + capped termini | three shapes at once |
| 1A3N | haemoglobin, 4 chains, haem, altlocs A/B | cofactor at the interface; chains dropped by the contig |
| 1BRS | barnase–barstar, 6 chains | five of six chains deleted silently |
| 5O45 | macrocyclic peptide, 5 modified components | HETATM filtering shortening a chain; altloc past B |
| 1UBQ | ubiquitin, nothing unusual | control — a failure here is the harness |
| 4ZQK | PD-1/PD-L1 | regression: the README's worked example |

## Known gaps

Expectations marked `"status": "open"` in `fixtures.json` report as `GAP` and do not fail
the run. They are real, reproduced, and not yet fixed:

- **#4** — non-polymer entities (HEM, NAG, NA) are dropped without being named. A haem at
  a chain interface changes what the interface is.
- **#10** — chains the contig does not preserve are deleted without being named. 1BRS
  loses five of six.

Fixing one means deleting its `status` field. The suite then enforces it, and a
regression shows up as `FAIL` rather than as silence.

## Adding a fixture

1. Download it and run `probe_entry.py` against it.
2. Add it to `fixtures.json` with the measured `evidence`, the `params` to run it under,
   and a `why` that says what a reader would get wrong without it.
3. Add `expect` entries naming existing checks. A new kind of check goes in the `CHECKS`
   table in `check_fixtures.py` and returns `(status, detail)` — the detail string should
   carry the numbers, so a failure reads as evidence rather than as an assertion.
