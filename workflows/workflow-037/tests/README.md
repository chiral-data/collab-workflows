# workflow-037 fixture suite

Ten real RCSB entries run through nodes 01 and 02, the inputs node 02 must refuse, and
nodes 04, 05 and 06 driven by structures built for the purpose -- all checked against a
measurement of the input rather than against an exit code.

## Why this exists

Every one of these structures exits 0 on nodes 01 and 02. Before this suite, that was
the whole test: the pipeline ran, so the pipeline worked.

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

It runs three sections:

- **fixtures** -- nodes 01 and 02 over real structures, output checked against the input.
- **rejections** -- inputs node 02 must refuse, asserting both the exit status and that
  the refusal gives the stated reason. One entry checks a legitimate override is still
  accepted, so a guard cannot pass by rejecting everything.
- **downstream** -- nodes 04, 05 and 06 in sequence, over structures from
  `make_designs.py`.

Needs `docker`, `curl`, `python3`, and the image from `../Dockerfile`
(override with `RFD4_IMAGE=...`). Structures are fetched once and cached.
Exit status is 0 when nothing in the `FAIL` column; known gaps do not fail the run.

Node 03 is out of scope: it needs a GPU and minutes per design. Everything after it is
reached anyway, because nodes 04, 05 and 06 read only chain ids, residue names, sequence
numbers and CA coordinates, so `make_designs.py` builds two-chain structures directly.
That exercises them exactly as a real design does, and deterministically, which a sampled
design never could.

## What this does not establish

The suite tests the plumbing, not the science. Passing means the pipeline does not
silently corrupt, drop, mislabel or misnumber its input, and that it notices output which
is not a folded polypeptide. It does not mean a design is good.

Nothing here refolds a sequence, computes self-consistency, or compares a result against
a published one. The fixtures are real entries chosen for the input shapes they carry,
not because the workflow is known to reproduce a literature result on them. There is no
positive control in the sense of a target with a known-good designed binder.

The `model-sanity` cases are the closest thing to one. A broken checkpoint would leave
every other check passing -- contacts, net charge, composition and radius of gyration are
all computable on an unfolded chain and all look ordinary -- so node 04 measures CA-CA
spacing and compares Rg against 2.2*N^0.38 A, and the suite feeds it a folded chain, an
unfolded one and a broken one to prove the check fires where it should and stays quiet
where it should. That catches a model that has stopped working. It will not catch a model
that works and produces designs which do not bind.

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
| 12E8-hotspots | the same Fab, hotspots named over CDR-H3 in author numbering | insertion-coded residues left out of a hotspot span with nothing said where it happens |
| 1IGT | intact IgG, insertion codes + N-glycans | same, plus branched non-polymer |
| 1G03 | NMR, 20 models | ensemble passed through whole |
| 2KOX | NMR, 640 models | same, at 787,840 atom rows for a 76-residue protein |
| 1ZNF | NMR + structural Zn + capped termini | three shapes at once |
| 1A3N | haemoglobin, 4 chains, haem, altlocs A/B | cofactor at the interface; chains dropped by the contig |
| 1BRS | barnase–barstar, 6 chains | five of six chains deleted silently |
| 5O45 | macrocyclic peptide, 5 modified components | HETATM filtering shortening a chain; altloc past B |
| 1UBQ | ubiquitin, nothing unusual | control — a failure here is the harness |
| 5O45-pdb | the same entry via the PDB conversion path | modified residues deleted from the chain; non-polymers sharing the polymer's asym id |
| 4ZQK | PD-1/PD-L1 | regression: the README's worked example |

## Known gaps

Expectations marked `"status": "open"` in `fixtures.json` report as `GAP` and do not fail
the run. They are real, reproduced, and not yet fixed.

**Nothing is marked open today.** Two findings were tracked here and both are closed:

- **#27** — node 04 identified the design as the first chain in the file. It now reads
  the target's chain ids out of `generation_manifest.json` and takes the design to be the
  chain that is not among them; node 06 exports the chain node 04 named. The
  `reversed-chain-order` case enforces both.
- **#51** — designs sorted lexicographically, so `design_10` landed between `design_1`
  and `design_2`. Node 04 sorts numerically and nodes 05 and 06 inherit that order.
  `ordering-at-12-designs` enforces it.

Marking a finding open is how a reproduced but unfixed defect gets recorded: add
`"status": "open"` and the suite reports it without failing. Deleting the field enforces
it, and a regression shows up as `FAIL` rather than as silence.

## Adding a fixture

1. Download it and run `probe_entry.py` against it.
2. Add it to `fixtures.json` with the measured `evidence`, the `params` to run it under,
   and a `why` that says what a reader would get wrong without it.
3. Add `expect` entries naming existing checks. A new kind of check goes in the `CHECKS`
   table in `check_fixtures.py` and returns `(status, detail)` — the detail string should
   carry the numbers, so a failure reads as evidence rather than as an assertion.
