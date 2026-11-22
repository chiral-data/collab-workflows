# workflow-038 fixture suite

Seven real RCSB entries run through nodes 01 to 03, the inputs the pipeline must refuse,
and nodes 04, 05 and 06 driven by ensembles built for the purpose — all checked against a
measurement of the input rather than against an exit code. Plus a GPU campaign that runs
the whole thing for real, on a card, against the public MSA server.

## Why this exists

Every one of these structures exits 0 on nodes 01 and 02. Before this suite that was the
whole test: the pipeline ran, so the pipeline worked.

It did not. The defect that started this suite is the clearest case of why an exit code is
not a test. Nodes 02 and 03 declared the alignment as `msa/*.a3m`; Silva matches an inputs
entry **by basename**, so the pattern matched nothing and the directory was never staged.
BioEmu2 builds its own MSA when handed `msa=None`, which means querying the public
ColabFold server from a node that sets `HF_HUB_OFFLINE=1`. Frames appeared, the run went
green, and the manifest recorded the alignment node 02 fetched rather than the one the
model used. Nothing downstream could tell. `node-03-refuses-to-sample-when-the-alignment-never-arrived`
is that case, and it now exits 1 with a message naming the basename rule.

The first runs of this suite found four more things, none of which moved an exit code, and
the GPU campaign then found two the offline suite structurally cannot see (a missing C
compiler and the real VRAM ceiling — both below, under *What the campaign measured*):

- **Any non-standard component needs the network, from inside `system.topology`.** An ion,
  a cofactor, a terminal cap or a modified residue is built from its CCD definition, which
  bioemu2 downloads from `files.rcsb.org` the first time it meets the code. That call sits
  four frames inside the library, in a function nodes 01, 03, 04 and 05 all reach. Five of
  the seven entries hit it: 1CLL needs `CA` and `EOH`, 1A3N `HEM`, 1AKE `AP5`, 5PTI `PO4`,
  `UNX` and `DOD`, 2ZTA `ACE`. Node 01 caught the exception and wrote `n_tokens: null` with a
  one-line warning, so every holo system silently lost its token count — and with it the batch-size
  advice and the >10k GPU-fit guard. Now 126 codes are pre-warmed into the image, each node
  that builds a topology names the component and the cache when a fetch fails, and node 01
  fails instead of warning, because nodes 03, 04 and 05 build the same structure and cannot
  get past what it could not.
- **The collapse warning could never fire.** It tested `max pairwise RMSD < 1e-3 nm`. An
  XTC quantises coordinates to 0.001 nm, so one conformation written twice reads about that
  far apart: measured 0.0018-0.0022 nm on the `identical` fixture across runs. The threshold
  is now 0.01 nm, an order of magnitude above the file format and an order below the smallest
  backbone difference that means anything.
- **Two components nobody would choose to sample.** Running offline made every CCD fetch a
  named failure, and the names were informative. 5PTI — a neutron structure — carries `DOD`,
  deuterated water: `System.from_cif` omits waters written as `HOH` and keeps those written as
  `DOD`, so a neutron entry's solvent is sampled as ligand tokens. It also carries `UNX`,
  "unknown atom or ion", a placeholder for density the depositor could not identify. Both now
  raise a warning naming the code and the `exclude_ligands` string that drops it. `UNL`, the
  unknown-ligand placeholder, cannot be built at all — it raises inside bioemu2's own template
  builder — so node 01 lists the system's non-polymer components and names the parameter instead
  of failing with an `AttributeError`. The advice was a sentence for a while before it was a
  measured fact: two 5PTI fixtures now run `exclude_ligands=UNX,DOD` and `include_ligands=false`
  and check the arithmetic, because advice that does not work is worse than none — a reader who
  tries it once stops reading the rest.
- **A single frame was reported as a collapsed ensemble.** `mode = structure` with
  `num_samples = 1` is the configuration the README recommends for predicting one folded
  structure, and it was never run with a model until this pass: no fixture set `mode`, and the
  campaign runs `equilibrium`. Run on a GPU, 1UBQ gave one frame at Rg 1.1616 nm and the
  bundled Trp-cage 1L2Y one at 0.7146 nm — both compact, both physical, and both announced by
  node 05 as *the ensemble has collapsed*, with advice to check the mode the user had set
  correctly. The cause was a statistic over zero pairs reported as `0.0`, which is exactly the
  collapse signature. `pairwise_rmsd` now carries `n_pairs: 0` with null statistics, node 05
  says a spread needs two frames, and the report's tile reads *n/a, needs 2+ frames*.

  Fixing it almost introduced the opposite defect. One frame from an **equilibrium** run is a
  wasted run — no spread, no clusters, an identically zero RMSF — and the spurious collapse
  warning was the only thing that had ever said so. So the message is now mode-aware, read
  from the mode node 03 recorded rather than guessed from the frame count, and two fixtures
  hold both halves: `single-frame-structure-run-has-no-spread` asserts no warning, and
  `single-frame-equilibrium-run-is-warned-about` asserts there is one. A node that printed
  nothing for either would pass the first fixture on its own.
- **"Polymer residues" was two different numbers.** Node 01 counts the entity sequence,
  which upstream documents as including residues the structure does not resolve; the probe
  counts resolved residues. 1CLL is a 148-residue protein whose crystal resolves 144, and
  the ensemble BioEmu2 produces has all 148. Neither number was wrong; the label was. The
  inventory now records where its number comes from and the report says *entity sequence*.

Two of the findings were in the suite itself, and they belong on this list because a test
that passes by luck is worse than no test.

- **The fixture ensembles were not deterministic.** `make_ensembles.py` seeded its generator
  from `abs(hash((shape, n_frames, n_atoms)))`, and Python salts the hash of a `str` per
  process unless `PYTHONHASHSEED` is set — so every run drew a different ensemble while the
  comment above the line claimed a rerun was byte-identical. It now seeds from a CRC of the
  same key, and the mode is part of it. Verified the way the claim is stated: rebuilding
  `1CLL-6-spread-passed-equilibrium` in a fresh directory reproduced the same
  `samples.xtc`, md5 `88e68a63`, as the suite's own copy.
- **The medoid check asserted something that is not true.** It required `rmsd_to_medoid_nm`
  to be zero at the frame named as the medoid, which is the obvious identity test and wrong:
  mdtraj's RMSD is a float32 QCP, and a frame measured against itself comes back at
  0.0013-0.0018 nm — the same order as the XTC quantisation, and the reason the collapse
  floor is where it is. Combined with the unseeded dice it passed twice and failed on the
  third run, on unchanged node code. The check now asserts what the quantity means: the
  medoid is the argmin, and its own distance sits below the numerical floor rather than at
  zero.

That comparison is the entire design. `probe_entry.py` parses `_atom_site` by hand rather
than through Biopython or biotite, so the measurement does not share a parser — or a bug —
with the code under test.

## Running it

```sh
./tests/run_fixtures.sh                 # offline, no GPU; writes to ./.fixture-run
./tests/run_fixtures.sh /some/dir       # or somewhere else
./tests/run_fixtures.sh /some/dir 1CLL  # one case, by slug
```

Four sections, all driven from `fixtures.json` through `plan.py`, so a case added to the
expectations is run by the driver with no change to the driver:

- **unit** — the three checks that need only the node code and mdtraj: the A3M path
  helpers, the streaming pairwise pass against a dense reference, and frame selection.
- **fixtures** — nodes 01 and 02 over real entries, then node 03 `--dry-run`, with the
  output checked against the measured input. A case may carry parameters, so two cases can
  run the same entry under different ones: 5PTI appears three times, unfiltered, with
  `exclude_ligands=UNX,DOD`, and with `include_ligands=false`. The dry run is three gates, not a parameter
  echo: the alignment is staged, the topology resolves, and each A3M's query row matches
  the sequence the system states for its chain.
- **rejections** — inputs the pipeline must refuse, asserting both the exit status and that
  the refusal gives the stated reason. Three entries check that a legitimate input is still
  accepted, so a guard cannot pass by refusing everything.
- **downstream** — nodes 04, 05 and 06 in sequence over ensembles from `make_ensembles.py`.
  A case declares its frame count, the shape of the spread, the physicality verdict and the
  sampling mode, and `make_ensembles.py` writes the manifest node 03 would have written — so
  a check on a mode-dependent message is driven by the mode, not by a guess from the data.

Every container runs with `--network none`. Node 02 is exercised on its uploaded-alignment
path, with one A3M per protein chain from `make_a3m.py` whose query row is the chain's
measured sequence — so the chain sets and sequences are real and only the server call is
absent. That one is covered for real by the campaign:

```sh
./tests/run_gpu_campaign.sh             # needs a GPU, the MSA server, minutes per system
```

The campaign runs all six nodes over each system, then runs node 04 again with
`propagate=true` on 1UBQ — the only place the time-lagged path executes with a model at all,
since the offline suite reaches node 04's identity pass-through and nothing else. That second
case gets only what the DAG gives node 04, which is node 03's outputs and nothing from nodes
01 or 02, so the checker reads its system facts from the case it was propagated from. Checking
a propagated run against its input is the whole point of running it; reporting four blanks
because the directory is laid out the way the DAG lays it out is not. `NUM_SAMPLES`
and `BATCH_SIZE` are overridable; `BATCH_SIZE=1` is how a system too large for the card gets
measured, and on this card the answer was that batch size is not what the ceiling is made of.

Each node's exit status is written next to its log, which is what lets the checker separate
*this node produced nothing* from *this node was never reached*. Two of the seven systems do
not fit on an 8 GB card, so their rows stop at node 03 and every check below that node reports
as skipped-with-the-reason rather than as a failure.

That concession is narrow on purpose, because a campaign that skips whatever a failed node
would have produced is one step from a campaign that passes because nothing ran.
`stopped only on a resource limit` is the check that holds the line, and it has three
outcomes: a node that stopped without writing an `ERROR` line **fails**, since a silent stop
is the defect class this whole suite exists for; a node that stopped on an out-of-memory or
no-device error **passes**, and only then are the checks below it skipped; a node that stopped
for any other stated reason **fails**, so a real defect cannot hide behind a wall of skips just
because the node described itself politely on the way out.

Needs `docker`, `curl`, `python3` and the image from `../Dockerfile` (override with
`BIOEMU2_IMAGE=...`). Entries are fetched once and cached. Exit status is 0 when nothing is
in the `FAIL` column; known gaps do not fail the run.

As of `bioemu2:2026_11_20` the offline suite reports **70 pass, 0 fail, 0 gap, 2 skipped**
and the campaign **86 pass, 0 fail, 42 skipped**. The two skips offline are both `node02_paired_expectation`
on the heteromers: the offline driver exercises node 02's user-supplied-A3M path, which is
unpaired by construction, so whether a paired alignment is produced when one is due can only
be checked against the live server — which the campaign does, on 1A3N. The campaign's skips
are each printed with their reason, and most belong to the two systems that do not fit the
card.

### What the campaign measured

The campaign samples in `equilibrium` mode. `structure` mode was measured separately, once,
on 1UBQ and on the bundled Trp-cage 1L2Y at `num_samples = 1`: Rg 1.1616 and 0.7146 nm, one
physical frame each, which is how the single-frame collapse warning above was found. It is not
in the campaign loop because one frame has no ensemble to check and the loop's measurements are
all ensemble measurements; the fixture carries it instead.

On an 8 GB RTX 3060 Ti (sm_86, driver 580.178.04), 8 samples per system in one batch, seed 42:
**86 pass, 0 fail, 42 skipped** over the seven systems and the propagated case. Two successive
campaigns on the same card returned every Rg, every max pairwise RMSD and every cluster count
identical to the digit, so the seed reproduces a run rather than merely labelling one — which
is worth knowing before reading any difference between two runs as a result.

The two counts below are different quantities and the table keeps them apart: tokens are what
the GPU pays for (one per polymer residue, one per *ligand atom*), while the Rg scaling is a
relation over **protein residues**. Reading one for the other is what put a wrong expectation
in this table — 5PTI is 65 tokens and 58 residues, 2ZTA 72 and 66. Where two residue counts
appear, the first is the standard amino acids and the second the entity sequence: 2ZTA's
entity is 68 residues, two of which are `ACE` caps rather than amino acids, and the scaling
relation is over the 66 that are. `check_radius_of_gyration` uses the standard count for
exactly that reason.

| System | Tokens | Protein residues | Result |
|---|---|---|---|
| 1UBQ | 76 | 76 | Rg 1.1597 nm (expected 1.14, experimental 1.17), max pairwise C-α RMSD 0.2204 nm, 1 cluster, 8/8 physical, 31.2 s |
| 5PTI | 65 | 58 | Rg 1.1214 nm (expected 1.03, ratio 1.09), max RMSD 0.1698 nm, 1 cluster, 8/8 physical |
| 2ZTA | 72 | 66 of 68 | Rg 1.5244 nm (expected 1.08), max RMSD 0.3342 nm, **4 clusters** — a parallel coiled coil is elongated, so 1.41× the globular scaling is the fold, not a failure |
| 1BNA | 24 | 0 | Rg 1.3226 nm, max RMSD 0.0623 nm, 492-atom all-atom selection — the no-protein fallback, exercised end to end. 24 tokens for a B-DNA dodecamer duplex: 12 nucleotides per strand, and no Rg expectation, because the globular relation is for proteins and the check skips it |
| 1CLL | 155 | 148 | Rg 1.5677 nm (expected 1.47, ratio 1.07), max RMSD **1.3682 nm** over **8 clusters**, 148-C-α selection with the 4 calcium ions excluded. Calmodulin's two lobes pivot about a flexible central helix, so a broad multi-state ensemble is the literature answer — and the contrast with ubiquitin's single 0.22 nm cluster is the closest thing here to a scientific control. It is also the system that needed a C compiler in the image, which 1UBQ, 5PTI and 2ZTA did not |
| 1AKE | 542 | 214 per chain, 428 modelled | **out of memory, at batch 8 and again at batch 1**, where upstream's estimate said 33 would fit. Both attempts died on the same 1.06 GiB allocation |
| 1A3N | 746 | 574 | **out of memory, at batch 8 and again at batch 1**, estimate 19; the same 1.46 GiB allocation both times. Node 02's paired-alignment path works against the live server even so: `paired: true`, unpaired depths 4522/5023 per chain and a paired alignment of 8794 rows that is the same depth on all four chains, as a positionally paired alignment must be, and node 01 measured 746 tokens |
| 1UBQ, propagated | 76 | 76 | The time-lagged path, run for real: 8 frames at lag 1.0 ns, 601 atoms unchanged, 8/8 physical, Rg 1.1628 nm against 1.1597 before propagation and spread 0.1943 nm against 0.2204. Propagating a rigid fold by 1 ns should move it slightly and should not unfold it, and that is what the numbers say. It is a weak test of the dynamics and a real test of the plumbing — nothing else in either suite runs `propagate_with_embeddings` |

Three of those are findings rather than results. The C compiler is a runtime dependency of
triton's kernel launcher that `python:3.12-slim` does not carry, and the first system to need
it failed after the model had loaded, with a message about `CC` that says nothing about size.
2ZTA is the reason the Rg tolerance stays at 1.5×: tightening it would flag a real coiled coil.

And the two OOMs are not a batch-size problem, which is the thing the first reading of them
got wrong. Re-run at `BATCH_SIZE=1`, both systems failed on the *same* allocation as at batch
8 — 1.06 GiB for 542 tokens, 1.46 GiB for 746 — so on this card the peak is set by the token
count, not by how many samples are denoised together. Upstream's formula is not merely
optimistic here, it is measuring the wrong variable: it scales the batch and the ceiling does
not move. The node's advice now says to try a lower batch but not to expect much, and names
the token count as what to come down on instead. 1CLL at 155 tokens fits with room to spare,
so somewhere between 155 and 542 tokens is where an 8 GB card stops being enough; the suite
brackets it rather than claiming the number.

The pieces run on their own too:

```sh
python3 -I tests/probe_entry.py .fixture-run/downloads        # measure the inputs
python3 -I tests/check_fixtures.py .fixture-run/results .fixture-run/downloads
python3 -I tests/check_campaign.py .gpu-campaign/results .gpu-campaign/downloads
```

Downloaded structures are untrusted input. They stay in their own directory, nothing is
executed from it, and every `python3` here runs with `-I` so a file dropped next to a
script cannot be imported in place of a stdlib module.

## What this does not establish

The offline suite tests the plumbing, not the science. Passing means the pipeline does not
silently corrupt, drop, mislabel or miscount its input, and that its guards fire where they
should and stay quiet where they should.

Node 03's sampling is not run offline: it needs a GPU and the 1.9 GiB checkpoint. The dry
run still exercises the staging, the parameter resolution, the alignment guard and — since
the topology is resolved before the model is loaded — the chemistry, which is where every
plumbing defect lived. Nodes 04, 05 and 06 are run in full, on synthetic ensembles, which
is what makes them deterministic: `make_ensembles.py` builds the XTC from the real
`System.topology`, so atom order matches by construction, but the coordinates are a spiral
plus noise. They are not protein-like, and they do not need to be — every quantity node 05
computes is a function of coordinates alone. Nothing offline exercises the code that *fills*
the sampling manifest, because the fixture supplies one.

The rotations in `make_ensembles.py` are the one place a synthetic ensemble tests something
a real one could not. A frame that differs from another only by a rigid rotation is not a
different conformation, so RMSD, RMSF and the collapse warning must all be invariant to it
— and `ensemble_top.pdb` must hold the *superposed* frames, or the report's 3D viewer
animates the molecule tumbling and a reader sees conformational change that is not there.
With the rotation present that bug is visible. Without it, it is not.

What the offline suite cannot see at all: whether the model works. Frame counts, cluster
counts, physicality records and chain labels are all computable on an unfolded or exploded
chain and all look ordinary. That is what `check_campaign.py` is for — it compares the
measured radius of gyration against 0.22·N^0.38 nm, the globular scaling, which is the one
number that stops being ordinary when the checkpoint, the precision or the alignment is
wrong. 1UBQ measured 1.1597 nm against 1.16 nm expected and an experimental 1.17 nm.

Nothing here reproduces a published ensemble, compares a free-energy surface, or validates
a state population. The entries were chosen for the input shapes they carry. 1AKE was picked
to be the scientific positive control — adenylate kinase is a textbook two-state enzyme and
upstream's own example — but it does not fit the reference card, so the control that actually
ran is 1CLL: calmodulin's two lobes pivot about a flexible central helix, and it comes out
multi-state (8 clusters, 1.37 nm) where ubiquitin comes out single-state (1 cluster, 0.22 nm).
Even there the claim is only that the ensembles are physically sensible and differ in the
direction the literature says they should, not that the populations are right. Frame counts
are not populations, and the report says so on every page.

Not every input path is covered here. Node 01 accepts an upload, a structure file, a raw
`sequence` or an RCSB `pdb_id`, and until this pass the suite staged a structure file every
time, so three of the four advertised ways in had never been taken. `sequence` is covered now
— `ubiquitin-from-a-bare-sequence` stages no file at all and asserts that the inventory
matches the length of the string, which on that path is the entire truth about the system.
`pdb_id` needs the network, so it was measured once by hand instead: 1UBQ resolves to
`RCSB:1UBQ` at 76 residues and 76 tokens, with the log recording the assembly and the ligand
policy. `rmsd_cluster_threshold_nm` is covered too, by
`cluster-count-follows-the-threshold-asked-for`, which re-derives the cluster count from the
published matrix under node 05's own greedy rule — tie-break included, since a check that
broke ties the other way would be re-deriving a different algorithm. Its threshold comes
from the fixture's own measured distance distribution, not from a round number: the 66
pairwise distances run 0.1146 to 0.1338 nm, so the count is 12 below 0.12 nm, 4 at 0.12 and
1 above 0.13. 0.12 is the only value where both ignoring the parameter and applying it at
the wrong scale change the answer. `msa_depth` and
`use_env` remain unexercised; both apply only to an alignment fetched from the server, which
the offline driver never does. The parameters that *are* covered are
in the fixture table, and the two ligand filters were added because node 01 prints their
names as advice — advice that nothing had ever run.

## Choosing a fixture

A fixture earns its place by **measurement**, never by reputation. Each entry in
`fixtures.json` carries an `evidence` string recording what `probe_entry.py` reported from
the downloaded file. Re-run the probe before trusting it; RCSB re-releases entries.

Measured on 2026-11-20:

```
entry        chains  prot  seqs    res    CA  Ca2+   SG
1A3N.cif         12     4     2    572   572     0    6
1AKE.cif          6     2     1    428   428     0    2
1BNA.cif          4     0     0     24     0     0    0
1CLL.cif          7     1     1    144   144     4    0
1UBQ.cif          2     1     1     76    76     0    0
2ZTA.cif          4     2     1     64    62     0    0
5PTI.cif          4     1     1     58    58     0    6
```

The `Ca2+` column is the reason the probe exists. A calcium ion is an atom named `CA` in a
residue named `CA`, and mdtraj's `name CA` matches **atom names**, so the obvious selection
returns 148 atoms where 1CLL has 144 C-alphas — shifting every RMSD, the RMSF, the cluster
count and the medoid frame, with nothing in the output to show it. Four atoms. No exit code
can see that; `len(residue_rmsf) == 144` can.

## What each fixture is for

| Entry | Shape | Catches |
|---|---|---|
| 1UBQ | ubiquitin, one chain, no ligand | control — a failure here is the harness; and the only entry that needs no CCD fetch |
| 1CLL | calmodulin, 4 structural Ca2+ and an ethanol | calcium counted as C-alpha; 148-vs-144 entity-vs-resolved residues; a cryoprotectant sampled as a ligand |
| 5PTI | BPTI by neutron diffraction: 3 disulfides, a phosphate, a `UNX` placeholder and a `DOD` | the sulfur physicality thresholds are live; three separate CCD fetches; a component of unknown chemistry, and heavy water that upstream does not recognise as water |
| 2ZTA | GCN4 leucine zipper, homodimer, N-terminal `ACE` | two chains and one sequence, so **no** paired alignment is due; per-chain RMSF traces; a modified-residue CCD |
| 1A3N | haemoglobin, 4 chains over 2 sequences, 4 haems | a paired alignment **is** due; the heteromer path; cofactors at the interface |
| 1AKE | adenylate kinase, 2 chains, AP5A inhibitor | upstream's own example and a real two-state enzyme, chosen as the case where the ensemble should be interesting rather than merely well formed. At 542 tokens it does not fit the 8 GB reference card at any batch size, so on this hardware it is an input-measurement and failure-reporting case, and 1CLL carries the interesting ensemble instead |
| 1BNA | B-DNA dodecamer, no protein chain | no alignment is due and none must be demanded — the counter-case to the staging guard |
| 5PTI with `exclude_ligands=UNX,DOD` | the same entry, placeholder and heavy water dropped | node 01's own warnings name this parameter as the remedy, and nothing had ever run it. Measured: 4 chains and 65 tokens become 2 chains and 63, with `PO4` surviving and the two excluded codes gone — one token each, which is the arithmetic the check asserts. One token each despite the file holding 63 D₂O molecules over 189 atom rows: the loader presents them as a single `DOD` component, so what gets sampled is not the solvent shell |
| 5PTI with `include_ligands=false` | the same entry as an apo protein | the stronger assertion: with every component dropped the token count has to fall to the polymer residue count, and it does — 58 and 58 |
| ubiquitin from a bare `sequence` | the raw-sequence input path, no file staged | three of node 01's four advertised inputs had never been run. With no file to parse and no entity tables to read, the length of the string is the whole truth: one chain, 76 residues, 76 tokens, no components, and the chain's stated sequence equal to the input |

## Known gaps

Expectations marked `"status": "open"` in `fixtures.json` report as `GAP` and do not fail
the run. They are real, reproduced, and not yet fixed.

**Nothing is marked open today.** The last one closed in this pass: an uploaded A3M whose
query row was not the system's sequence used to be used anyway, since neither `get_msa` nor
`load_a3m_msa` compares the two, so the alignment conditioned the model on a different
protein and the run looked entirely ordinary. Node 03 now compares every staged A3M's query
against the sequence the System states for that chain, before the model is loaded, and
`node 03 refuses an alignment whose query is not the system sequence` enforces it. A chain
whose stated sequence carries `X` — a modified residue or a cap, which an alignment cannot
represent — is reported rather than refused; 2ZTA is the fixture that exercises that path.

Everything tracked here is closed and enforced:

- the `msa/*.a3m` staging failure — `node-03-refuses-to-sample-when-the-alignment-never-arrived`,
  with `node-03-accepts-an-absent-alignment-when-none-was-due` as the counter-case.
- the paired alignment dropped between nodes 02 and 03 — `check_paired_alignment` in the
  campaign, keyed off the measured number of unique protein sequences.
- `medoid = 0` above 200 frames under the name `rmsd_to_medoid_nm` —
  `medoid-survives-above-the-matrix-cap`, which recomputes the argmin.
- the un-superposed viewer PDB — `check_viewer_frames_superposed` in the campaign, and the
  rotations in `make_ensembles.py` offline.
- `n_tokens` read from the wrong attribute — `inventory_tokens_measured`.
- `edm_fk` and `edm_fkc` handed to propagation, which cannot resample across differing
  frames — two refusals and `edm_gd` as the counter-case.
- the three propagation states that must not collapse into "propagation was disabled" —
  `node05_refuses_missing_manifest`, `node05_refuses_lost_ensemble`,
  `node05_passthrough_not_treated_as_lagged`.
- node 05 looping forever when `rmsd_cluster_threshold_nm` is finer than the RMSD
  measurement's own noise — `clustering-below-the-rmsd-noise-floor-terminates`.
- the campaign reporting a green run for a system whose first node died —
  `an-abort-never-skips-the-check-that-judges-it`.
- the campaign judging whichever directories happened to be in the results tree, rather than
  the cases the driver ran — `run_gpu_campaign.sh` now records them in
  `results/campaign_cases.txt` and the checker reads it.

Those last two came from a review of this suite rather than of the workflow, and both were
defects in the testing, which is the part nothing else tests:

- **`greedy_clusters` could not terminate below the RMSD noise floor.** `mdtraj.rmsd` is a
  float32 QCP, so a frame measured against itself reads 0.0013–0.0018 nm rather than 0. The
  membership test was `matrix[centre, i] <= threshold`, so at `rmsd_cluster_threshold_nm =
  0.001` — a value the parameter accepts, validated only as positive — the centre failed its
  own test, `absorbed` came back empty, the unassigned set never shrank, and node 05 spun on
  the GPU until something killed it. The centre is now absorbed unconditionally, because a
  frame is in its own cluster by definition, and the node warns when the threshold is below
  the floor instead of silently returning one cluster per frame.
- **the campaign checked directories, not cases.** `check_campaign.py` took every directory
  under the results tree as a case. That cut both ways: a scratch directory left behind
  reported as a system failing every check after node 01 — eight of them, from this session's
  own parameter probes, turned a 86/0/42 report into 104/54/98 — and a system the driver was
  asked for but never reached left no directory at all, so the report said nothing about it
  rather than saying it was missing. The driver now writes `campaign_cases.txt` and the
  checker reads that: an undeclared directory is named and ignored, a declared case with no
  directory is a failure, and the run exits non-zero. Same rule as the rest of the pipeline —
  the component that knows a fact records it, and the one downstream reads it instead of
  inferring it.
- **the campaign's abort handling skipped its own judge.** `check_campaign.py` converts FAIL
  to SKIP for any check needing a node at or after the one that aborted, since an output
  nobody wrote is not evidence. `check_failure_was_reported` is the check that decides
  whether a stop was a reported resource limit or a silent death — and it reads node 01, so
  an abort *at* node 01 satisfied the conversion test and skipped the judge along with all
  fifteen checks it was meant to judge: 0 pass, 0 fail, 16 skipped, exit 0. The judge is now
  excluded by identity rather than by node index. A comment at that line had asserted the
  invariant held; it was never tested, and it was false.

Marking a finding open is how a reproduced but unfixed defect gets recorded: add
`"status": "open"` and the suite reports it without failing. Deleting the field enforces
it, and a regression shows up as `FAIL` rather than as silence.

One gap is open and is not a defect, so it carries no fixture: **no heteromer finishes the
pipeline on the reference card.** 1A3N is the entry that would do it, and the 8 GB card
cannot hold 746 tokens at any batch size. What is measured on it is node 01's inventory,
node 02's paired alignment against the live server, and node 03's own log line attaching
that paired alignment to 4 chains — the attach happens, then the allocation fails. So the
paired-alignment path is verified up to the point the model is called and no further: no
multi-sequence ensemble, no per-chain RMSF over two different sequences, no clustering of a
heteromer. Closing it needs a larger card, not a change here, and `check_paired_alignment`
is already written to check it the moment one runs.

## Adding a fixture

1. Download it and run `probe_entry.py` against it.
2. Add it to `fixtures.json` with the measured `evidence`, the `params` to run it under,
   and a `why` that says what a reader would get wrong without it.
3. Add `expect` entries naming existing checks. A new kind of check goes in the `CHECKS`
   table in `check_fixtures.py` (or `check_campaign.py`) and returns `(status, detail)` —
   the detail string should carry the numbers, so a failure reads as evidence rather than
   as an assertion.
4. If it carries a component outside the 126 codes pre-warmed in the `Dockerfile`, either
   add the code there or expect the offline run to fail naming it.
