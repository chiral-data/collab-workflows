# Measured reference files

Real ColabFold output, captured on the GPU host so that node 03's flags and node 04's
parser are written against a **measured** file rather than a guess. Nothing in this
directory is synthetic.

037's `54dc7def` is why this exists: it read the real package before building and found
that several of the documented guesses the workflow had been built on were wrong. Until
these files existed, everything this repo knew about the scores-file shape came from
*reader* code in workflows 026 and 029.

## Provenance

| | |
|---|---|
| Host | `amd-c12-1`, NVIDIA GeForce RTX 3060 Ti 8 GB, driver 580.178.04 |
| Image | `colabfold:2026_10_08` (built from `workflow-026/Dockerfile`) |
| ColabFold | 1.6.1 |
| jax | 0.5.3, `CudaDevice(id=0)` |
| Captured | 2026-11-21 |
| Network | `docker run --network none` — **every file here was produced with no network at all** |

## The control run write-up

### `tier2-phase1-controls.md`

Not a captured artifact but the measured results of the Tier 2 Phase 1 control run: five
controls through all five nodes on the GPU, offline, twice. It is here rather than in
`tests/README.md` because it is evidence with provenance attached, and because two of its
three findings revise decisions in `validation_pipeline.md` rather than describing the
suite.

The short version: the scorer is right on a designed binder (0.93 A at 95.05 pLDDT,
exported), the native positive controls were the wrong controls for a workflow that
predicts at `single_sequence`, and `plddt_binder` cannot gate alone because 115 alanines
score 95.81.

## Files

### `monomer_scores_rank_001.json`

A 35-residue sequence, `--model-type alphafold2_ptm --msa-mode single_sequence
--num-models 1 --num-recycle 1`. 61 s wall clock, `rc=0`.

Keys, measured:

```
plddt     list[35]      flat, per-residue
pae       list[35][35]  full square matrix
max_pae   float         21.15625
ptm       float         0.57
```

**`iptm` is absent.** So is any interface quantity. This is the normal Phase 1 case, not
an edge case, and it is why `04-score` must follow 026's
`iptm = float(data["iptm"]) if "iptm" in data else None` rather than 029's
`scores.get('iptm', 0.0)` — the latter reports a measured-looking 0.0 for every design in
every monomer run.

Note this file is the **raw** ColabFold output. `workflow-026/sample_outputs/confidence_summary.json`
is 026's own post-processed summary, with different key names
(`plddt_per_residue`, `plddt_mean`, `interface_pae_mean`). Parsing the raw file against
the summary's names would fail.

### `multimer_scores_rank_001.json`

A 34+29 complex written `TARGET:BINDER` in one record, `--model-type
alphafold2_multimer_v3 --msa-mode single_sequence`. 80 s, `rc=0`.

```
plddt     list[63]      flat, spanning BOTH chains, concatenated
pae       list[63][63]  spans both chains; the interface blocks are the off-diagonals
max_pae   float         25.375
ptm       float         0.47
iptm      float         0.2      <- present only for multimer
```

The predicted PDB carried chain A with 34 residues and chain B with 29. So
**`TARGET:BINDER` puts the target in chain A and the binder in chain B**, and the binder is
the **tail** of the `plddt` array. The `plddt-slice-target-first` fixture exists for
exactly this: 029 slices `plddt[:binder_length]` because its binder is chain A, and a node
04 inheriting that would report the target's confidence as the binder's.

Neither file carries a chain boundary anywhere. The split has to come from the recorded
chain lengths, which is why node 02 records them and node 04 asserts
`len(plddt) == sum(chain_lengths)`.

### `complex_target_only_msa.a3m`

**The §11.1 answer, and the file that unblocks Phase 2.** A paired a3m in which the target
carries an alignment and the binder is present only as its own query sequence. Fed to
`colabfold_batch` it produced a genuine two-chain prediction (34 A / 29 B, `iptm` 0.14,
`coverage.png` written, so an MSA was really used), offline.

The format was **derived from `colabfold.batch.unserialize_msa()`'s source**, after two
hand-rolled guesses were rejected with `ValueError: MSA 0 must contain at least one
sequence`:

```
#<len1>,<len2>\t<card1>,<card2>
>101
<chain1 query><chain2 query>        <- CONCATENATED, no gaps. This is a3m_lines[2],
                                       sliced by the declared lengths. Residues in both
                                       segments, so it becomes the PAIRED row.
>101
<chain1 query><gaps>                <- single-chain row: seeds unpaired_chunks[0]
>102
<gaps><chain2 query>                <- single-chain row: seeds unpaired_chunks[1]
>homolog_1
<chain1 homolog><gaps>              <- chain 1's alignment
...
```

The rule that both guesses missed: `unserialize_msa` sends any row with residues in two or
more chain segments to the **paired** block for every chain, and only single-chain rows
populate `unpaired_chunks[j]`. **Every chain must therefore carry at least one single-chain
row**, including its own query, or `generate_input_feature` raises. A concatenated query
alone leaves every unpaired block empty.

The homologs in this captured file are synthetic variants, used to probe the mechanism
without sending anything to a third party. A real run gets the target's alignment properly
— that is `01b-align-target`'s job — but the *format* question this file answers was the
gate, and it is now closed.

### `monomer_unrelaxed_rank_001.pdb`

The rank-1 structure from the monomer run, kept so `04-score`'s PDB reader and
`make_predictions.py`'s synthetic output can be checked against a real one.

## Re-capturing

Re-run and re-measure rather than trusting this page; RCSB re-releases entries and images
get rebuilt. The commands are in `tests/README.md`. If a shape here changes, the files and
this page change in the same commit as the parser.
