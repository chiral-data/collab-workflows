# workflow-039 — Antibody Developability Triage (design doc)

Status: **implemented, never executed.** All seven nodes are written; the image has not been
built and no node has been run, because none of the dependencies are installed outside the
container. Every script compiles, and the pure-logic parts have unit tests (see
[Verification](#verification-plan)), but no number in this pipeline has been produced by running
it end to end. Treat the verification plan as outstanding work, not a record of results.

## Context

Nothing on `main` handles antibodies as antibodies. The repo can dock them
(workflow-016 DiffDock-PP, -017 LightDock), predict structures of them as generic proteins
(-012 Boltz-2, -026 ColabFold, -018 Boltz-2 vs Chai-1) and design binders de novo (-013, -033,
-037), but none of those knows what a CDR is, and none answers the question a biotech with a
panel in hand actually asks: **which of these clones do I express?**

There is one piece of prior art, off `main`: the unmerged `workflow-abodybuilder3` branch
(June 2026) rebuilds workflow-010 as an ABodyBuilder3 Fv structure predictor, with PLM embedding
and its own report. It overlaps squarely with node 02 here. It was not known to the author when
this workflow was designed; the overlap is recorded as open item 8 rather than quietly ignored.

That question is answered today by developability triage: number the sequences, model the Fv,
compute the handful of structure-derived risk metrics that correlate with the clinical-stage
distribution, scan for sequence liabilities, check humanness, and rank. It is the standard
first-pass filter in antibody discovery, it is cheap, and it is almost entirely CPU work.

## Corrections from the design recheck (2026-11-11)

The first draft of this document was written before the TAP paper's full text and the upstream
source were read closely. Five things were wrong, and all five are fixed in the implementation:

1. **TAP uses IMGT CDR definitions, not North.** The paper numbers with IMGT and splits by IMGT
   CDR definitions; North is used only for its canonical-form clustering. `cdr_definition`
   therefore defaults to `imgt`, and the CDR vicinity is built from IMGT CDRs plus anchors.
2. **TAP's absolute thresholds are published, in Table 2.** The first draft assumed only
   percentile language existed (Table 1) and left the thresholds untranscribed, defaulting
   `tap_threshold_set` to `none`. Table 2 gives numeric amber/red bounds over the full 242-CST
   set, so the thresholds are implemented with a citation and the default is now `clinical`.
3. **The metric formulas are fully specified in Methods** — PSH, PPC, PNC and SFvCSP, plus the
   surface-exposure, CDR-vicinity and salt-bridge definitions. Node 03 implements them as
   written rather than approximating.
4. **ImmuneBuilder's `save()` always refines.** There is no "skip refinement" flag, so
   `refine_structures = false` has to call `save_single_unrefined()` — which, unlike `save()`,
   does **not** write the error estimate into the B-factor column. Node 02 adds it explicitly on
   that path so the two paths produce comparable files.
5. **The exposure reference is a genuine deviation, not a detail.** TAP measures *side-chain*
   relative exposure against a Yang Zhu & Blundell reference we do not have. We use whole-residue
   relative ASA against Tien et al.'s theoretical Ala-X-Ala values. This shifts exposure calls,
   which shifts every metric built on them. See [the calibration
   problem](#the-calibration-problem).

## Scope decisions

| Decision | Choice | Consequence |
|---|---|---|
| Headline capability | **Panel triage, not single-antibody analysis** | Input is a panel of clones; every node is a loop; output is a ranked shortlist. A one-antibody run is just n=1. |
| Antigen | **Out of scope** | No epitope mapping, no paratope-from-complex analysis, no affinity prediction. Those need a complex predictor — point users at workflow-012/-016/-018. See [why the antigen is out of scope](#why-the-antigen-is-out-of-scope). |
| Compute | **CPU only** | No GPU node at all. The differentiator: a few hundred clones for the cost of a CPU job. |
| Structure model | **ImmuneBuilder (ABodyBuilder2 / NanoBodyBuilder2)** | Seconds per antibody vs minutes for a general folder. Nanobodies covered; TCRs deliberately not exposed although TCRBuilder2 ships in the same package. |
| Numbering | **Classic ANARCI, via conda** | Forced by ImmuneBuilder — see [the ANARCI / ANARCII split](#the-anarci--anarcii-split). |
| CDR definition | **IMGT by default**, matching TAP | Other definitions selectable, each with the numbering scheme it is defined on. |
| Developability metrics | **TAP's published formulas, reimplemented; thresholds from Table 2** | Still labelled "TAP-style", never "TAP" — see [TAP is a paper, not a package](#tap-is-a-paper-not-a-package). |
| Humanness | **promb (MIT), not BioPhi's OASis** | OASis needs a 22 GB database download; promb's `human-oas` DB is the >=10% prevalence peptide set, pip-installable. |
| Immunogenicity | **MHC-II prediction omitted entirely** | No open, commercially usable class-II predictor exists. See [the MHC-II hole](#the-mhc-ii-hole). Humanness is a proxy and is presented as one. |
| Liabilities and flags | **Reported, never enforced** | No clone is dropped, ever. Flags are columns; ranking is a sort. Mirrors workflow-038's treatment of physicality. |
| Mutation proposals | **AbLang2 probabilities, re-screened** | Proposals are hypotheses for the bench, explicitly not affinity predictions. |

## Upstream facts this design is built on

Every item was read from PyPI, the GitHub API, upstream source or the paper's full text on
**2026-11-11**, not recalled.

### Numbering

**ANARCI (classic)** — `oxpig/ANARCI`, BSD-3-Clause, 299 stars, last push **2024-05-21**. PyPI
`ANARCI` **2026.2.13.2** (2026-02-13), `requires-python >=3.8`, declares only `biopython`; the
HMMER binary and germline HMMs are not pip dependencies, which is why upstream recommends
conda/bioconda. Signature:

```python
anarci(sequences, scheme="imgt", database="ALL", output=False, outfile=None, csv=False,
       allow=set(["H","K","L","A","B","G","D"]), hmmerpath="", ncpu=None,
       assign_germline=False, allowed_species=['human','mouse'], bit_score_threshold=80)
```

Returns `(numbered, alignment_details, hit_tables)`. `numbered[i]` is a list of domains or
`None`; each domain is `(numbering, start_index, end_index)` where `numbering` is a list of
`((position, insertion), aa)` with `'-'` for gaps. `alignment_details[i][j]` carries `species`,
`chain_type`, `evalue`, `bitscore`, `scheme`, and — only when `assign_germline=True` —
`germlines`, as `{'v_gene': [[species, gene], identity], 'j_gene': [...]}`. Allowed species:
human, mouse, rat, rabbit, rhesus, pig, alpaca. Schemes: imgt, chothia, kabat, martin (aho via
the wrapper).

**ANARCII** — `oxpig/ANARCII`, BSD-3-Clause, last push **2026-11-09** (actively maintained),
PyPI `anarcii` **2.0.8** (2026-06-30), `requires-python >=3.11`. Language-model numbering, the
successor to ANARCI, no HMMER. `abnumber` 0.4.4 (MIT) has already migrated to it.

**CDR boundary tables** — from `prihoda/abnumber` (MIT), `abnumber/common.py`, as start
coordinates `[CDR1, FR2, CDR2, FR3, CDR3, FR4, end]`:

```text
imgt:      [27, 39, 56, 66, 105, 118, 129]
kabat_H:   [31, 36, 50, 66,  95, 103, 114]    kabat_K/L:   [24, 35, 50, 57, 89, 98, 108]
chothia_H: [26, 33, 52, 57,  95, 103, 114]    chothia_K/L: [24, 35, 50, 57, 89, 98, 108]
north_H:   [23, 36, 50, 59,  93, 103, 114]    north_K/L:   [24, 35, 49, 57, 89, 98, 108]
```

Each definition is expressed in its own numbering scheme — North is defined on Chothia
numbering, hence the `CDR_DEFINITION_SCHEME` map in node 01. Anchors are derived as the
positions immediately flanking each CDR, not hardcoded.

### Structure prediction

**ImmuneBuilder** — `brennanaba/ImmuneBuilder`, BSD-3-Clause, last push **2024-11-08**, PyPI
**1.2** (2024-06-25). `install_requires`: numpy, scipy>=1.6, einops>=0.3, torch>=1.8, requests.
Needs **OpenMM + pdbfixer** for refinement and ANARCI for trimming/numbering, both conda-first
per upstream. Weights **download at runtime** (`download_file`, `are_weights_ready` in
`util.py`), so they are baked at build time.

- `ABodyBuilder2(model_ids=[1,2,3,4], weights_dir=None, numbering_scheme='imgt')`,
  `.predict({'H': vh, 'L': vl})` → `Antibody`. `NanoBodyBuilder2` likewise with `{'H': vhh}`
  (it sets `L` to `""` internally).
- `sequence_checks.number_sequences` does `from anarci import validate_sequence, anarci,
  scheme_short_to_long`, asserts `len(sequence) > 70`, and defaults to
  `allowed_species=['human','mouse']` (ABodyBuilder2) or `None` (NanoBodyBuilder2).
- `Antibody.numbered_sequences[chain]` is a list of `((position, insertion), aa)` — the same
  shape ANARCI returns, in the predictor's own scheme, and exactly what `to_pdb` writes out. So
  PDB residue ids match node 01's numbering when both use the same scheme.
- `Antibody.error_estimates` is the squared deviation of the four-model ensemble about its mean;
  `error_estimates.mean(0).sqrt()` is the per-residue error, and `ranking` orders models by it.
- **`save()` always calls `refine()`** and then `add_errors_as_bfactors()`.
  `save_single_unrefined(filename, index)` skips refinement **and** the B-factor step.

### Language model

**AbLang2** — `oxpig/AbLang2`, PyPI `ablang2` **0.2.1** (2024-06-04), BSD-3-Clause. Deps:
torch>1.9, requests, einops, rotary-embedding-torch, numpy. Downloads weights on `pretrained()`.

```python
ablang = ablang2.pretrained(model_to_use='ablang2-paired', random_init=False, ncpu=1, device='cpu')
probs = ablang([f"{vh}|{vl}"], mode='probability')   # list of [len(seq), vocab] arrays
```

Valid modes: `rescoding`, `seqcoding`, `restore`, `likelihood`, `probability`,
`pseudo_log_likelihood`, `confidence`. `probability` is `logits.softmax(-1)`;
`stepwise_masking=True` masks each position in turn (one forward pass per position instead of
one per sequence). Output rows are trimmed by `res_to_list(logits, seq) = logits[:len(seq)]`, so
**row i is character i of the `"VH|VL"` input string**, with the separator occupying row
`len(vh)`. Vocabulary from `ablang.tokenizer.aa_to_token` (`ablang2/models/ablang2/vocab.py`):
`<`=0, `-`=21, `>`=22, `*`=23 (mask), `X`=24, `|`=25 (separator), and the 20 amino acids 1–20.
`align=True` is **paired-only** and asserts otherwise, so this workflow does not use it.

### Humanness

**promb** — `MSDLLCPapers/promb`, PyPI **1.0.2** (2025-05-17), MIT, `requires-python >=3.8`,
deps scipy/pandas/tqdm only. `init_db('human-oas')` is the antibody peptide set found in
**>=10% of human subjects**; also `human-reference` (UP000005640) and `human-swissprot` at a
chosen k. `db.compute_peptide_content(seq)` is the **OASis-like score**;
`db.chop_seq_peptides(seq)`, `db.contains(peptide)` and
`db.find_nearest_peptides(peptides, n)` (ranked by identity, ties by BLOSUM) are the rest of the
surface. Documented example: the VH in `input_files/README.md` scores
**0.7142857142857143**.

**BioPhi** — `Merck/BioPhi`, MIT, last push **2025-05-13**. **Not on PyPI**
(`pypi.org/simple/biophi/` → 404); bioconda only. Local OASis humanness requires the **22 GB
uncompressed** `OASis_9mers_v1.db` from Zenodo record 5164685. Upstream's own README points
users who only want the >=10% prevalence score at promb instead.

**mhcflurry** — PyPI **2.3.25** (2026-11-06), Apache-2.0. **MHC class I only.**

### TAP (Raybould et al., PNAS 2019, DOI `10.1073/pnas.1810576116`, CC BY)

Read from the Europe PMC full text (PMC6410772). Definitions, verbatim in substance:

- **Surface-exposed**: "≥ 7.5% relative exposure across side-chain atoms, compared with the
  open-chain form alanine-R-alanine, as calculated with the Shrake and Rupley algorithm".
- **CDR vicinity**: "every surface-exposed IMGT-defined CDR and anchor residue, and all other
  surface-exposed residues with a heavy atom within a 4-Å radius".
- **Salt bridges**: "pairs of lysines/arginines and aspartic acids/glutamic acids with a
  N⁺−O⁻ distance ≤ 3.2 Å".
- **PSH**: over pairs of surface-exposed residues with closest heavy-atom distance
  `r₁₂ < 7.5 Å`, `Σ H(R₁,S)·H(R₂,S) / r₁₂²`, where `H` is the hydrophobicity score normalized
  to **[1, 2]**. Kyte & Doolittle is the scale used for all reported values. **Salt-bridge
  residues take glycine's value.**
- **Charges**: D −1, E −1, K +1, R +1, H +0.1 (Henderson–Hasselbalch, pKa 6, pH 7.4). Tyrosine
  deprotonation not considered. **Salt-bridge residues take charge 0**, and buried side chains
  are ignored entirely.
- **PPC / PNC**: "analogous in form to PSH, with H(R,S) substituted for |Q(R)|" — PPC over
  positively charged residues, PNC over negatively charged ones.
- **SFvCSP**: `(Σ Q over surface-exposed V_H residues) × (Σ Q over surface-exposed V_L residues)`.
- **Total CDR length**: IMGT-numbered, IMGT CDR definitions, all six CDRs summed.
- **Thresholds, Table 2** ("as defined by the entire set of 242 CSTs"):

  | Metric | Amber | Red |
  |---|---|---|
  | Total CDR length | 54 ≤ L ≤ 60 | L > 60 |
  | PSH, CDR vicinity | 83.84 ≤ PSH ≤ 100.71 and 156.20 ≤ PSH ≤ 173.85 | PSH < 83.84 or PSH > 173.85 |
  | PPC, CDR vicinity | 1.25 ≤ PPC ≤ 3.16 | PPC > 3.16 |
  | PNC, CDR vicinity | 1.84 ≤ PNC ≤ 3.50 | PNC > 3.50 |
  | SFvCSP | −20.40 ≤ SFvCSP ≤ −6.30 | SFvCSP < −20.40 |

  Table 1 states the intent: amber means "within the extremes of the distribution", red means "a
  previously unobserved value".
- Reference CST distribution means, for sanity-checking our implementation: CDR vicinity PSH
  **123.30 ± 16.60** across 137 CST models; 88.32% of CSTs have PPC < 1 and 80.30% have PNC < 1;
  21.9% have negative SFvCSP.
- TAP's models were built with **ABodyBuilder v1**, and the paper notes modelling biases PSH
  **upward** relative to crystal structures.
- Distributed as a web application; a search on 2026-11-11 found **no public source repository**.

### Supporting data

**Tien et al. 2013**, *Maximum allowed solvent accessibilites of residues in proteins*
(PMC3836772), Table 1 — Ala-X-Ala theoretical **whole-residue** max ASA in Å², used as our
exposure denominator. These are whole-residue, not side-chain, values; that is the deviation in
[the calibration problem](#the-calibration-problem).

**Kyte & Doolittle scale** — taken from `Bio.SeqUtils.ProtParamData.kd` rather than transcribed,
so the numbers come from biopython.

**biotite** — `biotite.structure.sasa(array, probe_radius=1.4, atom_filter=None,
ignore_ions=True, point_number=1000, point_distr="Fibonacci", vdw_radii="ProtOr")` implements
Shrake-Rupley, which is the algorithm TAP specifies.

## Why the antigen is out of scope

Every antigen-dependent question — which epitope, how tight, which clones compete — needs an
antibody–antigen complex, and predicting one is a different workflow with a different image and
a GPU. Folding that in would make this a GPU workflow and destroy its one real advantage, which
is that a few hundred clones cost CPU-minutes.

The honest division: this workflow tells you which clones are *developable*, workflow-012/-016/
-018 tell you what they *bind*. The README says so plainly, because a user who expects affinity
ranking from a developability tool will read the output wrongly.

## The ANARCI / ANARCII split

ANARCII is the better tool on paper — maintained this month, no HMMER, pure Python, Python
>=3.11, and `abnumber` already uses it. ANARCI classic was last touched in May 2024.

We use **classic ANARCI anyway**, because ImmuneBuilder imports it by name
(`from anarci import validate_sequence, anarci, scheme_short_to_long`) in the code path that
runs before every structure prediction. Installing ANARCII instead would leave ImmuneBuilder
unable to number anything, and running both would mean two numbering engines that can disagree
on CDR boundaries — worse than one slightly stale engine.

So: node 01 and node 02 both number with classic ANARCI, and the CDR boundaries the whole
pipeline uses come from one engine. Open item 4 tracks revisiting this if ImmuneBuilder migrates.

## TAP is a paper, not a package

TAP's guidelines, formulas and thresholds are published under CC BY; the implementation is a web
service with no source release. Two consequences the implementation respects:

1. **Node 03 is a reimplementation from the paper**, not a wrapper. Every output column, the
   README and the report say "TAP-style". Claiming TAP output would be a correctness claim we
   cannot support and a provenance claim about someone else's tool.
2. **Our numbers will not match the web server's exactly**, for the reasons in the next section.
   The verification plan therefore checks *agreement in flag direction* on shared antibodies,
   not numerical identity.

Every implementation choice node 03 makes is recorded in `developability.json` under
`implementation`, so a disagreement can be attributed rather than argued about.

## The calibration problem

Node 03 applies thresholds calibrated on TAP's pipeline to values produced by ours. Three known
differences:

1. **Exposure reference.** TAP: side-chain relative exposure vs an Ala-R-Ala reference from Yang
   Zhu & Blundell (1996). Ours: whole-residue relative ASA vs Tien et al.'s theoretical
   Ala-X-Ala values. Including backbone area in both numerator and denominator makes small
   residues look relatively more exposed, so our exposed set is not identical to theirs, and
   every metric is built on that set.
2. **Structure source.** TAP's thresholds come from ABodyBuilder **v1** models, and the paper
   notes modelling biases PSH upward. We model with ABodyBuilder **2**, a different (better)
   predictor with its own bias.
3. **Numerical details.** SASA point density, van der Waals radii set, and the `[1,2]`
   normalization of the hydrophobicity scale are not pinned by the paper.

The implementation's response is not to hide this: `tap_threshold_set = clinical` is the default
because flags are what makes the output useful, but node 03 records a `calibration` block in
`developability.json`, the report carries the caveat beside the flag legend, and the summary
compares the panel's PSH mean against the published CST mean of 123.30 ± 16.60 so a gross offset
is visible rather than silent. Open item 2 tracks closing the gap properly.

## The MHC-II hole

T-cell epitope prediction is the scientifically right way to assess immunogenicity risk, and
this workflow does not do it. The reason is licensing, not difficulty:

- **NetMHCIIpan** — the field standard — is academic-use-only. It cannot ship in a commercial
  platform's container.
- **mhcflurry** is Apache-2.0 and would be fine, but it predicts **class I** binding, which is
  the wrong class for antibody immunogenicity.

Rather than ship a class-I predictor and let it be read as immunogenicity, or ship a licence
violation, node 04 reports **humanness** (promb peptide content + germline identity) and says in
as many words that humanness is a proxy and not a T-cell epitope prediction. Open item 3 tracks
a licensed NetMHCIIpan as an opt-in node for customers who hold their own licence.

This is a selling point worth stating: the whole image is BSD-3/MIT/Apache-2.0, with no
academic-only component and no 22 GB download.

## Node structure

```text
panel CSV / FASTA (VH + VL per clone)
   │
   ▼
[01-annotate]        ANARCI numbering, CDRs, germline, liability motifs
   │                 → annotations.json, annotations.csv, liabilities.csv, clones.fasta
   ├──────────────────────────────┬────────────────────────────┐
   ▼                              ▼                            │
[02-model-fv]                [04-humanness]                    │
 ABodyBuilder2 / NanoBody     promb + germline identity        │
 → models/*.pdb              → humanness.csv                   │
   │                              │                            │
   ▼                              │                            │
[03-developability]               │                            │
 TAP-style 5 metrics + SASA       │                            │
 → developability.csv, flags.csv  │                            │
   │                              │                            │
   ▼                              │                            │
[05-mutate-scan]  AbLang2 proposals, re-screened ◄─────────────┘
 → proposals.csv
   │
   ▼
[06-rank]  join, flag, sort → ranked_panel.csv, panel_summary.json
   │
   ▼
[07-report] → report.html
```

Nodes 02 and 04 are independent of each other; the `[dependencies]` table in
`.chiral/workflow.toml` expresses the real edges rather than a straight line, so Silva may
schedule them as it likes. Node 04 needs no structure, which is deliberate: it still produces
output if modelling fails for some clones.

The diagram shows where each measurement originates, not the full edge set: nodes 06 and 07 also
read upstream detail files directly (06 for model error and proposal counts, 07 for per-residue
exposure, the calibration block and the re-screen rules), so both list every upstream node whose
output they consume. A declared input that no dependency produces is a staging failure at run
time, so the two are checked against each other.

| Node | Cost driver | Scales as |
|---|---|---|
| 01-annotate | ANARCI HMM search | clones × chains, milliseconds each |
| 02-model-fv | ImmuneBuilder + OpenMM refinement | clones, ~seconds each; refinement dominates |
| 03-developability | SASA + pairwise distances per model | clones, fast |
| 04-humanness | exact 9-mer lookups | clones, very fast (nearest-peptide search is the slow path) |
| 05-mutate-scan | AbLang2 forward passes | clones (one pass), or clones × positions with stepwise masking |
| 06-rank / 07-report | table joins, HTML | negligible |

## Conventions followed

Taken from workflow-037 and -038, which are the current house style:

- **Every declared output is always written.** A declared-but-unwritten output is a silent
  failure in Silva. Node 05 with `mutate_scan = false` writes an empty `proposals.csv` and a
  manifest with `enabled: false`; nodes 06/07 read that flag and omit the section rather than
  presenting an empty table as a result.
- **A manifest per node that makes a decision**, recording resolved parameters, tool versions
  and anything skipped.
- **Weights and databases baked into the image**, every node offline.
- **`--dry-run` on the expensive nodes** (02 and 05): resolve and print every setting, touch no
  model.
- **Per-clone failures are recorded, not fatal.** One unnumberable sequence in a 500-clone panel
  must not kill the run; it gets a row with a reason and the panel continues. A node fails hard
  only when *every* clone fails, or when an input is missing.
- **Nothing is filtered.** Flags are columns.
- **The report leads with the caveats** rather than burying them.

## Open items

1. **Conda install resolvability (blocking the first build).** The image installs `anarci`,
   `hmmer`, `openmm` and `pdbfixer` from conda-forge/bioconda at python=3.11. Untested. If
   bioconda's `anarci` lags 3.11, pin python=3.10 and check promb still installs, or install
   ANARCI from the GitHub tarball over a conda HMMER.
2. **Close the calibration gap.** Either obtain a side-chain Ala-R-Ala exposure reference and
   switch the numerator to side-chain atoms, or re-derive thresholds from a CST set modelled
   with ABodyBuilder2. Until then the flags carry a documented offset of unknown size. See
   [the calibration problem](#the-calibration-problem).
3. **MHC-II as an opt-in node.** Design a node 08 that activates only when the operator mounts
   their own licensed NetMHCIIpan. Needs a story for how Silva passes a licensed binary in.
4. **ANARCII migration.** Revisit when ImmuneBuilder either migrates or gains a pluggable
   numbering backend.
5. **Aggregation propensity.** Deliberately absent. The credible open tools (Aggrescan3D,
   CamSol) are web services or licence-restricted; a hand-rolled hydrophobic-patch score called
   "aggregation propensity" would overclaim. PSH is the honest partial answer.
6. **Expression/titre prediction.** Out of scope, no credible open tool.
7. **scFv linker detection.** Node 01 rejects scFvs rather than splitting them. A linker
   detector is a plausible v2 convenience and a plausible source of silent mis-splits.
8. **Reconcile with `workflow-abodybuilder3`.** That unmerged branch predicts antibody Fv
   structures with ABodyBuilder3; node 02 here does the same with ABodyBuilder2. If it lands,
   decide which modeller this workflow consumes rather than shipping two. ABodyBuilder3 is the
   newer model; the reason to keep ImmuneBuilder is that it covers nanobodies in the same call
   and is already a dependency of the numbering path.
9. **Nanobody thresholds.** The TAP thresholds are derived from paired VH/VL therapeutics.
   Node 03 computes the metrics for nanobodies but SFvCSP is undefined (no light chain) and the
   remaining bounds are not calibrated for VHH. Nanobody rows are flagged
   `thresholds_not_calibrated` rather than silently compared.

## Verification plan

**Done:** every script byte-compiles, and `tests/test_logic.py` covers the dependency-free
logic — CDR/anchor derivation from the border tables, liability motif scanning (including the
unpaired-cysteine rule and position reporting), hydrophobicity normalization to [1,2], the
charge table, PSH/PPC/PNC summation on a hand-built synthetic geometry with a known answer,
salt-bridge detection, flag boundary conditions at every Table 2 threshold, and panel ranking
order. Run with `python3 tests/test_logic.py`.

**Outstanding — nothing below has been run:**

1. **Image builds.** `docker build` completes; the four bake steps print their ready messages,
   including the ANARCI numbering smoke test.
2. **Versions pinned.** `pip list` shows ImmuneBuilder 1.2, ablang2 0.2.1, promb 1.0.2;
   `conda list` shows anarci, openmm, pdbfixer. Record resolved versions in the README.
3. **Node 04 against a published value.** `compute_peptide_content` on promb's documented VH
   returns 0.7142857142857143. The only check with an externally fixed expected number.
4. **Node 02 against the published benchmark.** Model antibodies with solved structures and
   confirm CDR-H3 RMSD lands near upstream's reported 2.81 Å average. Catches a broken weights
   bake rather than validating their model.
5. **Node 01 numbering sanity.** Known VH/VL pairs get the expected chain types and CDR
   boundaries under IMGT; a truncated sequence and a full IgG are both rejected with a reason.
6. **Node 03 against the CST distribution.** Run the bundled panel and compare CDR-vicinity PSH
   against the published CST mean of 123.30 ± 16.60. A mean far outside that range means the
   exposure reference or the normalization is wrong, not that the antibodies are unusual.
7. **Node 03 flag agreement.** Pick antibodies with published TAP web-server output and compare
   **flag direction**, not values. Amber-band disagreement is expected and must be documented;
   green-vs-red disagreement is a bug.
8. **Liability scan against hand annotation.** A sequence with a planted N-glycosylation sequon,
   an NG and a free cysteine is flagged for exactly those and nothing else. (Covered
   synthetically in the unit tests; still worth doing against a real sequence.)
9. **Node 05 re-screening actually screens.** Construct a case where the top proposal introduces
   an NG motif and confirm it is dropped with a recorded reason. (Unit-tested synthetically.)
10. **Panel-scale run.** A few hundred clones end to end on CPU; record wall clock per node to
    substantiate the cost claim in the README. The claim must be deleted or corrected if the
    measurement disagrees.
11. **Failure isolation.** Inject one unnumberable sequence into a panel and confirm the run
    completes with that clone marked failed.

## Scientific caveats this workflow must not paper over

- **Developability flags are risk indicators, not predictions of clinical success.** The TAP
  thresholds describe where clinical-stage antibodies sit; sitting outside them is a reason to
  look, not a reason to kill a clone. Approved therapeutics violate these guidelines.
- **Our metrics are TAP-style, not TAP**, and carry a calibration offset of unknown size.
- **Humanness is not immunogenicity.** No T-cell epitope prediction is performed.
- **Predicted structures carry error**, and CDR-H3 — the loop that matters most and is modelled
  worst — carries the most. The per-residue error estimate is in the B-factors and the report
  shows it next to any structure-derived metric.
- **AbLang2 proposals are likelihood-based, not affinity-based.** A high-probability mutant is
  one that looks like a natural antibody, which is not the same as one that binds better. They
  are bench hypotheses.
- **No antigen is considered anywhere**, so nothing here speaks to potency, specificity or
  epitope coverage.
- **Germline bias.** Both the language model and the humanness reference are built from observed
  repertoires, so unusual-but-fine clones score worse than common-but-mediocre ones.
- **Nanobody metrics are uncalibrated** — see open item 9.

## References

- [Raybould et al., *Five computational developability guidelines for therapeutic antibody
  profiling*, PNAS 2019](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6410772/) — TAP's
  guidelines, formulas and Table 2 thresholds. CC BY. DOI `10.1073/pnas.1810576116`.
- [Tien et al., *Maximum allowed solvent accessibilites of residues in proteins*, PLoS ONE
  2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3836772/) — the Ala-X-Ala exposure
  reference.
- [brennanaba/ImmuneBuilder](https://github.com/brennanaba/ImmuneBuilder) — ABodyBuilder2 /
  NanoBodyBuilder2, BSD-3.
- [oxpig/ANARCI](https://github.com/oxpig/ANARCI) — numbering, BSD-3.
- [oxpig/ANARCII](https://github.com/oxpig/ANARCII) — the successor, BSD-3; see open item 4.
- [prihoda/abnumber](https://github.com/prihoda/abnumber) — the CDR border tables, MIT.
- [oxpig/AbLang2](https://github.com/oxpig/AbLang2) — paired antibody language model, BSD-3.
- [MSDLLCPapers/promb](https://github.com/MSDLLCPapers/promb) — human peptide prevalence, MIT.
- [Merck/BioPhi](https://github.com/Merck/BioPhi) — Sapiens humanization and OASis, MIT; not
  used here, see the scope table.
- [README.md](README.md) — the user-facing workflow documentation.
