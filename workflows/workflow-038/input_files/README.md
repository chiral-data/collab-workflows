# Sample input files

## trp_cage_system.yaml (default)

Trp-cage, PDB [1L2Y](https://www.rcsb.org/structure/1L2Y) — a 20-residue miniprotein given as a
BioEmu2 **system YAML**. This is the default `input_file` and the workflow's smoke test: it
samples in well under a minute on any supported GPU.

A 20-residue query returns a shallow alignment, and BioEmu2's ensemble quality tracks MSA depth,
so treat this as a plumbing test rather than a scientific result.

## Supplying your own input

Node 01 accepts any one of the following, in this order of precedence:

1. **A file in the workflow's input directory**, detected by extension:
   - `.yaml` / `.yml` — a BioEmu2 system YAML, used as-is (validated by round-tripping through
     `System.from_yaml()`). A Boltz-format YAML is detected and converted with
     `bioemu2 convert-boltz`.
   - `.cif` / `.mmcif` — an mmCIF structure, read with `System.from_cif()`. Needs the
     `_entity_poly` metadata that RCSB files carry; a coordinates-only CIF will not import.
   - `.pdb` — converted to mmCIF first, then read as above.
   - `.fasta` / `.fa` / `.faa` — one entity per record, chain IDs assigned A, B, C… in file
     order.
2. **`pdb_id`** — an RCSB accession, fetched as biological assembly 1 including all non-water
   ligands. The CLI exposes no assembly, chain or ligand selection; for that, build a system
   YAML instead.
3. **`sequence`** — a raw one-letter protein sequence, single chain A.

## Suggested second target

`pdb_id = "1AKE"` — adenylate kinase, upstream's own worked example and a genuine two-state
system. An equilibrium ensemble should contain both the open and the closed lobe conformation;
if it shows only one, sampling is broken in a way Trp-cage will not reveal. It is ~214 residues
per chain, so expect minutes rather than seconds, and a real MSA.
