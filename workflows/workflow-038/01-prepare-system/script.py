#!/usr/bin/env python3
"""Node 01 — Prepare System.

Resolve the molecular system to sample and write it as a canonical BioEmu2 system YAML,
plus system_inventory.json describing what is in it.

Resolution order (first match wins):
  1. a file in ./inputs/ — .yaml/.yml (BioEmu2 system YAML, or Boltz YAML, auto-detected),
     .cif/.mmcif, .pdb, or .fasta/.fa/.faa
  2. PARAM_PDB_ID  — fetched by System.from_pdb_id() as biological assembly 1, ligands included
  3. PARAM_SEQUENCE — a single protein chain A

Everything is funnelled through `System` and written back out with `system.to_yaml()`, so
node 03 always receives one canonical format no matter what the user supplied, and any
malformed input fails here, on CPU, in seconds.

The inventory exists because BioEmu2's cost and memory both scale with **token** count, not
residue count: one token per polymer residue plus one per ligand atom. The token total is the
number that predicts whether a run will fit on the GPU, so it is computed once, here, and
carried forward into every manifest and the report.
"""
import glob
import json
import os
import subprocess
import sys
from pathlib import Path

STRUCTURE_EXT = (".cif", ".mmcif")
FASTA_EXT = (".fasta", ".fa", ".faa")
YAML_EXT = (".yaml", ".yml")

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYX")

# Searched in this order so an explicit system YAML always beats a structure file.
SEARCH_ORDER = YAML_EXT + STRUCTURE_EXT + (".pdb",) + FASTA_EXT


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def find_upload():
    for ext in SEARCH_ORDER:
        matches = sorted(glob.glob(f"./inputs/*{ext}"))
        if matches:
            return matches[0]
    return None


def is_boltz_yaml(path):
    """Boltz input YAML uses a top-level `sequences:` list; BioEmu2 uses `entities:`."""
    try:
        import yaml

        with open(path) as handle:
            loaded = yaml.safe_load(handle)
    except Exception as exc:
        fail(f"Could not parse YAML '{path}': {exc}")
    if not isinstance(loaded, dict):
        fail(f"'{path}' is not a YAML mapping; expected a BioEmu2 system YAML.")
    if "entities" in loaded:
        return False
    if "sequences" in loaded:
        return True
    fail(
        f"'{path}' has neither an `entities:` key (BioEmu2 system YAML) nor a `sequences:` "
        "key (Boltz YAML). See input_files/README.md for the accepted formats."
    )


def convert_boltz(path):
    """Delegate to `bioemu2 convert-boltz`, which also reports fields it cannot carry over."""
    converted = "./_converted_system.yaml"
    argv = ["bioemu2", "convert-boltz", path, "--out", converted]
    print(f"Boltz YAML detected; converting: {argv}", flush=True)
    result = subprocess.run(argv)
    if result.returncode != 0:
        fail(f"bioemu2 convert-boltz exited {result.returncode} (output above)")
    return converted


def pdb_to_mmcif(src):
    dest = "./_converted.cif"
    from Bio.PDB import MMCIFIO, PDBParser

    print(f"Converting PDB to mmCIF: {src}", flush=True)
    parser = PDBParser(QUIET=True)
    try:
        structure = parser.get_structure("system", src)
    except Exception as exc:
        fail(f"Could not parse PDB '{src}': {exc}")
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(dest)
    return dest


def system_from_fasta(path):
    """One entity per record, chain IDs assigned A, B, C... in file order."""
    from Bio import SeqIO

    from bioemu2 import System

    records = list(SeqIO.parse(path, "fasta"))
    if not records:
        fail(f"No sequences found in '{path}'")

    sequences = []
    for record in records:
        sequence = str(record.seq).upper().replace("*", "")
        if not sequence:
            fail(f"Record '{record.id}' in '{path}' is empty")
        bad = sorted(set(sequence) - VALID_AA)
        if bad:
            fail(
                f"Record '{record.id}' contains non-amino-acid characters {bad}. "
                "Use X for an unidentified residue."
            )
        sequences.append(sequence)

    print(f"Building system from {len(sequences)} FASTA record(s)", flush=True)
    return System.from_description(proteins=sequences)


def system_from_sequence(sequence):
    from bioemu2 import System

    sequence = sequence.strip().upper().replace(" ", "").replace("*", "")
    bad = sorted(set(sequence) - VALID_AA)
    if bad:
        fail(
            f"sequence contains non-amino-acid characters {bad}. "
            "Use X for an unidentified residue."
        )
    print(f"Building system from a {len(sequence)}-residue sequence parameter", flush=True)
    return System.from_description(proteins=[sequence])


def load_system(upload, pdb_id, sequence):
    """Return (System, source_description)."""
    from bioemu2 import System

    if upload:
        lower = upload.lower()
        if lower.endswith(YAML_EXT):
            path = convert_boltz(upload) if is_boltz_yaml(upload) else upload
            print(f"Loading system YAML: {path}", flush=True)
            try:
                return System.from_yaml(Path(path)), f"system YAML: {os.path.basename(upload)}"
            except Exception as exc:
                fail(f"Could not load system YAML '{path}': {exc}")
        if lower.endswith(".pdb"):
            upload = pdb_to_mmcif(upload)
            lower = upload.lower()
        if lower.endswith(STRUCTURE_EXT):
            print(f"Loading structure: {upload}", flush=True)
            try:
                return System.from_cif(Path(upload)), f"structure: {os.path.basename(upload)}"
            except Exception as exc:
                fail(
                    f"Could not import '{upload}': {exc}. System.from_cif() reads molecular "
                    "composition from mmCIF metadata (_entity_poly and friends), not from "
                    "coordinates alone, so a coordinates-only CIF cannot be imported. "
                    "Use an RCSB file, or write a system YAML instead."
                )
        if lower.endswith(FASTA_EXT):
            return system_from_fasta(upload), f"FASTA: {os.path.basename(upload)}"
        fail(f"Unsupported input file type: {upload}")

    if pdb_id:
        print(f"Fetching RCSB entry {pdb_id.upper()} (assembly 1, ligands included)", flush=True)
        try:
            return System.from_pdb_id(pdb_id.upper()), f"RCSB:{pdb_id.upper()}"
        except Exception as exc:
            fail(f"Could not load PDB ID '{pdb_id}': {exc}")

    if sequence:
        return system_from_sequence(sequence), "sequence parameter"

    fail(
        "No system. Upload a system YAML, structure or FASTA file, or set pdb_id, "
        "or set sequence."
    )


def build_inventory(system, source):
    """Per-chain inventory plus the token total that predicts GPU fit.

    One token per polymer residue, one per ligand atom. Chain inspection goes through the
    public System API; anything it does not expose is reported as null rather than guessed.
    """
    inventory = {
        "source": source,
        "chain_ids": list(getattr(system, "chain_ids", []) or []),
        "chains": [],
        "n_chains": 0,
        "n_polymer_residues": 0,
        "n_tokens": None,
        "has_ligands": False,
        "has_nucleic_acid": False,
        "protein_chains": [],
    }

    try:
        from bioemu2.structure import AllAtomStructure

        structure = AllAtomStructure.from_system(system)
        n_tokens = getattr(structure, "n_tokens", None)
        if n_tokens is None:
            tokens = getattr(structure, "tokens", None)
            n_tokens = len(tokens) if tokens is not None else None
        if n_tokens is not None:
            inventory["n_tokens"] = int(n_tokens)
    except Exception as exc:
        print(f"WARNING: could not compute token count: {exc}", flush=True)

    for entity in getattr(system, "entities", []) or []:
        entity_type = getattr(entity, "entity_type", None) or type(entity).__name__.lower()
        entity_type = str(entity_type).lower()
        asym_ids = list(getattr(entity, "asym_ids", []) or [])
        sequence = getattr(entity, "sequence", None)
        length = len(sequence) if isinstance(sequence, str) else None

        if "protein" in entity_type:
            inventory["protein_chains"].extend(asym_ids)
        if "rna" in entity_type or "dna" in entity_type:
            inventory["has_nucleic_acid"] = True
        if "ligand" in entity_type or "branched" in entity_type:
            inventory["has_ligands"] = True

        for asym_id in asym_ids:
            inventory["chains"].append(
                {
                    "asym_id": asym_id,
                    "entity_type": entity_type,
                    "length": length,
                    "sequence": sequence if isinstance(sequence, str) else None,
                }
            )
            if length:
                inventory["n_polymer_residues"] += length

    inventory["n_chains"] = len(inventory["chains"])
    if not inventory["chain_ids"]:
        inventory["chain_ids"] = [chain["asym_id"] for chain in inventory["chains"]]
    return inventory


def main():
    os.makedirs("./outputs", exist_ok=True)

    pdb_id = os.environ.get("PARAM_PDB_ID", "")
    sequence = os.environ.get("PARAM_SEQUENCE", "")

    upload = find_upload()
    if upload and (pdb_id or sequence):
        print(
            f"NOTE: uploaded file '{upload}' takes precedence over pdb_id/sequence.",
            flush=True,
        )

    system, source = load_system(upload, pdb_id, sequence)

    out_yaml = "./outputs/system.yaml"
    try:
        system.to_yaml(Path(out_yaml))
    except Exception as exc:
        fail(f"Could not serialise the system to YAML: {exc}")

    # Round-trip as a validation step: if node 03 cannot read this back, fail here instead.
    try:
        from bioemu2 import System

        System.from_yaml(Path(out_yaml))
    except Exception as exc:
        fail(f"The system YAML written to {out_yaml} does not load back: {exc}")

    inventory = build_inventory(system, source)
    with open("./outputs/system_inventory.json", "w") as handle:
        json.dump(inventory, handle, indent=2)

    if not inventory["protein_chains"]:
        print(
            "WARNING: no protein chain in this system. BioEmu2 computes no MSA for a "
            "non-protein system and performs no nucleic-acid or ligand physicality checks.",
            flush=True,
        )

    print(f"System source: {source}", flush=True)
    print(f"Wrote {out_yaml}", flush=True)
    for chain in inventory["chains"]:
        length = chain["length"] if chain["length"] is not None else "?"
        print(
            f"  chain {chain['asym_id']}: {chain['entity_type']}, {length} residues",
            flush=True,
        )
    tokens = inventory["n_tokens"]
    print(
        f"{inventory['n_chains']} chain(s), {inventory['n_polymer_residues']} polymer "
        f"residues, {tokens if tokens is not None else 'unknown'} tokens",
        flush=True,
    )
    if isinstance(tokens, int) and tokens > 10000:
        print(
            "WARNING: over 10k tokens. Upstream notes that physically plausible structures "
            "may be hard to produce above this size, and it may not fit on one GPU; "
            "context-parallel sampling is out of scope for this workflow.",
            flush=True,
        )


if __name__ == "__main__":
    main()
