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
import re
import subprocess
import sys
from pathlib import Path

STRUCTURE_EXT = (".cif", ".mmcif")
FASTA_EXT = (".fasta", ".fa", ".faa")
YAML_EXT = (".yaml", ".yml")

VALID_AA = set("ACDEFGHIKLMNPQRSTVWYX")

# Searched in this order so an explicit system YAML always beats a structure file.
SEARCH_ORDER = YAML_EXT + STRUCTURE_EXT + (".pdb",) + FASTA_EXT


# Codes a depositor uses for density that could not be identified. Not chemistry.
PLACEHOLDER_COMPONENTS = {"UNX", "UNK", "UNL", "DUM"}

# Solvent that upstream does not recognise as solvent. `System.from_cif` omits waters, but it
# omits HOH: a neutron structure writes its waters as DOD, deuterated water, and those come
# through as ligands and are sampled as ligand tokens. 5PTI is a neutron structure and carries
# one. The difference matters more the more of them an entry has.
SOLVENT_COMPONENTS = {"DOD"}


# Anything that is not a standard residue — an ion, a cofactor, a terminal cap, a modified
# residue — is built from its CCD definition, which bioemu2 downloads from files.rcsb.org the
# first time it meets the code and caches under BIOEMU2_CCD_CACHE / BIOEMU2_TEMPLATE_CACHE.
# The image pre-warms the common codes, so an air-gapped run works for those and only those.
# Without this explanation the failure surfaces as a DNS error raised four frames inside the
# library, which reads as a broken image rather than as a missing chemical component.
def ccd_advice(exc):
    """The sentence that turns a connection error into an actionable one, or '' if unrelated."""
    text = str(exc)
    if "ligands/download" not in text and "rcsb.org" not in text:
        return ""
    match = re.search(r"ligands/download/([A-Za-z0-9]+)\.cif", text)
    code = match.group(1) if match else None
    return (
        " This is bioemu2 fetching the CCD definition for "
        + (f"component '{code}'" if code else "a non-standard component")
        + " from files.rcsb.org, which it does the first time it meets a code that is not a "
        "standard residue. The image pre-warms the common ions, caps, cofactors and modified "
        "residues into BIOEMU2_CCD_CACHE and BIOEMU2_TEMPLATE_CACHE; this code is not among "
        "them, so this node needs network access, or those cache directories have to be "
        "supplied already holding it."
    )


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


def load_system(upload, pdb_id, sequence, include_ligands, exclude_ligands):
    """Return (System, source_description).

    `include_ligands` and `exclude_ligands` are upstream's own arguments to `System.from_cif`
    and `System.from_pdb_id`, surfaced because the default keeps every non-water component in
    the file -- cryoprotectants and buffer molecules included. 1CLL, for example, carries an
    ethanol: by default BioEmu2 samples it along with the protein and its four calcium ions.
    That is the right default (dropping a component silently is worse than sampling one), but
    it has to be visible and it has to be refusable, which is what these two give.

    A system YAML, a FASTA or a bare sequence declares its own components, so the two
    arguments do not apply there and are reported as ignored rather than silently dropped.
    """
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
                return (
                    System.from_cif(
                        Path(upload),
                        include_ligands=include_ligands,
                        exclude_ligands=exclude_ligands or None,
                    ),
                    f"structure: {os.path.basename(upload)}",
                )
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
        print(
            f"Fetching RCSB entry {pdb_id.upper()} (assembly 1, "
            + ("ligands included" if include_ligands else "ligands excluded")
            + (f", excluding {exclude_ligands}" if exclude_ligands else "")
            + ")",
            flush=True,
        )
        try:
            return (
                System.from_pdb_id(
                    pdb_id.upper(),
                    include_ligands=include_ligands,
                    exclude_ligands=exclude_ligands or None,
                ),
                f"RCSB:{pdb_id.upper()}",
            )
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
        # Which residues that counts, recorded because the two plausible readings differ on
        # most real entries: this is the length of the entity sequence BioEmu2 will model, not
        # the number of residues resolved in the input file. 1CLL is a 148-residue protein with
        # 144 resolved, and the ensemble BioEmu2 produces has all 148. Reported as the second
        # number, a report says "148 residues" about a 144-residue crystal structure.
        "n_polymer_residues_from": "entity sequence (includes residues not resolved in the input)",
        "n_tokens": None,
        "n_tokens_from": None,
        "has_ligands": False,
        "ligand_components": [],
        "has_nucleic_acid": False,
        "protein_chains": [],
    }

    # `num_tokens` is the field name on AllAtomStructure (bioemu2/structure/all_atom_structure.py:
    # `num_tokens: int  # Number of tokens (different from num_residues in presence of ligands)`).
    # This read was `n_tokens` and so always returned None, which left the token total null in the
    # inventory, the >10k guard below unreachable, and the README claiming a check that never ran.
    # The attribute that answered is recorded, so a renamed field shows up as a changed provenance
    # string rather than as a silently absent number.
    try:
        from bioemu2.structure import AllAtomStructure

        structure = AllAtomStructure.from_system(system)
        for attribute in ("num_tokens", "n_tokens"):
            value = getattr(structure, attribute, None)
            if value is not None:
                inventory["n_tokens"] = int(value)
                inventory["n_tokens_from"] = f"AllAtomStructure.{attribute}"
                break
        else:
            print(
                "WARNING: AllAtomStructure exposed neither num_tokens nor n_tokens, so the "
                "token total is unknown and the GPU-fit warning cannot be applied. This is a "
                "bioemu2 API change, not a problem with the input.",
                flush=True,
            )
    except Exception as exc:
        # Not a warning. Nodes 03, 04 and 05 build this same structure from this same system,
        # so a failure here is a failure of the whole run, and the only thing continuing buys
        # is that it happens later, on a GPU, after the model has loaded.
        #
        # Naming the non-polymer components is the useful part: whatever the immediate error,
        # the component that cannot be resolved or modelled is almost always one of them, and
        # dropping it is a parameter away. 5PTI carries a UNX — "unknown atom or ion" — and
        # UNL, the unknown-ligand placeholder, raises inside bioemu2's own template builder,
        # so an entry can be unmodellable for reasons that have nothing to do with the network.
        components = sorted(
            {
                code
                for code in (
                    getattr(entity, "ccd_code", None)
                    for entity in (getattr(system, "entities", []) or [])
                )
                if code
            }
        )
        fail(
            f"Could not build the all-atom structure for this system: {exc}."
            f"{ccd_advice(exc)}"
            + (
                f" The non-polymer components in this system are {components}. If one of them "
                "cannot be modelled, drop it with exclude_ligands (e.g. "
                f"exclude_ligands={components[0]}), or set include_ligands=false to sample the "
                "polymers alone."
                if components
                else ""
            )
            + " Nodes 03, 04 and 05 construct the same structure, so this run cannot proceed "
            "past it."
        )

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
        ccd_code = getattr(entity, "ccd_code", None)
        if "ligand" in entity_type or "branched" in entity_type:
            inventory["has_ligands"] = True
            # Named, not just counted. BioEmu2 samples every one of these as ligand tokens, so
            # the difference between four calcium ions and one ethanol is the difference between
            # the holo protein and the crystallisation buffer -- and the report can only say
            # which if this node writes down what it kept.
            for asym_id in asym_ids:
                inventory["ligand_components"].append(
                    {"asym_id": asym_id, "ccd_code": ccd_code}
                )

        for asym_id in asym_ids:
            inventory["chains"].append(
                {
                    "asym_id": asym_id,
                    "entity_type": entity_type,
                    "length": length,
                    "sequence": sequence if isinstance(sequence, str) else None,
                    "ccd_code": ccd_code,
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
    include_ligands = os.environ.get("PARAM_INCLUDE_LIGANDS", "true").lower() == "true"
    exclude_ligands = [
        code.strip().upper()
        for code in os.environ.get("PARAM_EXCLUDE_LIGANDS", "").replace(",", " ").split()
        if code.strip()
    ]
    # Upstream raises on the combination rather than ignoring one of them; say so here, where
    # the two parameters are visible, instead of surfacing it from inside from_cif.
    if exclude_ligands and not include_ligands:
        fail(
            f"exclude_ligands={exclude_ligands} cannot be combined with include_ligands=false. "
            "include_ligands=false already drops every non-water component; exclude_ligands is "
            "for dropping named ones while keeping the rest."
        )

    upload = find_upload()
    if upload and (pdb_id or sequence):
        print(
            f"NOTE: uploaded file '{upload}' takes precedence over pdb_id/sequence.",
            flush=True,
        )

    system, source = load_system(upload, pdb_id, sequence, include_ligands, exclude_ligands)

    if (not include_ligands or exclude_ligands) and not source.startswith(
        ("structure:", "RCSB:")
    ):
        print(
            f"NOTE: include_ligands/exclude_ligands apply to a structure or a PDB ID. This "
            f"system came from {source}, which declares its own components, so they were not "
            "applied.",
            flush=True,
        )

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

    placeholders = sorted(
        {
            item["ccd_code"]
            for item in inventory["ligand_components"]
            if item["ccd_code"] in PLACEHOLDER_COMPONENTS
        }
    )
    if placeholders:
        print(
            f"WARNING: this system contains placeholder component(s) {placeholders}. These are "
            "the codes a depositor uses for density that could not be identified — UNX is an "
            "unknown atom or ion, UNK an unknown residue, UNL an unknown ligand, DUM a dummy "
            "atom — so BioEmu2 will sample a molecule whose chemistry nobody has established. "
            f"Drop them with exclude_ligands={','.join(placeholders)} unless you specifically "
            "want them modelled.",
            flush=True,
        )

    solvent = sorted(
        {
            item["ccd_code"]
            for item in inventory["ligand_components"]
            if item["ccd_code"] in SOLVENT_COMPONENTS
        }
    )
    if solvent:
        count = sum(
            1 for item in inventory["ligand_components"] if item["ccd_code"] in SOLVENT_COMPONENTS
        )
        print(
            f"WARNING: this system contains {count} solvent component(s) {solvent}. BioEmu2 omits "
            "waters written as HOH, but a neutron structure writes them as DOD — deuterated "
            "water — and those are kept and sampled as ligand tokens. That count is of "
            "components, not of molecules: 5PTI's 63 D2O molecules arrive here as one DOD "
            "component worth one token, so this is a stray component in the system rather than "
            "a solvent shell around it. Drop them with "
            f"exclude_ligands={','.join(solvent)} if you did not mean to sample solvent.",
            flush=True,
        )

    if not inventory["protein_chains"]:
        print(
            "WARNING: no protein chain in this system. BioEmu2 computes no MSA for a "
            "non-protein system and performs no nucleic-acid or ligand physicality checks.",
            flush=True,
        )

    print(f"System source: {source}", flush=True)
    print(f"Wrote {out_yaml}", flush=True)
    for chain in inventory["chains"]:
        if chain["ccd_code"]:
            print(
                f"  chain {chain['asym_id']}: {chain['entity_type']} {chain['ccd_code']}",
                flush=True,
            )
            continue
        length = chain["length"] if chain["length"] is not None else "?"
        print(
            f"  chain {chain['asym_id']}: {chain['entity_type']}, {length} residues",
            flush=True,
        )
    tokens = inventory["n_tokens"]
    print(
        f"{inventory['n_chains']} chain(s), {inventory['n_polymer_residues']} polymer "
        f"residues in the entity sequence (which is what gets modelled, so it includes any "
        f"residue the input structure does not resolve), "
        f"{tokens if tokens is not None else 'unknown'} tokens "
        f"(from {inventory['n_tokens_from'] or 'nothing — unavailable'})",
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
