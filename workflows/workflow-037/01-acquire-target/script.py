#!/usr/bin/env python3
"""Node 01 — Acquire Target.

Resolve the design target: prefer a structure uploaded into ./inputs/ (.cif, .pdb or
.cif.zst), otherwise fetch PARAM_PDB_ID from RCSB. Normalise to mmCIF as target.cif and
emit target_inventory.json.

The inventory carries BOTH numberings per chain. RFD4 contigs and sel() selectors use the
original mmCIF label_seq_id, while RCSB files are author-numbered, so a hotspot copied
from a paper ("A45") is usually author numbering and needs translating in node 02.

An uploaded mmCIF is passed through byte-for-byte rather than re-serialised: rewriting it
would renumber label_seq_id and silently invalidate every contig and selector downstream.
Only PDB input is converted, and conversion assigns label_seq_id sequentially per chain,
which is what any converter does and what RFD4 will then see.
"""
import glob
import json
import os
import shutil
import sys

RCSB_URL = "https://files.rcsb.org/download/{pdb_id}.cif"

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def resolve_upload():
    uploads = sorted(
        glob.glob("./inputs/*.cif.zst")
        + glob.glob("./inputs/*.cif")
        + glob.glob("./inputs/*.pdb")
    )
    return uploads[0] if uploads else None


def fetch_from_rcsb(pdb_id, dest):
    import requests

    url = RCSB_URL.format(pdb_id=pdb_id.upper())
    print(f"Fetching {url}", flush=True)
    try:
        response = requests.get(url, timeout=120)
    except Exception as exc:
        print(f"ERROR: RCSB fetch failed for '{pdb_id}': {exc}", flush=True)
        sys.exit(1)
    if response.status_code != 200:
        print(
            f"ERROR: RCSB returned HTTP {response.status_code} for '{pdb_id}'. "
            "Check the accession.",
            flush=True,
        )
        sys.exit(1)
    with open(dest, "wb") as handle:
        handle.write(response.content)


def decompress_zst(src, dest):
    try:
        import zstandard
    except ImportError:
        print(
            "ERROR: zstandard is not installed, cannot read a .cif.zst input. "
            "Decompress it before uploading.",
            flush=True,
        )
        sys.exit(1)
    with open(src, "rb") as fin, open(dest, "wb") as fout:
        zstandard.ZstdDecompressor().copy_stream(fin, fout)


def pdb_to_mmcif(src, dest):
    """Convert PDB to mmCIF, then canonicalise the label numbering it invents."""
    from Bio.PDB import MMCIFIO, PDBParser

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("target", src)
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(dest)
    canonicalise_label_numbering(dest)


def canonicalise_label_numbering(path):
    """Give every author chain exactly one label_asym_id and a sequential label_seq_id.

    Biopython's mmCIF writer starts a new label_asym_id at each non-standard residue, so a
    peptide containing one modified residue comes out as several one- and two-residue
    "chains" (5O45 chain B becomes 14 of them). Node 02 would then resolve chain B to a
    one-residue fragment and design against it without complaint. PDB input carries no
    label numbering of its own, so we are free to define it — and define it per author
    chain, which is what the contigs built downstream assume.

    Only files written by pdb_to_mmcif reach here: a single data block holding one
    _atom_site loop, which is why it is safe to re-emit the file wholesale.
    """
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    cif = MMCIF2Dict(path)
    columns = [key for key in cif if key.startswith("_atom_site.")]
    if not columns:
        fail(f"Converted mmCIF '{path}' has no _atom_site loop.")

    n_atoms = len(cif[columns[0]])
    auth_asym = cif.get("_atom_site.auth_asym_id", ["A"] * n_atoms)
    auth_seq = cif.get("_atom_site.auth_seq_id", ["."] * n_atoms)
    ins_code = cif.get("_atom_site.pdbx_PDB_ins_code", ["?"] * n_atoms)
    comp = cif.get("_atom_site.label_comp_id", ["UNK"] * n_atoms)
    group = cif.get("_atom_site.group_PDB", ["ATOM"] * n_atoms)

    counters = {}
    assigned = {}
    new_asym = []
    new_seq = []
    for i in range(n_atoms):
        chain = auth_asym[i]
        new_asym.append(chain)
        is_polymer = group[i] == "ATOM" or comp[i] in THREE_TO_ONE
        if not is_polymer:
            new_seq.append(".")  # ligands and waters carry no label_seq_id
            continue
        key = (chain, auth_seq[i], ins_code[i])
        if key not in assigned:
            counters[chain] = counters.get(chain, 0) + 1
            assigned[key] = counters[chain]
        new_seq.append(str(assigned[key]))

    cif["_atom_site.label_asym_id"] = new_asym
    cif["_atom_site.label_seq_id"] = new_seq

    def quote(value):
        text = "" if value is None else str(value)
        if text == "":
            return "."
        if any(ch.isspace() for ch in text) or text.startswith("_"):
            return f"'{text}'"
        return text

    with open(path, "w") as handle:
        handle.write("data_target\n#\nloop_\n")
        for key in columns:
            handle.write(f"{key}\n")
        for i in range(n_atoms):
            handle.write(" ".join(quote(cif[key][i]) for key in columns) + "\n")
        handle.write("#\n")


def normalise(src, dest):
    lower = src.lower()
    if lower.endswith(".cif.zst"):
        print(f"Decompressing {src}", flush=True)
        decompress_zst(src, dest)
    elif lower.endswith(".pdb"):
        print(f"Converting PDB to mmCIF: {src}", flush=True)
        print(
            "NOTE: PDB input has no label_seq_id; conversion assigns it sequentially "
            "per chain over standard residues only, so modified residues are skipped and "
            "the numbering can differ from the same entry's RCSB mmCIF. Prefer mmCIF "
            "input when the target contains modified residues.",
            flush=True,
        )
        pdb_to_mmcif(src, dest)
    else:
        print(f"Passing mmCIF through unchanged: {src}", flush=True)
        shutil.copyfile(src, dest)


def build_inventory(cif_path):
    """Per-chain residue inventory with both numberings, read from _atom_site."""
    from Bio.PDB.MMCIF2Dict import MMCIF2Dict

    try:
        cif = MMCIF2Dict(cif_path)
    except Exception as exc:
        print(f"ERROR: Could not parse mmCIF '{cif_path}': {exc}", flush=True)
        sys.exit(1)

    required = ["_atom_site.label_asym_id", "_atom_site.label_comp_id"]
    for key in required:
        if key not in cif:
            print(f"ERROR: mmCIF is missing {key}; not a usable structure file.", flush=True)
            sys.exit(1)

    n_atoms = len(cif["_atom_site.label_asym_id"])

    def column(key, default="."):
        values = cif.get(key)
        if values is None:
            return [default] * n_atoms
        return values

    label_asym = column("_atom_site.label_asym_id")
    auth_asym = column("_atom_site.auth_asym_id")
    label_seq = column("_atom_site.label_seq_id")
    auth_seq = column("_atom_site.auth_seq_id")
    comp_id = column("_atom_site.label_comp_id")
    group = column("_atom_site.group_PDB", "ATOM")

    chains = {}
    seen = set()
    for i in range(n_atoms):
        if group[i] == "HETATM" and comp_id[i] not in THREE_TO_ONE:
            continue
        if label_seq[i] in (".", "?", ""):
            continue  # non-polymer: no label_seq_id, cannot appear in a contig
        key = (label_asym[i], label_seq[i])
        if key in seen:
            continue
        seen.add(key)
        chain = chains.setdefault(
            label_asym[i],
            {
                "label_asym_id": label_asym[i],
                "auth_asym_id": auth_asym[i],
                "residues": [],
            },
        )
        chain["residues"].append(
            {
                "label_seq_id": int(label_seq[i]),
                "auth_seq_id": auth_seq[i],
                "comp_id": comp_id[i],
                "one_letter": THREE_TO_ONE.get(comp_id[i], "X"),
            }
        )

    if not chains:
        print("ERROR: No polymer residues found in the target structure.", flush=True)
        sys.exit(1)

    inventory = {"source_file": os.path.basename(cif_path), "chains": []}
    for label_id in sorted(chains):
        chain = chains[label_id]
        residues = sorted(chain["residues"], key=lambda r: r["label_seq_id"])
        auth_to_label = {}
        label_to_auth = {}
        ambiguous = []
        for residue in residues:
            auth_key = str(residue["auth_seq_id"])
            if auth_key in auth_to_label:
                # Insertion codes make one author number cover several residues.
                ambiguous.append(auth_key)
                continue
            auth_to_label[auth_key] = residue["label_seq_id"]
            label_to_auth[str(residue["label_seq_id"])] = auth_key
        if ambiguous:
            print(
                f"WARNING: chain {chain['label_asym_id']} has repeated author numbers "
                f"({', '.join(sorted(set(ambiguous)))}), probably insertion codes. "
                "Author->label translation uses the first occurrence; give those "
                "residues in label numbering (residue_numbering=label) to be exact.",
                flush=True,
            )
        inventory["chains"].append(
            {
                "label_asym_id": chain["label_asym_id"],
                "auth_asym_id": chain["auth_asym_id"],
                "n_residues": len(residues),
                "label_seq_range": [residues[0]["label_seq_id"], residues[-1]["label_seq_id"]],
                "auth_seq_range": [residues[0]["auth_seq_id"], residues[-1]["auth_seq_id"]],
                "sequence": "".join(r["one_letter"] for r in residues),
                "residues": residues,
                "auth_to_label": auth_to_label,
                "label_to_auth": label_to_auth,
            }
        )
    return inventory


def main():
    os.makedirs("./outputs", exist_ok=True)

    pdb_id = os.environ.get("PARAM_PDB_ID", "")
    target_chain = os.environ.get("PARAM_TARGET_CHAIN", "A")

    upload = resolve_upload()
    target_cif = "./outputs/target.cif"

    if upload:
        if pdb_id:
            print(
                f"NOTE: Uploaded structure '{upload}' takes precedence over pdb_id='{pdb_id}'.",
                flush=True,
            )
        normalise(upload, target_cif)
        source = upload
    elif pdb_id:
        fetched = "./outputs/_fetched.cif"
        fetch_from_rcsb(pdb_id, fetched)
        normalise(fetched, target_cif)
        os.remove(fetched)
        source = f"RCSB:{pdb_id.upper()}"
    else:
        print(
            "ERROR: No target structure. Upload a .cif/.pdb/.cif.zst file or set pdb_id.",
            flush=True,
        )
        sys.exit(1)

    inventory = build_inventory(target_cif)
    inventory["source"] = source
    inventory["target_chain_param"] = target_chain

    label_ids = [c["label_asym_id"] for c in inventory["chains"]]
    auth_ids = [c["auth_asym_id"] for c in inventory["chains"]]
    matched = [
        c
        for c in inventory["chains"]
        if target_chain in (c["label_asym_id"], c["auth_asym_id"])
    ]
    if not matched:
        print(
            f"ERROR: target_chain '{target_chain}' is not in the structure. "
            f"label_asym_id present: {', '.join(label_ids)}; "
            f"auth_asym_id present: {', '.join(sorted(set(auth_ids)))}",
            flush=True,
        )
        sys.exit(1)
    inventory["target_chain_resolved"] = {
        "label_asym_id": matched[0]["label_asym_id"],
        "auth_asym_id": matched[0]["auth_asym_id"],
    }

    with open("./outputs/target_inventory.json", "w") as handle:
        json.dump(inventory, handle, indent=2)

    print(f"Target source: {source}", flush=True)
    print(f"Wrote {target_cif}", flush=True)
    for chain in inventory["chains"]:
        print(
            f"  chain label={chain['label_asym_id']} auth={chain['auth_asym_id']}: "
            f"{chain['n_residues']} residues, "
            f"label {chain['label_seq_range'][0]}-{chain['label_seq_range'][1]}, "
            f"auth {chain['auth_seq_range'][0]}-{chain['auth_seq_range'][1]}",
            flush=True,
        )
    print(
        f"Target chain '{target_chain}' resolved to "
        f"label_asym_id={inventory['target_chain_resolved']['label_asym_id']}",
        flush=True,
    )


if __name__ == "__main__":
    main()
