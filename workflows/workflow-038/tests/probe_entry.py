#!/usr/bin/env python3
"""Measure an mmCIF by hand, so the suite does not check the code with the code.

`_atom_site` is parsed directly here — no Biopython, no biotite, no bioemu2. The whole
point of the fixture suite is to compare what a node *reported* against what the file
*contains*, and a measurement taken through the same parser as the node would share the
node's bugs and agree with it for the wrong reason.

The quantities are chosen for the defects this workflow has actually had:

* `n_residues_by_chain` and `n_atoms_by_chain` — what node 01's inventory must match.
* `ca_atoms` versus `calcium_atoms` — the two are both named CA. mdtraj's `name CA`
  matches atom *names*, so a calcium ion counts as a C-alpha unless the selection says
  `protein`. Node 05 shifted every RMSD, every RMSF and the medoid on any system with
  structural calcium, and no output revealed it. These two numbers are what make that
  visible: the C-alpha count node 05 uses must equal `ca_atoms`, never `ca_atoms +
  calcium_atoms`.
* `unique_protein_sequences` — a heteromer (more than one) is due a positionally paired
  alignment from ColabFold; a monomer or homomer is not. Node 02 wrote only the unpaired
  half, and the model substitutes zeros for the paired feature it then never receives.
* `sg_atoms` — cysteine SG atoms, so a fixture with disulfides can be shown to exercise
  BioEmu2's `bad_bonds_sulfur` check rather than merely being described as doing so.

Run standalone to print the table:

    python3 -I tests/probe_entry.py <downloads-dir>
"""
import os
import sys

# The 20 standard amino acids plus the modified residues common enough to appear in a
# polymer chain of a fixture. Anything else is counted as non-standard, which is a
# measurement and not a judgement.
STANDARD_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
NUCLEIC = {
    "A", "C", "G", "U", "T", "DA", "DC", "DG", "DT", "DU", "I", "DI",
}
WATER = {"HOH", "DOD", "WAT"}

# Residue names that are a single metal ion. The ion's atom carries the same name, which
# is exactly the collision this probe exists to measure for calcium.
MONATOMIC_IONS = {
    "CA", "MG", "ZN", "FE", "FE2", "MN", "CU", "CU1", "NA", "K", "CL", "CD",
    "CO", "NI", "HG", "BA", "SR", "CS", "RB", "LI",
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def atom_site(path):
    """Yield one dict per `_atom_site` row, keyed by the column names the file declares.

    mmCIF loop columns are positional, so the header order is read rather than assumed.
    Rows are returned for the first model only when a model number column is present:
    an NMR ensemble repeats every atom per model and would otherwise multiply every
    count by the model count.
    """
    columns, rows, in_loop = [], [], False
    first_model = None
    with open(path, errors="replace") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("_atom_site."):
                columns.append(line.split(".", 1)[1].strip().split()[0])
                in_loop = True
                continue
            if not in_loop:
                continue
            if not line or line.startswith("#") or line.startswith("loop_"):
                if rows:
                    break
                continue
            if line.startswith("_"):
                break
            fields = _split_cif(line)
            if len(fields) != len(columns):
                continue
            row = dict(zip(columns, fields))
            model = row.get("pdbx_PDB_model_num")
            if model is not None:
                if first_model is None:
                    first_model = model
                elif model != first_model:
                    continue
            rows.append(row)
    return rows


def _split_cif(line):
    """Split one mmCIF data line, honouring single and double quoting."""
    out, token, quote = [], [], None
    for char in line:
        if quote:
            if char == quote:
                quote = None
            else:
                token.append(char)
            continue
        if char in "'\"":
            quote = char
            continue
        if char.isspace():
            if token:
                out.append("".join(token))
                token = []
            continue
        token.append(char)
    if token:
        out.append("".join(token))
    return out


def probe(path):
    """Measure one mmCIF. Every value here is counted, never inferred from the header."""
    rows = atom_site(path)
    if not rows:
        return {"error": "no _atom_site rows parsed", "path": path}

    def col(row, *names):
        for name in names:
            value = row.get(name)
            if value not in (None, "", ".", "?"):
                return value
        return None

    chains = {}
    ca_atoms = 0
    calcium_atoms = 0
    sg_atoms = 0
    hetero_rows = 0

    for row in rows:
        chain = col(row, "label_asym_id", "auth_asym_id") or "?"
        comp = (col(row, "label_comp_id", "auth_comp_id") or "?").upper()
        atom = (col(row, "label_atom_id", "auth_atom_id") or "?").upper()
        seq = col(row, "label_seq_id", "auth_seq_id") or "."
        group = (row.get("group_PDB") or "").upper()

        info = chains.setdefault(
            chain,
            {"residues": set(), "n_atoms": 0, "comps": set(), "sequence": []},
        )
        info["n_atoms"] += 1
        info["comps"].add(comp)
        key = (comp, seq)
        if key not in info["residues"]:
            info["residues"].add(key)
            if comp in STANDARD_AA:
                info["sequence"].append((seq, THREE_TO_ONE[comp]))
        if group == "HETATM":
            hetero_rows += 1

        # The collision this probe exists for: both are the atom name CA.
        if atom == "CA":
            if comp in MONATOMIC_IONS:
                calcium_atoms += 1 if comp == "CA" else 0
            elif comp in STANDARD_AA or comp not in NUCLEIC | WATER:
                ca_atoms += 1
        if atom == "SG" and comp == "CYS":
            sg_atoms += 1

    protein_chains, nucleic_chains, other_chains = [], [], []
    sequences = {}
    for chain, info in sorted(chains.items()):
        comps = info["comps"] - WATER
        standard = comps & STANDARD_AA
        nucleotides = comps & NUCLEIC
        if len(info["residues"]) > 1 and standard:
            protein_chains.append(chain)
            ordered = sorted(
                info["sequence"], key=lambda pair: _as_int(pair[0])
            )
            sequences[chain] = "".join(letter for _, letter in ordered)
        elif nucleotides and len(info["residues"]) > 1:
            nucleic_chains.append(chain)
        else:
            other_chains.append(chain)

    return {
        "path": path,
        "n_chains": len(chains),
        "chain_ids": sorted(chains),
        "protein_chains": protein_chains,
        "nucleic_chains": nucleic_chains,
        "other_chains": other_chains,
        "n_residues_by_chain": {c: len(i["residues"]) for c, i in sorted(chains.items())},
        "n_atoms_by_chain": {c: i["n_atoms"] for c, i in sorted(chains.items())},
        "n_polymer_residues": sum(
            len(chains[c]["residues"]) for c in protein_chains + nucleic_chains
        ),
        # Every component on a chain that is neither protein nor nucleic acid, measured from
        # the file. This is what `include_ligands` and `exclude_ligands` act on, so a check
        # on those parameters can compare the inventory against the entry rather than against
        # another run of the same code.
        "nonpolymer_components": sorted(
            {
                comp
                for chain in other_chains
                for comp in chains[chain]["comps"]
                if comp not in WATER
            }
        ),
        "ca_atoms": ca_atoms,
        "calcium_atoms": calcium_atoms,
        "sg_atoms": sg_atoms,
        "hetero_rows": hetero_rows,
        "sequences": sequences,
        "unique_protein_sequences": len({s for s in sequences.values() if s}),
        "has_nucleic_acid": bool(nucleic_chains),
    }


def _as_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 2
    directory = sys.argv[1]
    names = sorted(
        name for name in os.listdir(directory) if name.endswith((".cif", ".pdb"))
    )
    if not names:
        print(f"no structures in {directory}")
        return 1
    header = f"{'entry':<12} {'chains':>6} {'prot':>5} {'seqs':>5} {'res':>6} {'CA':>5} {'Ca2+':>5} {'SG':>4}"
    print(header)
    print("-" * len(header))
    for name in names:
        measured = probe(os.path.join(directory, name))
        if "error" in measured:
            print(f"{name:<12} {measured['error']}")
            continue
        print(
            f"{name:<12} {measured['n_chains']:>6} {len(measured['protein_chains']):>5} "
            f"{measured['unique_protein_sequences']:>5} "
            f"{measured['n_polymer_residues']:>6} {measured['ca_atoms']:>5} "
            f"{measured['calcium_atoms']:>5} {measured['sg_atoms']:>4}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
