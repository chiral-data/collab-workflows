#!/usr/bin/env python3
"""Node 01 — Annotate Sequences.

Turn a panel of raw VH/VL (or VHH) sequences into a numbered, region-delimited,
liability-scanned annotation that every downstream node keys off.

Input: one CSV (clone_id,vh,vl) or FASTA (ids suffixed _VH/_VL) in ./inputs/.

Why this node exists separately from node 02, which also numbers: node 02 numbers *inside*
ImmuneBuilder, for modelling. This node owns the annotation — chain type, germline, CDR
boundaries, liabilities — and it runs in milliseconds per clone, so a malformed panel fails
here, on CPU, before any model is loaded.

Two numbering passes may happen. CDR definitions are each expressed in the numbering scheme
they were defined on (North on Chothia numbering, Kabat on Kabat numbering, and so on), so when
`cdr_definition` does not match `numbering_scheme` the sequence is numbered twice: once in the
reporting scheme and once in the definition's own scheme. Both are recorded per residue, keyed
on the residue's index in the input sequence, which is what lets node 03 map structure residues
back to regions no matter which scheme the PDB carries.

Per-clone failures are recorded, never fatal: one unnumberable sequence in a 500-clone panel
yields a row with a reason and the panel continues. The node fails hard only if the panel
cannot be read at all, or if *every* clone failed.

TAP uses IMGT numbering with IMGT CDR definitions, which is why `cdr_definition` defaults to
imgt rather than north — see the pipeline doc, "Corrections from the design recheck".
"""
import csv
import glob
import json
import os
import re
import sys

VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

CSV_EXT = (".csv",)
FASTA_EXT = (".fasta", ".fa", ".faa", ".fas")

HEAVY_SUFFIXES = ("_VH", "_H", "_HEAVY", "_VHH")
LIGHT_SUFFIXES = ("_VL", "_L", "_LIGHT", "_K", "_KAPPA", "_LAMBDA")

# CDR region start coordinates, from prihoda/abnumber (MIT) abnumber/common.py:
#   [CDR1, FR2, CDR2, FR3, CDR3, FR4, end]
# Each definition is expressed in its own numbering scheme; CDR_DEFINITION_SCHEME maps it.
SCHEME_BORDERS = {
    "imgt": [27, 39, 56, 66, 105, 118, 129],
    "kabat_H": [31, 36, 50, 66, 95, 103, 114],
    "kabat_K": [24, 35, 50, 57, 89, 98, 108],
    "kabat_L": [24, 35, 50, 57, 89, 98, 108],
    "chothia_H": [26, 33, 52, 57, 95, 103, 114],
    "chothia_K": [24, 35, 50, 57, 89, 98, 108],
    "chothia_L": [24, 35, 50, 57, 89, 98, 108],
    "north_H": [23, 36, 50, 59, 93, 103, 114],
    "north_K": [24, 35, 49, 57, 89, 98, 108],
    "north_L": [24, 35, 49, 57, 89, 98, 108],
}

# North is defined on Chothia numbering; the rest on their own.
CDR_DEFINITION_SCHEME = {
    "imgt": "imgt",
    "kabat": "kabat",
    "chothia": "chothia",
    "north": "chothia",
}

NUMBERING_SCHEMES = {"imgt", "chothia", "kabat", "martin", "aho"}

# Liability motifs. Each entry is (name, regex, severity, note). These are pattern matches on
# sequence and nothing more: they over-call by design (most N-glycosylation sequons are never
# glycosylated), and node 03's solvent exposure is what separates a real risk from a hit.
DEFAULT_MOTIFS = [
    {
        "name": "n_glycosylation",
        "pattern": "N[^P][ST]",
        "severity": "high",
        "note": "N-X-S/T sequon (X not proline); a candidate N-linked glycosylation site",
    },
    {
        "name": "deamidation_ng",
        "pattern": "NG",
        "severity": "high",
        "note": "Asn-Gly, the fastest deamidation motif",
    },
    {
        "name": "deamidation_ns",
        "pattern": "N[SHAT]",
        "severity": "medium",
        "note": "Asn followed by Ser/His/Ala/Thr; slower deamidation",
    },
    {
        "name": "isomerization_dg",
        "pattern": "DG",
        "severity": "high",
        "note": "Asp-Gly, isomerization/succinimide formation",
    },
    {
        "name": "isomerization_ds",
        "pattern": "D[SATH]",
        "severity": "medium",
        "note": "Asp followed by Ser/Ala/Thr/His; slower isomerization",
    },
    {
        "name": "acid_cleavage_dp",
        "pattern": "D[PG]",
        "severity": "medium",
        "note": "Asp-Pro/Asp-Gly, acid-labile peptide bond",
    },
    {
        "name": "oxidation_met",
        "pattern": "M",
        "severity": "low",
        "note": "Methionine; oxidation-prone when solvent-exposed",
    },
    {
        "name": "oxidation_trp",
        "pattern": "W",
        "severity": "low",
        "note": "Tryptophan; oxidation-prone when solvent-exposed",
    },
    {
        "name": "fragmentation_cleavage",
        "pattern": "[ST]C",
        "severity": "low",
        "note": "Ser/Thr-Cys; reported fragmentation motif",
    },
]

# Cysteine parity is a count rule, not a regex, so it is handled separately.
CYS_RULE = {
    "name": "unpaired_cysteine",
    "severity": "high",
    "note": "Odd number of cysteines in the domain; an unpaired thiol risks scrambling",
}


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def warn(message):
    print(f"WARNING: {message}", flush=True)


def choice(name, value, allowed):
    if value not in allowed:
        fail(f"{name} must be one of {sorted(allowed)}, got '{value}'")
    return value


def find_panel():
    for ext in CSV_EXT + FASTA_EXT:
        matches = sorted(glob.glob(f"./inputs/*{ext}"))
        if matches:
            return matches[0]
    return None


def clean_sequence(raw):
    return re.sub(r"[\s\-\*]", "", (raw or "")).upper()


def read_csv_panel(path):
    """clone_id + vh/vl columns; unknown columns are carried through untouched."""
    with open(path, newline="") as handle:
        sniffed = csv.DictReader(handle)
        if not sniffed.fieldnames:
            fail(f"'{path}' has no header row.")
        lowered = {name.lower().strip(): name for name in sniffed.fieldnames}

        id_col = next((lowered[k] for k in ("clone_id", "id", "name") if k in lowered), None)
        vh_col = next((lowered[k] for k in ("vh", "vhh", "heavy", "h") if k in lowered), None)
        vl_col = next((lowered[k] for k in ("vl", "light", "l") if k in lowered), None)
        if not vh_col:
            fail(
                f"'{path}' has no heavy-chain column. Expected one of vh/vhh/heavy/h; "
                f"found {sniffed.fieldnames}. See input_files/README.md."
            )

        clones = []
        for row_number, row in enumerate(sniffed, start=2):
            extras = {
                key: value
                for key, value in row.items()
                if key not in {id_col, vh_col, vl_col} and key and value not in (None, "")
            }
            clone_id = (row.get(id_col) or "").strip() if id_col else ""
            clones.append(
                {
                    "clone_id": clone_id or f"clone_{row_number - 1:04d}",
                    "vh": clean_sequence(row.get(vh_col)),
                    "vl": clean_sequence(row.get(vl_col)) if vl_col else "",
                    "extras": extras,
                }
            )
    return clones


def read_fasta_panel(path):
    """One record per chain; the chain is taken from the record id's suffix."""
    records = []
    name, chunks = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name, chunks = line[1:].split()[0], []
            elif line:
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    if not records:
        fail(f"No FASTA records found in '{path}'.")

    clones, order = {}, []
    for record_id, sequence in records:
        upper = record_id.upper()
        chain = None
        for suffix in HEAVY_SUFFIXES:
            if upper.endswith(suffix):
                chain, base = "vh", record_id[: -len(suffix)]
                break
        if chain is None:
            for suffix in LIGHT_SUFFIXES:
                if upper.endswith(suffix):
                    chain, base = "vl", record_id[: -len(suffix)]
                    break
        if chain is None:
            fail(
                f"FASTA record '{record_id}' has no recognised chain suffix. Expected one of "
                f"{list(HEAVY_SUFFIXES + LIGHT_SUFFIXES)} (case-insensitive). This is an error "
                "rather than a guess: mis-assigning a chain would silently mis-model the clone."
            )
        base = base.rstrip("_-") or record_id
        if base not in clones:
            clones[base] = {"clone_id": base, "vh": "", "vl": "", "extras": {}}
            order.append(base)
        if clones[base][chain]:
            fail(f"FASTA has two {chain.upper()} records for clone '{base}'.")
        clones[base][chain] = clean_sequence(sequence)
    return [clones[key] for key in order]


def borders_for(definition, chain_type):
    """Border list for a CDR definition and ANARCI chain type (H, K or L)."""
    if definition == "imgt":
        return SCHEME_BORDERS["imgt"]
    key = f"{definition}_{'H' if chain_type == 'H' else chain_type}"
    if key not in SCHEME_BORDERS:
        key = f"{definition}_{'H' if chain_type == 'H' else 'K'}"
    return SCHEME_BORDERS[key]


def regions_from_borders(borders):
    """{position -> region} plus the anchor positions, derived rather than hardcoded.

    Anchors are the positions immediately flanking each CDR, which is what TAP's CDR-vicinity
    definition ("every surface-exposed IMGT-defined CDR and anchor residue") refers to.
    """
    cdr1, fr2, cdr2, fr3, cdr3, fr4, end = borders
    spans = [
        ("FR1", 1, cdr1 - 1),
        ("CDR1", cdr1, fr2 - 1),
        ("FR2", fr2, cdr2 - 1),
        ("CDR2", cdr2, fr3 - 1),
        ("FR3", fr3, cdr3 - 1),
        ("CDR3", cdr3, fr4 - 1),
        ("FR4", fr4, end - 1),
    ]
    position_to_region = {}
    for region, start, stop in spans:
        for position in range(start, stop + 1):
            position_to_region[position] = region

    anchors = set()
    for _, start, stop in [span for span in spans if span[0].startswith("CDR")]:
        anchors.update({start - 1, stop + 1})
    return position_to_region, sorted(anchors)


def number_chain(sequence, scheme, allowed_species, allow_chains):
    """ANARCI one chain. Returns (numbering, details) or raises ValueError with a reason."""
    from anarci import anarci

    try:
        numbered, details, _ = anarci(
            [("query", sequence)],
            scheme=scheme,
            output=False,
            allow=set(allow_chains),
            assign_germline=True,
            allowed_species=allowed_species,
        )
    except Exception as exc:  # ANARCI raises on malformed input and on a missing HMM library
        raise ValueError(f"ANARCI failed: {exc}") from exc

    if not numbered or not numbered[0]:
        raise ValueError(
            f"ANARCI did not recognise this as an antibody variable domain of type "
            f"{sorted(allow_chains)}. Fabs, full IgGs, scFvs and CDR fragments are rejected "
            "here by design; trim to the variable domain first."
        )
    if len(numbered[0]) > 1:
        raise ValueError(
            f"ANARCI found {len(numbered[0])} variable domains in one sequence. This looks like "
            "an scFv or a multi-domain construct; split it at the linker first (node 01 does "
            "not cut linkers, because a silent mis-split is worse than a rejection)."
        )
    return numbered[0][0], details[0][0]


def residue_table(domain, sequence):
    """Map numbered residues onto their index in the input sequence.

    ANARCI returns the numbering for the domain it found between start_index and end_index,
    with '-' for scheme positions the sequence does not fill. Walking the numbering in order
    and advancing through the input sequence gives each numbered residue its original index,
    which is the key that lets two different numbering schemes be joined.
    """
    numbering, start_index, _end_index = domain
    rows, cursor = [], start_index
    for (position, insertion), aa in numbering:
        if aa == "-":
            continue
        if cursor >= len(sequence) or sequence[cursor] != aa:
            # Resynchronise rather than trust the cursor: ANARCI should be in lockstep with
            # the sequence here, and a mismatch means the assumption is broken.
            found = sequence.find(aa, cursor)
            if found < 0:
                raise ValueError(
                    "Could not align ANARCI numbering back onto the input sequence "
                    f"at position {position}{insertion.strip()}"
                )
            cursor = found
        rows.append(
            {
                "seq_index": cursor,
                "aa": aa,
                "position": int(position),
                "insertion": insertion.strip(),
            }
        )
        cursor += 1
    return rows


def scan_liabilities(clone_id, chain_label, rows, motifs, position_to_region):
    """Motif hits over the numbered domain, reported with region and numbered position."""
    domain_sequence = "".join(row["aa"] for row in rows)
    hits = []

    for motif in motifs:
        for match in re.finditer(f"(?={motif['pattern']})", domain_sequence):
            start = match.start()
            width = len(re.match(motif["pattern"], domain_sequence[start:]).group(0))
            covered = rows[start : start + width]
            if not covered:
                continue
            regions = sorted(
                {position_to_region.get(row["position"], "FR?") for row in covered}
            )
            hits.append(
                {
                    "clone_id": clone_id,
                    "chain": chain_label,
                    "motif": motif["name"],
                    "severity": motif["severity"],
                    "residues": domain_sequence[start : start + width],
                    "seq_index": covered[0]["seq_index"],
                    "position": f"{covered[0]['position']}{covered[0]['insertion']}",
                    "region": ",".join(regions),
                    "in_cdr": any(region.startswith("CDR") for region in regions),
                    "note": motif["note"],
                }
            )

    cysteines = [row for row in rows if row["aa"] == "C"]
    if len(cysteines) % 2 == 1:
        hits.append(
            {
                "clone_id": clone_id,
                "chain": chain_label,
                "motif": CYS_RULE["name"],
                "severity": CYS_RULE["severity"],
                "residues": "C",
                "seq_index": cysteines[-1]["seq_index"],
                "position": f"{cysteines[-1]['position']}{cysteines[-1]['insertion']}",
                "region": position_to_region.get(cysteines[-1]["position"], "FR?"),
                "in_cdr": str(
                    position_to_region.get(cysteines[-1]["position"], "")
                ).startswith("CDR"),
                "note": f"{len(cysteines)} cysteines in the domain. {CYS_RULE['note']}",
            }
        )
    return hits


def load_motifs(spec):
    if spec in ("", "default"):
        return DEFAULT_MOTIFS, "default"
    if not os.path.exists(spec):
        fail(f"liability_motifs='{spec}' is neither 'default' nor an existing file.")
    try:
        with open(spec) as handle:
            loaded = json.load(handle)
    except Exception as exc:
        fail(f"Could not parse the custom motif file '{spec}': {exc}")
    if not isinstance(loaded, list) or not loaded:
        fail(f"'{spec}' must contain a non-empty JSON list of motif objects.")
    for entry in loaded:
        missing = {"name", "pattern"} - set(entry or {})
        if missing:
            fail(f"Motif entry {entry} is missing {sorted(missing)}.")
        try:
            re.compile(entry["pattern"])
        except re.error as exc:
            fail(f"Motif '{entry['name']}' has an invalid regex: {exc}")
        entry.setdefault("severity", "medium")
        entry.setdefault("note", "custom motif")
    return loaded, os.path.basename(spec)


def annotate_chain(clone_id, chain_label, sequence, settings):
    """Number, region-delimit and liability-scan one chain. Raises ValueError on rejection."""
    bad = sorted(set(sequence) - VALID_AA)
    if bad:
        raise ValueError(
            f"{chain_label} contains characters that are not standard amino acids: {bad}. "
            "Ambiguity codes are rejected rather than modelled."
        )
    if len(sequence) <= 70:
        raise ValueError(
            f"{chain_label} is {len(sequence)} residues. ImmuneBuilder asserts a variable "
            "domain is longer than 70, so this would fail in node 02; give the whole domain."
        )

    allow = {"H"} if chain_label == "VH" else {"K", "L"}
    domain, details = number_chain(
        sequence, settings["numbering_scheme"], settings["allowed_species"], allow
    )
    rows = residue_table(domain, sequence)
    chain_type = details.get("chain_type", "H" if chain_label == "VH" else "K")

    # CDR definition may live in another numbering scheme; number again if so and join on
    # seq_index, which is scheme-independent.
    definition_scheme = CDR_DEFINITION_SCHEME[settings["cdr_definition"]]
    if definition_scheme == settings["numbering_scheme"]:
        definition_rows = rows
    else:
        definition_domain, _ = number_chain(
            sequence, definition_scheme, settings["allowed_species"], allow
        )
        definition_rows = residue_table(definition_domain, sequence)

    borders = borders_for(settings["cdr_definition"], chain_type)
    position_to_region, anchors = regions_from_borders(borders)
    region_by_index = {
        row["seq_index"]: position_to_region.get(row["position"], "FR?")
        for row in definition_rows
    }
    anchor_indices = {
        row["seq_index"] for row in definition_rows if row["position"] in set(anchors)
    }

    for row in rows:
        row["region"] = region_by_index.get(row["seq_index"], "outside")
        row["is_anchor"] = row["seq_index"] in anchor_indices

    cdr_lengths = {}
    for cdr in ("CDR1", "CDR2", "CDR3"):
        cdr_lengths[cdr] = sum(1 for row in rows if row["region"] == cdr)

    germlines = details.get("germlines") or {}

    def germline_entry(key):
        entry = germlines.get(key)
        if not entry:
            return None, None
        try:
            (species, gene), identity = entry
            return f"{species}:{gene}" if gene else None, round(float(identity), 4)
        except Exception:
            return None, None

    v_gene, v_identity = germline_entry("v_gene")
    j_gene, j_identity = germline_entry("j_gene")

    liabilities = scan_liabilities(
        clone_id, chain_label, rows, settings["motifs"], position_to_region
    )

    return {
        "chain": chain_label,
        "chain_type": chain_type,
        "species": details.get("species"),
        "bitscore": details.get("bitscore"),
        "numbering_scheme": settings["numbering_scheme"],
        "cdr_definition": settings["cdr_definition"],
        "cdr_definition_scheme": definition_scheme,
        "input_length": len(sequence),
        "domain_length": len(rows),
        "domain_start_index": domain[1],
        "domain_end_index": domain[2],
        "domain_sequence": "".join(row["aa"] for row in rows),
        "v_gene": v_gene,
        "v_identity": v_identity,
        "j_gene": j_gene,
        "j_identity": j_identity,
        "cdr_lengths": cdr_lengths,
        "cdrs": {
            cdr: "".join(row["aa"] for row in rows if row["region"] == cdr)
            for cdr in ("CDR1", "CDR2", "CDR3")
        },
        "n_cysteines": sum(1 for row in rows if row["aa"] == "C"),
        "residues": rows,
        "anchor_positions": anchors,
        "liabilities": liabilities,
    }


def main():
    os.makedirs("./outputs", exist_ok=True)

    receptor_type = choice(
        "receptor_type",
        os.environ.get("PARAM_RECEPTOR_TYPE", "auto"),
        {"auto", "antibody", "nanobody"},
    )
    numbering_scheme = choice(
        "numbering_scheme",
        os.environ.get("PARAM_NUMBERING_SCHEME", "imgt"),
        NUMBERING_SCHEMES,
    )
    cdr_definition = choice(
        "cdr_definition",
        os.environ.get("PARAM_CDR_DEFINITION", "imgt"),
        set(CDR_DEFINITION_SCHEME),
    )
    allowed_species = [
        name.strip().lower()
        for name in os.environ.get("PARAM_ALLOWED_SPECIES", "human,mouse").split(",")
        if name.strip()
    ] or None
    motifs, motif_source = load_motifs(os.environ.get("PARAM_LIABILITY_MOTIFS", "default"))

    if cdr_definition != "imgt":
        print(
            f"NOTE: cdr_definition={cdr_definition}. TAP itself uses IMGT CDR definitions, so "
            "the Table 2 thresholds in node 03 were calibrated against IMGT boundaries. CDR "
            "length in particular is not comparable across definitions.",
            flush=True,
        )

    panel_path = find_panel()
    if not panel_path:
        fail(
            "No panel found in ./inputs/. Supply a CSV (clone_id,vh,vl) or a FASTA with "
            "_VH/_VL record suffixes. See input_files/README.md."
        )
    print(f"Reading panel: {panel_path}", flush=True)
    clones = (
        read_csv_panel(panel_path)
        if panel_path.lower().endswith(CSV_EXT)
        else read_fasta_panel(panel_path)
    )
    if not clones:
        fail(f"'{panel_path}' contains no clones.")

    seen = set()
    for clone in clones:
        if clone["clone_id"] in seen:
            fail(f"Duplicate clone_id '{clone['clone_id']}'. Clone ids must be unique.")
        seen.add(clone["clone_id"])

    settings = {
        "receptor_type": receptor_type,
        "numbering_scheme": numbering_scheme,
        "cdr_definition": cdr_definition,
        "cdr_definition_scheme": CDR_DEFINITION_SCHEME[cdr_definition],
        "allowed_species": allowed_species,
        "liability_motifs": motif_source,
        "panel_file": os.path.basename(panel_path),
    }
    print(f"Resolved settings: {json.dumps(settings, indent=2)}", flush=True)

    annotations, liability_rows, failures = [], [], []
    for clone in clones:
        clone_id = clone["clone_id"]
        vh, vl = clone["vh"], clone["vl"]

        if receptor_type == "antibody":
            kind = "antibody"
        elif receptor_type == "nanobody":
            kind = "nanobody"
        else:
            kind = "antibody" if (vh and vl) else "nanobody"

        if not vh:
            failures.append(
                {"clone_id": clone_id, "stage": "input", "reason": "no heavy chain supplied"}
            )
            continue
        if kind == "antibody" and not vl:
            failures.append(
                {
                    "clone_id": clone_id,
                    "stage": "input",
                    "reason": (
                        "receptor_type=antibody but no light chain supplied. Set "
                        "receptor_type=nanobody for VHH panels."
                    ),
                }
            )
            continue
        if kind == "nanobody" and vl:
            warn(
                f"{clone_id}: receptor_type=nanobody, so the light chain supplied "
                "is ignored."
            )
            vl = ""

        record = {
            "clone_id": clone_id,
            "receptor_type": kind,
            "extras": clone["extras"],
            "chains": {},
        }
        try:
            record["chains"]["VH"] = annotate_chain(clone_id, "VH", vh, {**settings, "motifs": motifs})
            if kind == "antibody":
                record["chains"]["VL"] = annotate_chain(
                    clone_id, "VL", vl, {**settings, "motifs": motifs}
                )
        except ValueError as exc:
            failures.append({"clone_id": clone_id, "stage": "annotate", "reason": str(exc)})
            print(f"  {clone_id}: REJECTED — {exc}", flush=True)
            continue

        total_cdr_length = sum(
            sum(chain["cdr_lengths"].values()) for chain in record["chains"].values()
        )
        record["total_cdr_length"] = total_cdr_length
        for chain in record["chains"].values():
            liability_rows.extend(chain["liabilities"])
        annotations.append(record)

        summary = "/".join(
            f"{name}:{chain['chain_type']}" for name, chain in record["chains"].items()
        )
        print(
            f"  {clone_id}: {kind}, {summary}, total CDR length {total_cdr_length}, "
            f"{sum(len(chain['liabilities']) for chain in record['chains'].values())} "
            "liability hit(s)",
            flush=True,
        )

    if not annotations:
        fail(
            f"Every clone in the panel failed annotation ({len(failures)} of {len(clones)}). "
            "Check that the sequences are whole variable domains; see the reasons above."
        )

    payload = {
        "settings": settings,
        "motifs": motifs,
        "cys_rule": CYS_RULE,
        "clones": annotations,
        "failures": failures,
    }
    with open("./outputs/annotations.json", "w") as handle:
        json.dump(payload, handle, indent=2)

    summary_fields = [
        "clone_id",
        "receptor_type",
        "vh_chain_type",
        "vl_chain_type",
        "vh_length",
        "vl_length",
        "total_cdr_length",
        "cdrh1",
        "cdrh2",
        "cdrh3",
        "cdrl1",
        "cdrl2",
        "cdrl3",
        "cdrh3_length",
        "vh_v_gene",
        "vh_v_identity",
        "vl_v_gene",
        "vl_v_identity",
        "n_liabilities",
        "n_cdr_liabilities",
        "status",
        "reason",
    ]
    with open("./outputs/annotations.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, extrasaction="ignore")
        writer.writeheader()
        for record in annotations:
            heavy = record["chains"]["VH"]
            light = record["chains"].get("VL")
            hits = [hit for chain in record["chains"].values() for hit in chain["liabilities"]]
            writer.writerow(
                {
                    "clone_id": record["clone_id"],
                    "receptor_type": record["receptor_type"],
                    "vh_chain_type": heavy["chain_type"],
                    "vl_chain_type": light["chain_type"] if light else "",
                    "vh_length": heavy["domain_length"],
                    "vl_length": light["domain_length"] if light else 0,
                    "total_cdr_length": record["total_cdr_length"],
                    "cdrh1": heavy["cdrs"]["CDR1"],
                    "cdrh2": heavy["cdrs"]["CDR2"],
                    "cdrh3": heavy["cdrs"]["CDR3"],
                    "cdrl1": light["cdrs"]["CDR1"] if light else "",
                    "cdrl2": light["cdrs"]["CDR2"] if light else "",
                    "cdrl3": light["cdrs"]["CDR3"] if light else "",
                    "cdrh3_length": heavy["cdr_lengths"]["CDR3"],
                    "vh_v_gene": heavy["v_gene"] or "",
                    "vh_v_identity": heavy["v_identity"] if heavy["v_identity"] is not None else "",
                    "vl_v_gene": (light["v_gene"] or "") if light else "",
                    "vl_v_identity": (
                        light["v_identity"] if light and light["v_identity"] is not None else ""
                    ),
                    "n_liabilities": len(hits),
                    "n_cdr_liabilities": sum(1 for hit in hits if hit["in_cdr"]),
                    "status": "ok",
                    "reason": "",
                }
            )
        for failure in failures:
            writer.writerow(
                {
                    "clone_id": failure["clone_id"],
                    "status": "failed",
                    "reason": failure["reason"],
                }
            )

    liability_fields = [
        "clone_id",
        "chain",
        "motif",
        "severity",
        "residues",
        "position",
        "region",
        "in_cdr",
        "seq_index",
        "note",
    ]
    with open("./outputs/liabilities.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=liability_fields, extrasaction="ignore")
        writer.writeheader()
        for hit in liability_rows:
            writer.writerow(hit)

    with open("./outputs/clones.fasta", "w") as handle:
        for record in annotations:
            for name, chain in record["chains"].items():
                handle.write(f">{record['clone_id']}_{name}\n{chain['domain_sequence']}\n")

    manifest = {
        "settings": settings,
        "n_clones_in": len(clones),
        "n_annotated": len(annotations),
        "n_failed": len(failures),
        "failures": failures,
        "n_liability_hits": len(liability_rows),
        "receptor_types": {
            kind: sum(1 for record in annotations if record["receptor_type"] == kind)
            for kind in ("antibody", "nanobody")
        },
        "anarci_version": anarci_version(),
        "caveat": (
            "Liability motifs are sequence pattern matches: they over-call, and node 03's "
            "solvent exposure is what distinguishes a real risk from a hit. CDR lengths are a "
            "function of cdr_definition and are not comparable across definitions."
        ),
    }
    with open("./outputs/panel_manifest.json", "w") as handle:
        json.dump(manifest, handle, indent=2)

    print(
        f"Annotated {len(annotations)}/{len(clones)} clone(s); {len(liability_rows)} liability "
        f"hit(s); {len(failures)} failure(s)",
        flush=True,
    )
    if failures:
        warn(
            f"{len(failures)} clone(s) failed and carry a reason in annotations.csv. They are "
            "excluded from modelling but still present in the panel record."
        )
    nanobodies = manifest["receptor_types"]["nanobody"]
    if nanobodies:
        warn(
            f"{nanobodies} nanobody clone(s). TAP's thresholds were derived from paired VH/VL "
            "therapeutics; node 03 marks VHH rows as not calibrated (open item 9)."
        )


def anarci_version():
    try:
        from importlib.metadata import version

        return version("ANARCI")
    except Exception:
        return "unknown"


if __name__ == "__main__":
    main()
