#!/usr/bin/env python3
"""Node 03 — Developability Metrics (TAP-style).

The five developability metrics of Raybould et al., PNAS 2019 (DOI 10.1073/pnas.1810576116,
CC BY), reimplemented from the paper's Methods because TAP itself ships as a web service with
no redistributable source. Everything this node emits is labelled TAP-*style*: see the pipeline
doc, "TAP is a paper, not a package".

Definitions implemented, as the paper states them:

  Surface-exposed   >= 7.5% relative exposure, Shrake-Rupley algorithm.
  CDR vicinity      every surface-exposed IMGT-defined CDR and anchor residue, plus all other
                    surface-exposed residues with a heavy atom within a 4 A radius.
  Salt bridges      K/R vs D/E pairs with an N+-O- distance <= 3.2 A.
  PSH               sum over surface-exposed residue pairs with closest heavy-atom distance
                    r12 < 7.5 A of H(R1,S)*H(R2,S) / r12^2, H normalized to [1, 2].
                    Salt-bridge residues take glycine's hydrophobicity value.
  PPC / PNC         analogous in form to PSH with H replaced by |Q(R)|, over positively and
                    negatively charged residues respectively.
  SFvCSP            (sum of Q over surface-exposed VH residues) * (same for VL).
  Charges           D -1, E -1, K +1, R +1, H +0.1. Salt-bridge residues take 0; buried side
                    chains are ignored entirely.

Two deliberate deviations, both recorded in developability.json under `implementation` and
`calibration`:

  1. TAP measures *side-chain* relative exposure against an Ala-R-Ala reference from Yang Zhu &
     Blundell (1996), which is not available to us. We use whole-residue relative ASA against
     Tien et al. (2013) theoretical Ala-X-Ala values. Including backbone area in numerator and
     denominator shifts which residues count as exposed, and every metric is built on that set.
  2. TAP's thresholds were calibrated on ABodyBuilder v1 models; ours come from ABodyBuilder2.

The pair geometry is deliberately plain Python rather than vectorised: the sums run over the
CDR vicinity (tens of residues), a centroid pre-filter removes most pairs, and keeping it
dependency-free means the metric maths is unit-testable without numpy or biotite.

The structure arrives from node 02 numbered in the user's scheme, so residues are joined to
node 01's regions on (chain, residue id, insertion code).
"""
import csv
import json
import math
import os
import sys

# Tien et al. 2013 (PLoS ONE, PMC3836772) Table 1, "Theoretical" column: Ala-X-Ala
# whole-residue maximum ASA in square Angstrom. Used as the exposure denominator.
TIEN_THEORETICAL_MAX_ASA = {
    "A": 129.0, "R": 274.0, "N": 195.0, "D": 193.0, "C": 167.0,
    "E": 223.0, "Q": 225.0, "G": 104.0, "H": 224.0, "I": 197.0,
    "L": 201.0, "K": 236.0, "M": 224.0, "F": 240.0, "P": 159.0,
    "S": 155.0, "T": 172.0, "W": 285.0, "Y": 263.0, "V": 174.0,
}

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # Protonation variants OpenMM/pdbfixer may write after refinement.
    "HID": "H", "HIE": "H", "HIP": "H", "CYX": "C", "CYM": "C",
    "LYN": "K", "ASH": "D", "GLH": "E",
}

# Paper: "aspartic acid, -1; glutamic acid, -1; lysine, +1; arginine, +1; and histidine, +0.1
# (Henderson-Hasselbalch equation applied: pKa 6, pH 7.4, and rounded up to one decimal place)."
RESIDUE_CHARGE = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1}

# Salt-bridge donor/acceptor atoms: "pairs of lysines/arginines and aspartic acids/glutamic
# acids with a N+-O- distance <= 3.2 A".
CATION_ATOMS = {"K": {"NZ"}, "R": {"NE", "NH1", "NH2"}}
ANION_ATOMS = {"D": {"OD1", "OD2"}, "E": {"OE1", "OE2"}}
SALT_BRIDGE_CUTOFF = 3.2

EXPOSURE_THRESHOLD = 0.075
VICINITY_RADIUS = 4.0
PAIR_CUTOFF = 7.5

# Raybould et al. 2019, Table 2: "TAP amber and red flag regions, as defined by the entire set
# of 242 CSTs". Amber bands are inclusive; red is strictly beyond them.
TAP_THRESHOLDS = {
    "total_cdr_length": {
        "amber": [(54.0, 60.0)],
        "red_above": 60.0,
        "red_below": None,
    },
    "psh_cdr_vicinity": {
        "amber": [(83.84, 100.71), (156.20, 173.85)],
        "red_above": 173.85,
        "red_below": 83.84,
    },
    "ppc_cdr_vicinity": {
        "amber": [(1.25, 3.16)],
        "red_above": 3.16,
        "red_below": None,
    },
    "pnc_cdr_vicinity": {
        "amber": [(1.84, 3.50)],
        "red_above": 3.50,
        "red_below": None,
    },
    "sfvcsp": {
        "amber": [(-20.40, -6.30)],
        "red_above": None,
        "red_below": -20.40,
    },
}

# Published reference values, for detecting a gross calibration offset rather than for flagging.
CST_REFERENCE = {
    "psh_cdr_vicinity_mean": 123.30,
    "psh_cdr_vicinity_sd": 16.60,
    "n_models": 137,
    "source": "Raybould et al. 2019, Fig. 2A (137 CST models)",
}

CHAIN_FOR = {"VH": "H", "VL": "L"}


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def warn(message):
    print(f"WARNING: {message}", flush=True)


def read_json(path, label):
    if not os.path.exists(path):
        fail(f"{label} not found: {path}. The upstream node must run first.")
    try:
        with open(path) as handle:
            return json.load(handle)
    except Exception as exc:
        fail(f"Could not read {label} ('{path}'): {exc}")


def normalized_hydrophobicity():
    """Kyte & Doolittle, normalized to [1, 2] as the paper specifies.

    The scale itself comes from biopython (Bio.SeqUtils.ProtParamData.kd) rather than being
    transcribed here, so the numbers are the library's.
    """
    from Bio.SeqUtils.ProtParamData import kd

    low, high = min(kd.values()), max(kd.values())
    span = high - low
    return {aa: 1.0 + (value - low) / span for aa, value in kd.items()}


def load_structure(path):
    """Heavy atoms of model 1, as a list of per-residue dicts.

    SASA is computed on the heavy-atom structure: Shrake-Rupley with ProtOr radii is defined
    for heavy atoms, and refined models carry hydrogens that would otherwise occlude.
    """
    import biotite.structure as struc
    from biotite.structure.io.pdb import PDBFile

    atoms = PDBFile.read(path).get_structure(model=1)
    atoms = atoms[struc.filter_amino_acids(atoms)]
    atoms = atoms[atoms.element != "H"]
    if atoms.array_length() == 0:
        raise ValueError("no heavy protein atoms found")

    per_atom_sasa = struc.sasa(
        atoms, probe_radius=1.4, point_number=1000, vdw_radii="ProtOr", ignore_ions=True
    )
    # NaN marks atoms SASA could not be computed for; they contribute nothing.
    per_atom_sasa = [0.0 if value != value else float(value) for value in per_atom_sasa]

    residues = {}
    order = []
    for index in range(atoms.array_length()):
        key = (
            str(atoms.chain_id[index]),
            int(atoms.res_id[index]),
            str(atoms.ins_code[index]).strip(),
        )
        if key not in residues:
            one_letter = THREE_TO_ONE.get(str(atoms.res_name[index]).upper().strip())
            residues[key] = {
                "chain": key[0],
                "res_id": key[1],
                "insertion": key[2],
                "aa": one_letter,
                "sasa": 0.0,
                "atoms": [],
                "named_atoms": {},
                "b_factors": [],
            }
            order.append(key)
        residue = residues[key]
        coordinate = (
            float(atoms.coord[index][0]),
            float(atoms.coord[index][1]),
            float(atoms.coord[index][2]),
        )
        atom_name = str(atoms.atom_name[index]).strip()
        residue["sasa"] += per_atom_sasa[index]
        residue["atoms"].append(coordinate)
        residue["named_atoms"].setdefault(atom_name, coordinate)

    return [residues[key] for key in order]


def centroid(coordinates):
    count = len(coordinates)
    return (
        sum(point[0] for point in coordinates) / count,
        sum(point[1] for point in coordinates) / count,
        sum(point[2] for point in coordinates) / count,
    )


def max_atom_radius(residue):
    origin = residue["centroid"]
    return max(math.dist(origin, point) for point in residue["atoms"])


def min_atom_distance(first, second, cutoff):
    """Closest heavy-atom distance, with a centroid pre-filter and early exit."""
    gap = math.dist(first["centroid"], second["centroid"]) - first["radius"] - second["radius"]
    if gap > cutoff:
        return None
    best = None
    for point in first["atoms"]:
        for other in second["atoms"]:
            distance = math.dist(point, other)
            if best is None or distance < best:
                best = distance
                if best < 1e-6:
                    return best
    return best


def find_salt_bridges(residues):
    """Indices of residues locked in a salt bridge, plus the pairs themselves."""
    cations = [
        index
        for index, residue in enumerate(residues)
        if residue["aa"] in CATION_ATOMS
    ]
    anions = [
        index for index, residue in enumerate(residues) if residue["aa"] in ANION_ATOMS
    ]

    bridged, pairs = set(), []
    for cation_index in cations:
        cation = residues[cation_index]
        cation_points = [
            point
            for name, point in cation["named_atoms"].items()
            if name in CATION_ATOMS[cation["aa"]]
        ]
        if not cation_points:
            continue
        for anion_index in anions:
            anion = residues[anion_index]
            anion_points = [
                point
                for name, point in anion["named_atoms"].items()
                if name in ANION_ATOMS[anion["aa"]]
            ]
            if not anion_points:
                continue
            closest = min(
                math.dist(point, other) for point in cation_points for other in anion_points
            )
            if closest <= SALT_BRIDGE_CUTOFF:
                bridged.add(cation_index)
                bridged.add(anion_index)
                pairs.append(
                    {
                        "cation": f"{cation['chain']}{cation['res_id']}{cation['insertion']}"
                        f"({cation['aa']})",
                        "anion": f"{anion['chain']}{anion['res_id']}{anion['insertion']}"
                        f"({anion['aa']})",
                        "distance": round(closest, 2),
                    }
                )
    return bridged, pairs


def build_cdr_vicinity(residues, exposed_indices):
    """Exposed CDR and anchor residues, plus exposed residues within 4 A of them."""
    seeds = [
        index
        for index in exposed_indices
        if residues[index]["region"].startswith("CDR") or residues[index]["is_anchor"]
    ]
    vicinity = set(seeds)
    seed_set = set(seeds)
    for index in exposed_indices:
        if index in seed_set:
            continue
        residue = residues[index]
        for seed in seeds:
            distance = min_atom_distance(residue, residues[seed], VICINITY_RADIUS)
            if distance is not None and distance <= VICINITY_RADIUS:
                vicinity.add(index)
                break
    return sorted(vicinity), sorted(seed_set)


def pair_sum(residues, indices, weight_key, require_positive=None):
    """The TAP pair sum: Sum w(R1)*w(R2)/r12^2 over pairs closer than 7.5 A.

    `require_positive` restricts contributions by charge sign for PPC/PNC; PSH passes None.
    """
    members = []
    for index in indices:
        residue = residues[index]
        weight = residue[weight_key]
        if require_positive is True and residue["charge"] <= 0:
            continue
        if require_positive is False and residue["charge"] >= 0:
            continue
        if weight == 0:
            continue
        members.append(index)

    total = 0.0
    contributions = 0
    for position, first_index in enumerate(members):
        first = residues[first_index]
        for second_index in members[position + 1 :]:
            second = residues[second_index]
            distance = min_atom_distance(first, second, PAIR_CUTOFF)
            if distance is None or distance >= PAIR_CUTOFF or distance <= 0:
                continue
            total += (first[weight_key] * second[weight_key]) / (distance * distance)
            contributions += 1
    return total, contributions


def flag_for(metric, value):
    """'green', 'amber' or 'red' against the Table 2 bounds; '' when value is unavailable."""
    if value is None:
        return ""
    bounds = TAP_THRESHOLDS[metric]
    red_above, red_below = bounds["red_above"], bounds["red_below"]
    if red_above is not None and value > red_above:
        return "red"
    if red_below is not None and value < red_below:
        return "red"
    for low, high in bounds["amber"]:
        if low <= value <= high:
            return "amber"
    return "green"


def main():
    os.makedirs("./outputs", exist_ok=True)

    cdr_definition = os.environ.get("PARAM_CDR_DEFINITION", "imgt")
    threshold_set = os.environ.get("PARAM_TAP_THRESHOLD_SET", "clinical")
    if threshold_set not in {"clinical", "none"}:
        fail(f"tap_threshold_set must be 'clinical' or 'none', got '{threshold_set}'")

    annotations = read_json("./inputs/annotations.json", "annotations.json")
    clones = {record["clone_id"]: record for record in annotations.get("clones") or []}
    if not clones:
        fail("annotations.json contains no annotated clones.")

    annotation_definition = (annotations.get("settings") or {}).get("cdr_definition", "imgt")
    if annotation_definition != cdr_definition:
        warn(
            f"cdr_definition is '{cdr_definition}' here but node 01 annotated with "
            f"'{annotation_definition}'. The regions in annotations.json win, because they are "
            "what the residues are labelled with. Set both to the same value."
        )
        cdr_definition = annotation_definition

    liabilities = []
    if os.path.exists("./inputs/liabilities.csv"):
        with open("./inputs/liabilities.csv", newline="") as handle:
            liabilities = list(csv.DictReader(handle))

    hydrophobicity = normalized_hydrophobicity()
    glycine_hydrophobicity = hydrophobicity["G"]

    implementation = {
        "sasa_algorithm": "Shrake-Rupley (biotite.structure.sasa)",
        "sasa_probe_radius_angstrom": 1.4,
        "sasa_point_number": 1000,
        "sasa_vdw_radii": "ProtOr",
        "sasa_atoms": "heavy atoms only (hydrogens removed before the calculation)",
        "exposure_numerator": "whole-residue SASA",
        "exposure_reference": (
            "Tien et al. 2013 theoretical Ala-X-Ala whole-residue max ASA (PMC3836772, Table 1)"
        ),
        "exposure_threshold_relative": EXPOSURE_THRESHOLD,
        "hydrophobicity_scale": "Kyte & Doolittle via Bio.SeqUtils.ProtParamData.kd",
        "hydrophobicity_normalization": "linear to [1, 2] across the scale's own min and max",
        "salt_bridge_rule": f"K/R N+ to D/E O- distance <= {SALT_BRIDGE_CUTOFF} A",
        "salt_bridge_hydrophobicity": "glycine's normalized value",
        "salt_bridge_charge": 0.0,
        "cdr_vicinity_rule": (
            f"exposed CDR and anchor residues, plus exposed residues with a heavy atom within "
            f"{VICINITY_RADIUS} A"
        ),
        "pair_cutoff_angstrom": PAIR_CUTOFF,
        "cdr_definition": cdr_definition,
    }
    calibration = {
        "threshold_source": (
            "Raybould et al. 2019, PNAS, Table 2 (242 clinical-stage therapeutics). "
            "DOI 10.1073/pnas.1810576116, CC BY."
        ),
        "threshold_set": threshold_set,
        "known_deviations": [
            "TAP uses side-chain relative exposure against a Yang Zhu & Blundell (1996) "
            "Ala-R-Ala reference; this node uses whole-residue relative ASA against Tien et al. "
            "(2013). The exposed residue set, and therefore every metric, differs.",
            "TAP's thresholds were calibrated on ABodyBuilder v1 models; these models come from "
            "ABodyBuilder2.",
            "SASA point density, van der Waals radii and the [1, 2] hydrophobicity "
            "normalization are not pinned by the paper.",
        ],
        "cst_reference": CST_REFERENCE,
        "interpretation": (
            "Flags are risk indicators, not predictions of clinical success. Amber means the "
            "value sits in the extremes of the clinical-stage distribution; red means a value "
            "not previously observed among those therapeutics. Approved therapeutics violate "
            "these guidelines."
        ),
    }
    print(f"CDR definition: {cdr_definition}; threshold set: {threshold_set}", flush=True)

    rows, detail, flag_rows, failures = [], {}, [], []
    psh_values = []

    for clone_id, record in clones.items():
        path = f"./inputs/models/{clone_id}.pdb"
        if not os.path.exists(path):
            failures.append({"clone_id": clone_id, "reason": "no model from node 02"})
            continue

        try:
            residues = load_structure(path)
        except Exception as exc:
            failures.append({"clone_id": clone_id, "reason": f"could not read model: {exc}"})
            print(f"  {clone_id}: FAILED — {exc}", flush=True)
            continue

        # Join structure residues to node 01's regions on (chain, position, insertion).
        region_lookup = {}
        for chain_label, chain in (record.get("chains") or {}).items():
            letter = CHAIN_FOR.get(chain_label, chain_label)
            for residue_row in chain.get("residues") or []:
                key = (letter, int(residue_row["position"]), residue_row.get("insertion", ""))
                region_lookup[key] = (
                    residue_row.get("region", "outside"),
                    bool(residue_row.get("is_anchor")),
                )

        unmatched = 0
        for residue in residues:
            key = (residue["chain"], residue["res_id"], residue["insertion"])
            region, is_anchor = region_lookup.get(key, ("outside", False))
            residue["region"] = region
            residue["is_anchor"] = is_anchor
            if region == "outside":
                unmatched += 1
            residue["centroid"] = centroid(residue["atoms"])
            residue["radius"] = max_atom_radius(residue)
            reference = TIEN_THEORETICAL_MAX_ASA.get(residue["aa"])
            residue["relative_exposure"] = (
                residue["sasa"] / reference if reference else None
            )
            residue["exposed"] = bool(
                residue["relative_exposure"] is not None
                and residue["relative_exposure"] >= EXPOSURE_THRESHOLD
            )

        if unmatched > 0.2 * len(residues):
            failures.append(
                {
                    "clone_id": clone_id,
                    "reason": (
                        f"{unmatched} of {len(residues)} model residues could not be matched to "
                        "node 01's numbering. The model and the annotation disagree, which "
                        "means the CDR assignment would be wrong."
                    ),
                }
            )
            print(f"  {clone_id}: FAILED — numbering mismatch with node 01", flush=True)
            continue

        bridged, bridge_pairs = find_salt_bridges(residues)
        for index, residue in enumerate(residues):
            in_bridge = index in bridged
            residue["salt_bridged"] = in_bridge
            if in_bridge:
                residue["charge"] = 0.0
                residue["hydrophobicity"] = glycine_hydrophobicity
            else:
                residue["charge"] = RESIDUE_CHARGE.get(residue["aa"], 0.0)
                residue["hydrophobicity"] = hydrophobicity.get(residue["aa"], 1.0)
            residue["abs_charge"] = abs(residue["charge"])

        exposed_indices = [
            index for index, residue in enumerate(residues) if residue["exposed"]
        ]
        vicinity, seeds = build_cdr_vicinity(residues, exposed_indices)

        psh, psh_pairs = pair_sum(residues, vicinity, "hydrophobicity")
        ppc, ppc_pairs = pair_sum(residues, vicinity, "abs_charge", require_positive=True)
        pnc, pnc_pairs = pair_sum(residues, vicinity, "abs_charge", require_positive=False)

        is_nanobody = record.get("receptor_type") == "nanobody"
        heavy_charge = sum(
            residues[index]["charge"]
            for index in exposed_indices
            if residues[index]["chain"] == "H"
        )
        light_charge = sum(
            residues[index]["charge"]
            for index in exposed_indices
            if residues[index]["chain"] == "L"
        )
        sfvcsp = None if is_nanobody else heavy_charge * light_charge

        total_cdr_length = record.get("total_cdr_length")

        # Isoelectric point of the modelled domains, for context rather than as a TAP metric.
        try:
            from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint

            sequence = "".join(
                chain["domain_sequence"] for chain in (record.get("chains") or {}).values()
            )
            isoelectric_point = round(IsoelectricPoint(sequence).pi(), 2)
        except Exception as exc:
            warn(f"{clone_id}: could not compute the isoelectric point: {exc}")
            isoelectric_point = None

        metrics = {
            "total_cdr_length": float(total_cdr_length) if total_cdr_length else None,
            "psh_cdr_vicinity": psh,
            "ppc_cdr_vicinity": ppc,
            "pnc_cdr_vicinity": pnc,
            "sfvcsp": sfvcsp,
        }
        flags = (
            {name: flag_for(name, value) for name, value in metrics.items()}
            if threshold_set == "clinical"
            else {name: "" for name in metrics}
        )
        if is_nanobody:
            flags["sfvcsp"] = ""
        if psh is not None:
            psh_values.append(psh)

        # Exposure of each liability hit: a buried sequon is a smaller problem than an exposed
        # one, and this is the only place the two sources of evidence meet.
        exposure_by_key = {
            (residue["chain"], residue["res_id"], residue["insertion"]): residue
            for residue in residues
        }
        exposed_liabilities = []
        for hit in liabilities:
            if hit.get("clone_id") != clone_id:
                continue
            letter = CHAIN_FOR.get(hit.get("chain"), hit.get("chain"))
            raw = (hit.get("position") or "").strip()
            digits = "".join(character for character in raw if character.isdigit())
            insertion = "".join(character for character in raw if character.isalpha())
            if not digits:
                continue
            residue = exposure_by_key.get((letter, int(digits), insertion))
            if residue is None:
                continue
            exposed_liabilities.append(
                {
                    "chain": hit.get("chain"),
                    "motif": hit.get("motif"),
                    "severity": hit.get("severity"),
                    "position": raw,
                    "region": hit.get("region"),
                    "relative_exposure": (
                        round(residue["relative_exposure"], 3)
                        if residue["relative_exposure"] is not None
                        else None
                    ),
                    "exposed": residue["exposed"],
                }
            )

        n_exposed_liabilities = sum(1 for hit in exposed_liabilities if hit["exposed"])
        n_exposed_high = sum(
            1
            for hit in exposed_liabilities
            if hit["exposed"] and hit["severity"] == "high"
        )

        rows.append(
            {
                "clone_id": clone_id,
                "receptor_type": record.get("receptor_type"),
                "total_cdr_length": total_cdr_length if total_cdr_length else "",
                "psh_cdr_vicinity": round(psh, 2),
                "ppc_cdr_vicinity": round(ppc, 2),
                "pnc_cdr_vicinity": round(pnc, 2),
                "sfvcsp": round(sfvcsp, 2) if sfvcsp is not None else "",
                "isoelectric_point": isoelectric_point if isoelectric_point is not None else "",
                "n_residues": len(residues),
                "n_exposed": len(exposed_indices),
                "n_cdr_vicinity": len(vicinity),
                "n_salt_bridges": len(bridge_pairs),
                "n_exposed_liabilities": n_exposed_liabilities,
                "n_exposed_high_liabilities": n_exposed_high,
                "flag_total_cdr_length": flags["total_cdr_length"],
                "flag_psh": flags["psh_cdr_vicinity"],
                "flag_ppc": flags["ppc_cdr_vicinity"],
                "flag_pnc": flags["pnc_cdr_vicinity"],
                "flag_sfvcsp": flags["sfvcsp"],
                "thresholds_calibrated": "no" if is_nanobody else "approximate",
                "status": "ok",
                "reason": "",
            }
        )

        detail[clone_id] = {
            "metrics": metrics,
            "flags": flags,
            "pair_counts": {"psh": psh_pairs, "ppc": ppc_pairs, "pnc": pnc_pairs},
            "chain_surface_charge": {"H": round(heavy_charge, 2), "L": round(light_charge, 2)},
            "salt_bridges": bridge_pairs,
            "cdr_vicinity_size": len(vicinity),
            "cdr_vicinity_seeds": len(seeds),
            "n_exposed": len(exposed_indices),
            "n_residues": len(residues),
            "unmatched_residues": unmatched,
            "exposed_liabilities": exposed_liabilities,
            "isoelectric_point": isoelectric_point,
            "nanobody_thresholds_not_calibrated": is_nanobody,
            "per_residue": [
                {
                    "chain": residue["chain"],
                    "position": f"{residue['res_id']}{residue['insertion']}",
                    "aa": residue["aa"],
                    "region": residue["region"],
                    "is_anchor": residue["is_anchor"],
                    "relative_exposure": (
                        round(residue["relative_exposure"], 3)
                        if residue["relative_exposure"] is not None
                        else None
                    ),
                    "exposed": residue["exposed"],
                    "salt_bridged": residue["salt_bridged"],
                    "charge": residue["charge"],
                    "in_cdr_vicinity": index in set(vicinity),
                }
                for index, residue in enumerate(residues)
            ],
        }

        flag_rows.append(
            {
                "clone_id": clone_id,
                **{f"flag_{name}": flags[name] for name in metrics},
                "n_amber": sum(1 for value in flags.values() if value == "amber"),
                "n_red": sum(1 for value in flags.values() if value == "red"),
            }
        )
        print(
            f"  {clone_id}: PSH {psh:.1f} ({flags['psh_cdr_vicinity'] or 'unflagged'}), "
            f"PPC {ppc:.2f}, PNC {pnc:.2f}, "
            f"SFvCSP {f'{sfvcsp:.1f}' if sfvcsp is not None else 'n/a'}, "
            f"CDR length {total_cdr_length}, {len(vicinity)} residues in the CDR vicinity",
            flush=True,
        )

    for failure in failures:
        rows.append(
            {
                "clone_id": failure["clone_id"],
                "status": "failed",
                "reason": failure["reason"],
            }
        )
        flag_rows.append({"clone_id": failure["clone_id"], "n_amber": "", "n_red": ""})

    fields = [
        "clone_id",
        "receptor_type",
        "total_cdr_length",
        "psh_cdr_vicinity",
        "ppc_cdr_vicinity",
        "pnc_cdr_vicinity",
        "sfvcsp",
        "isoelectric_point",
        "n_residues",
        "n_exposed",
        "n_cdr_vicinity",
        "n_salt_bridges",
        "n_exposed_liabilities",
        "n_exposed_high_liabilities",
        "flag_total_cdr_length",
        "flag_psh",
        "flag_ppc",
        "flag_pnc",
        "flag_sfvcsp",
        "thresholds_calibrated",
        "status",
        "reason",
    ]
    with open("./outputs/developability.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    flag_fields = [
        "clone_id",
        "flag_total_cdr_length",
        "flag_psh_cdr_vicinity",
        "flag_ppc_cdr_vicinity",
        "flag_pnc_cdr_vicinity",
        "flag_sfvcsp",
        "n_amber",
        "n_red",
    ]
    with open("./outputs/flags.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=flag_fields, extrasaction="ignore")
        writer.writeheader()
        for row in flag_rows:
            writer.writerow(row)

    panel_psh_mean = sum(psh_values) / len(psh_values) if psh_values else None
    offset_note = None
    if panel_psh_mean is not None:
        deviation = (panel_psh_mean - CST_REFERENCE["psh_cdr_vicinity_mean"]) / CST_REFERENCE[
            "psh_cdr_vicinity_sd"
        ]
        offset_note = (
            f"Panel mean CDR-vicinity PSH is {panel_psh_mean:.1f}, which is "
            f"{deviation:+.1f} standard deviations from the published CST mean of "
            f"{CST_REFERENCE['psh_cdr_vicinity_mean']} +/- {CST_REFERENCE['psh_cdr_vicinity_sd']}."
        )

    payload = {
        "implementation": implementation,
        "calibration": {**calibration, "panel_psh_offset": offset_note},
        "clones": detail,
        "failures": failures,
        "n_scored": len(detail),
        "metric_definitions": {
            "psh_cdr_vicinity": (
                "Sum over CDR-vicinity surface-exposed residue pairs closer than 7.5 A of "
                "H(R1)*H(R2)/r12^2, Kyte & Doolittle normalized to [1,2]. TAP-style."
            ),
            "ppc_cdr_vicinity": "As PSH with |Q(R)| in place of H, positive residues only.",
            "pnc_cdr_vicinity": "As PSH with |Q(R)| in place of H, negative residues only.",
            "sfvcsp": (
                "Product of the summed surface charges of the VH and VL chains. Undefined for "
                "nanobodies."
            ),
            "total_cdr_length": "Sum of all six CDR lengths under the chosen CDR definition.",
        },
    }
    with open("./outputs/developability.json", "w") as handle:
        json.dump(payload, handle, indent=2)

    if not detail:
        fail(
            f"No clone could be scored ({len(failures)} failure(s)). See developability.csv "
            "for the per-clone reasons."
        )

    print(f"Scored {len(detail)} clone(s); {len(failures)} failure(s)", flush=True)
    if offset_note:
        print(offset_note, flush=True)
        if abs(panel_psh_mean - CST_REFERENCE["psh_cdr_vicinity_mean"]) > 3 * CST_REFERENCE[
            "psh_cdr_vicinity_sd"
        ]:
            warn(
                "The panel's PSH mean is more than 3 SD from the published CST mean. That is "
                "either a genuinely unusual panel or a calibration problem in this "
                "implementation (most likely the exposure reference). Treat the PSH flags as "
                "unreliable until it is explained — see 'the calibration problem' in the "
                "pipeline doc."
            )
    if threshold_set == "clinical":
        print(
            "Flags are TAP-style risk indicators against Table 2 of Raybould et al. 2019, not "
            "predictions of clinical success, and they carry a documented calibration offset.",
            flush=True,
        )


if __name__ == "__main__":
    main()
