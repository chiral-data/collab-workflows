#!/usr/bin/env python3
"""Dependency-free unit tests for workflow-039's node logic.

These cover the parts that can be checked without ANARCI, ImmuneBuilder, AbLang2, promb,
biotite or biopython installed: CDR and anchor derivation, liability motif scanning, the
hydrophobicity normalization, the TAP pair sums on hand-built geometry, salt-bridge detection,
every Table 2 flag boundary, the mutation re-screen rules, and the ranking order.

They are not a substitute for running the pipeline — see the verification plan in
antibody_triage_pipeline.md for what still needs a built image. What they do catch is the
arithmetic and the boundary conditions, which is where a reimplementation of someone else's
published metric goes wrong quietly.

Run: python3 tests/test_logic.py
"""
import importlib.util
import math
import os
import sys
import types

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)

FAILURES = []
CHECKS = [0]


def load(node, name):
    """Import a node's script.py under its own module name."""
    path = os.path.join(ROOT, node, "script.py")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def check(label, actual, expected):
    CHECKS[0] += 1
    ok = actual == expected
    if not ok:
        FAILURES.append(f"{label}: expected {expected!r}, got {actual!r}")
    print(f"{'ok  ' if ok else 'FAIL'} {label}")


def check_close(label, actual, expected, tolerance=1e-9):
    CHECKS[0] += 1
    ok = actual is not None and abs(actual - expected) <= tolerance
    if not ok:
        FAILURES.append(f"{label}: expected {expected!r} +/- {tolerance}, got {actual!r}")
    print(f"{'ok  ' if ok else 'FAIL'} {label}")


def install_fake_biopython():
    """Minimal stand-ins so node 03/05 logic can be imported without biopython.

    The Kyte & Doolittle values here are only a fixture for testing the normalization maths;
    the real nodes read the scale from Bio.SeqUtils.ProtParamData so the numbers are the
    library's, not ours.
    """
    bio = types.ModuleType("Bio")
    sequtils = types.ModuleType("Bio.SeqUtils")
    protparamdata = types.ModuleType("Bio.SeqUtils.ProtParamData")
    protparamdata.kd = {
        "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
        "Q": -3.5, "E": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
        "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
        "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
    }
    isoelectric = types.ModuleType("Bio.SeqUtils.IsoelectricPoint")

    class IsoelectricPoint:
        def __init__(self, sequence):
            self.sequence = sequence

        def pi(self):
            return 7.0

    isoelectric.IsoelectricPoint = IsoelectricPoint
    sequtils.ProtParamData = protparamdata
    sequtils.IsoelectricPoint = isoelectric
    bio.SeqUtils = sequtils
    sys.modules.setdefault("Bio", bio)
    sys.modules.setdefault("Bio.SeqUtils", sequtils)
    sys.modules.setdefault("Bio.SeqUtils.ProtParamData", protparamdata)
    sys.modules.setdefault("Bio.SeqUtils.IsoelectricPoint", isoelectric)


def residue(aa, coordinates, **extra):
    """A node-03-shaped residue dict with the derived geometry filled in."""
    entry = {
        "chain": extra.get("chain", "H"),
        "res_id": extra.get("res_id", 1),
        "insertion": extra.get("insertion", ""),
        "aa": aa,
        "atoms": list(coordinates),
        "named_atoms": extra.get("named_atoms", {}),
        "region": extra.get("region", "FR1"),
        "is_anchor": extra.get("is_anchor", False),
        "exposed": extra.get("exposed", True),
        "charge": extra.get("charge", 0.0),
        "hydrophobicity": extra.get("hydrophobicity", 1.0),
        "salt_bridged": extra.get("salt_bridged", False),
    }
    entry["abs_charge"] = abs(entry["charge"])
    count = len(entry["atoms"])
    entry["centroid"] = (
        sum(point[0] for point in entry["atoms"]) / count,
        sum(point[1] for point in entry["atoms"]) / count,
        sum(point[2] for point in entry["atoms"]) / count,
    )
    entry["radius"] = max(math.dist(entry["centroid"], point) for point in entry["atoms"])
    return entry


# --------------------------------------------------------------------------- node 01


def test_annotate(annotate):
    print("\n== node 01: regions, anchors, motifs ==")

    imgt = annotate.borders_for("imgt", "H")
    check("imgt borders are scheme-wide", imgt, [27, 39, 56, 66, 105, 118, 129])
    check(
        "kabat borders differ by chain",
        (annotate.borders_for("kabat", "H"), annotate.borders_for("kabat", "K")),
        ([31, 36, 50, 66, 95, 103, 114], [24, 35, 50, 57, 89, 98, 108]),
    )
    check(
        "north is defined on chothia numbering",
        annotate.CDR_DEFINITION_SCHEME["north"],
        "chothia",
    )

    regions, anchors = annotate.regions_from_borders(imgt)
    check("IMGT CDR1 starts at 27", regions[27], "CDR1")
    check("IMGT CDR1 ends at 38", regions[38], "CDR1")
    check("IMGT 39 is framework", regions[39], "FR2")
    check("IMGT CDR3 spans 105-117", (regions[105], regions[117]), ("CDR3", "CDR3"))
    check("IMGT 118 is framework", regions[118], "FR4")
    # TAP's CDR vicinity includes "CDR and anchor residues"; anchors flank each CDR.
    check("IMGT anchors are the flanking positions", anchors, [26, 39, 55, 66, 104, 118])

    rows = []
    sequence = "QVQLNGSCDPMWNAT"
    for index, aa in enumerate(sequence):
        rows.append(
            {"seq_index": index, "aa": aa, "position": 27 + index, "insertion": ""}
        )
    hits = annotate.scan_liabilities(
        "c1", "VH", rows, annotate.DEFAULT_MOTIFS, regions
    )
    names = sorted({hit["motif"] for hit in hits})
    # The fixture deliberately contains overlapping motifs: NGS and NAT are sequons, NG is
    # a deamidation site, NA is a slower one, SC is a fragmentation motif, DP is acid-labile,
    # and the single cysteine is unpaired. All of them should be reported.
    check(
        "planted motifs are all found, including the overlapping ones",
        names,
        sorted(
            {
                "n_glycosylation",
                "deamidation_ng",
                "deamidation_ns",
                "fragmentation_cleavage",
                "acid_cleavage_dp",
                "oxidation_met",
                "oxidation_trp",
                "unpaired_cysteine",
            }
        ),
    )
    glyco = [hit for hit in hits if hit["motif"] == "n_glycosylation"]
    check("two N-X-S/T sequons (NGS and NAT)", len(glyco), 2)
    check("sequon is reported at its numbered position", glyco[0]["position"], "31")
    check("sequon is reported in its region", glyco[0]["region"], "CDR1")
    check("CDR membership is flagged", glyco[0]["in_cdr"], True)

    even = [
        {"seq_index": i, "aa": aa, "position": 27 + i, "insertion": ""}
        for i, aa in enumerate("CAAC")
    ]
    paired = annotate.scan_liabilities("c2", "VH", even, [], regions)
    check("an even cysteine count is not flagged", paired, [])
    odd = [
        {"seq_index": i, "aa": aa, "position": 27 + i, "insertion": ""}
        for i, aa in enumerate("CAA")
    ]
    unpaired = annotate.scan_liabilities("c3", "VH", odd, [], regions)
    check("an odd cysteine count is flagged once", len(unpaired), 1)
    check("and flagged as unpaired", unpaired[0]["motif"], "unpaired_cysteine")

    # The sequence-index mapping is what joins two numbering schemes and the structure.
    table = annotate.residue_table(
        ([((1, " "), "E"), ((2, " "), "V"), ((3, " "), "-"), ((4, " "), "Q")], 2, 5),
        "XXEVQ",
    )
    check(
        "numbering maps back onto input indices, skipping gaps",
        [(row["seq_index"], row["aa"], row["position"]) for row in table],
        [(2, "E", 1), (3, "V", 2), (4, "Q", 4)],
    )


# --------------------------------------------------------------------------- node 03


def test_developability(dev):
    print("\n== node 03: TAP-style metrics ==")

    hydrophobicity = dev.normalized_hydrophobicity()
    check_close("most hydrophilic residue normalizes to 1.0", hydrophobicity["R"], 1.0)
    check_close("most hydrophobic residue normalizes to 2.0", hydrophobicity["I"], 2.0)
    check("every normalized value is inside [1, 2]",
          all(1.0 <= value <= 2.0 for value in hydrophobicity.values()), True)
    # Glycine's value is what salt-bridged residues inherit, per the paper.
    check_close(
        "glycine normalizes as expected",
        hydrophobicity["G"],
        1.0 + (-0.4 - -4.5) / (4.5 - -4.5),
    )

    check("charges match the paper's table",
          (dev.RESIDUE_CHARGE["D"], dev.RESIDUE_CHARGE["E"], dev.RESIDUE_CHARGE["K"],
           dev.RESIDUE_CHARGE["R"], dev.RESIDUE_CHARGE["H"]),
          (-1.0, -1.0, 1.0, 1.0, 0.1))

    # PSH on a hand-built geometry: two single-atom residues 3 A apart.
    first = residue("I", [(0.0, 0.0, 0.0)], hydrophobicity=2.0)
    second = residue("L", [(3.0, 0.0, 0.0)], hydrophobicity=1.5)
    total, pairs = dev.pair_sum([first, second], [0, 1], "hydrophobicity")
    check_close("PSH of one pair is H1*H2/r^2", total, (2.0 * 1.5) / 9.0)
    check("and counts one contributing pair", pairs, 1)

    far = residue("L", [(8.0, 0.0, 0.0)], hydrophobicity=1.5)
    total_far, pairs_far = dev.pair_sum([first, far], [0, 1], "hydrophobicity")
    check("pairs beyond 7.5 A contribute nothing", (total_far, pairs_far), (0.0, 0))

    # Three residues: the sum runs over all pairs within the cutoff.
    third = residue("V", [(0.0, 4.0, 0.0)], hydrophobicity=1.8)
    total_three, pairs_three = dev.pair_sum(
        [first, second, third], [0, 1, 2], "hydrophobicity"
    )
    expected = (2.0 * 1.5) / 9.0 + (2.0 * 1.8) / 16.0 + (1.5 * 1.8) / 25.0
    check_close("PSH sums every pair inside the cutoff", total_three, expected)
    check("and counts three pairs", pairs_three, 3)

    # PPC/PNC use |Q| and are restricted by sign.
    lysine = residue("K", [(0.0, 0.0, 0.0)], charge=1.0)
    arginine = residue("R", [(3.0, 0.0, 0.0)], charge=1.0)
    aspartate = residue("D", [(0.0, 3.0, 0.0)], charge=-1.0)
    residues = [lysine, arginine, aspartate]
    ppc, ppc_pairs = dev.pair_sum(residues, [0, 1, 2], "abs_charge", require_positive=True)
    check_close("PPC counts only the positive pair", ppc, 1.0 / 9.0)
    check("PPC pair count", ppc_pairs, 1)
    pnc, pnc_pairs = dev.pair_sum(residues, [0, 1, 2], "abs_charge", require_positive=False)
    check("PNC needs two negatives, so a lone aspartate scores 0", (pnc, pnc_pairs), (0.0, 0))

    glutamate = residue("E", [(0.0, 5.0, 0.0)], charge=-1.0)
    pnc_two, pnc_two_pairs = dev.pair_sum(
        [aspartate, glutamate], [0, 1], "abs_charge", require_positive=False
    )
    check_close("PNC of two negatives 2 A apart", pnc_two, 1.0 / 4.0)
    check("PNC pair count", pnc_two_pairs, 1)

    # Salt bridges: K/R N+ to D/E O- within 3.2 A.
    bridged_k = residue("K", [(0.0, 0.0, 0.0)], named_atoms={"NZ": (0.0, 0.0, 0.0)})
    bridged_d = residue("D", [(3.0, 0.0, 0.0)], named_atoms={"OD1": (3.0, 0.0, 0.0)})
    indices, bridges = dev.find_salt_bridges([bridged_k, bridged_d])
    check("a 3.0 A N+-O- contact is a salt bridge", sorted(indices), [0, 1])
    check_close("and its distance is recorded", bridges[0]["distance"], 3.0)

    distant_d = residue("D", [(3.3, 0.0, 0.0)], named_atoms={"OD1": (3.3, 0.0, 0.0)})
    indices_far, bridges_far = dev.find_salt_bridges([bridged_k, distant_d])
    check("3.3 A is beyond the 3.2 A cutoff", (sorted(indices_far), bridges_far), ([], []))

    # Wrong atoms: a lysine backbone N does not make a bridge.
    wrong_atom = residue("K", [(0.0, 0.0, 0.0)], named_atoms={"N": (0.0, 0.0, 0.0)})
    indices_wrong, _ = dev.find_salt_bridges([wrong_atom, bridged_d])
    check("only side-chain N+ atoms count", sorted(indices_wrong), [])

    # CDR vicinity: exposed CDR/anchor residues plus exposed residues within 4 A.
    cdr = residue("Y", [(0.0, 0.0, 0.0)], region="CDR3", exposed=True)
    near = residue("S", [(3.0, 0.0, 0.0)], region="FR3", exposed=True)
    far_fr = residue("S", [(9.0, 0.0, 0.0)], region="FR3", exposed=True)
    buried = residue("L", [(3.2, 0.0, 0.0)], region="FR3", exposed=False)
    anchor = residue("W", [(20.0, 0.0, 0.0)], region="FR2", exposed=True, is_anchor=True)
    pool = [cdr, near, far_fr, buried, anchor]
    exposed_indices = [0, 1, 2, 4]
    vicinity, seeds = dev.build_cdr_vicinity(pool, exposed_indices)
    check("vicinity seeds are the CDR and anchor residues", seeds, [0, 4])
    check("vicinity adds the exposed neighbour within 4 A", vicinity, [0, 1, 4])
    check("a buried neighbour is excluded", 3 in vicinity, False)
    check("a distant framework residue is excluded", 2 in vicinity, False)

    print("\n== node 03: Table 2 flag boundaries ==")
    # Total CDR length: amber 54-60 inclusive, red above 60.
    for value, expected in [
        (53.0, "green"), (54.0, "amber"), (60.0, "amber"), (60.1, "red"), (61.0, "red")
    ]:
        check(f"total_cdr_length {value} -> {expected}",
              dev.flag_for("total_cdr_length", value), expected)
    # PSH has two amber bands and red on both sides.
    for value, expected in [
        (83.83, "red"), (83.84, "amber"), (100.71, "amber"), (100.72, "green"),
        (123.30, "green"), (156.19, "green"), (156.20, "amber"), (173.85, "amber"),
        (173.86, "red"),
    ]:
        check(f"psh {value} -> {expected}", dev.flag_for("psh_cdr_vicinity", value), expected)
    for value, expected in [(1.24, "green"), (1.25, "amber"), (3.16, "amber"), (3.17, "red")]:
        check(f"ppc {value} -> {expected}", dev.flag_for("ppc_cdr_vicinity", value), expected)
    for value, expected in [(1.83, "green"), (1.84, "amber"), (3.50, "amber"), (3.51, "red")]:
        check(f"pnc {value} -> {expected}", dev.flag_for("pnc_cdr_vicinity", value), expected)
    # SFvCSP is flagged on the low side only.
    for value, expected in [
        (-20.41, "red"), (-20.40, "amber"), (-6.30, "amber"), (-6.29, "green"), (10.0, "green")
    ]:
        check(f"sfvcsp {value} -> {expected}", dev.flag_for("sfvcsp", value), expected)
    check("an unavailable value is not flagged", dev.flag_for("sfvcsp", None), "")

    check(
        "Tien reference covers all twenty residues",
        len(dev.TIEN_THEORETICAL_MAX_ASA),
        20,
    )
    check(
        "exposure threshold is the paper's 7.5%",
        dev.EXPOSURE_THRESHOLD,
        0.075,
    )
    check("vicinity radius is 4 A", dev.VICINITY_RADIUS, 4.0)
    check("pair cutoff is 7.5 A", dev.PAIR_CUTOFF, 7.5)
    check("salt-bridge cutoff is 3.2 A", dev.SALT_BRIDGE_CUTOFF, 3.2)


# --------------------------------------------------------------------------- node 05


def test_mutate_scan(scan):
    print("\n== node 05: proposal re-screening ==")

    motifs = [
        {"name": "n_glycosylation", "pattern": "N[^P][ST]"},
        {"name": "deamidation_ng", "pattern": "NG"},
    ]
    # G -> N completes N-S-T, which matches N[^P][ST].
    check(
        "a substitution completing a sequon is caught",
        scan.introduced_motifs("AAGSTA", "AANSTA", motifs),
        ["n_glycosylation"],
    )
    check(
        "a harmless substitution introduces nothing",
        scan.introduced_motifs("AAGSAA", "AAGSAV", motifs),
        [],
    )
    check(
        "a motif already present is not counted as introduced",
        scan.introduced_motifs("ANGSAA", "ANGSAV", motifs),
        [],
    )
    check(
        "NG created by the mutation is caught",
        scan.introduced_motifs("AAGAA", "ANGAA", motifs),
        ["deamidation_ng"],
    )

    check(
        "introducing a lone cysteine is caught",
        scan.creates_unpaired_cysteine("AAA", "ACA"),
        True,
    )
    check(
        "pairing up an already-unpaired cysteine is allowed, not blocked",
        scan.creates_unpaired_cysteine("ACA", "CCA"),
        False,
    )
    check(
        "removing a cysteine from an even count is caught",
        scan.creates_unpaired_cysteine("ACCA", "ACAA"),
        True,
    )

    kd = {"K": -3.9, "A": 1.8, "D": -3.5, "I": 4.5, "S": -0.8}
    check(
        "adding positive charge is refused on a PPC-flagged clone",
        scan.developability_direction("A", "K", {"ppc_cdr_vicinity": "amber"}, kd) is not None,
        True,
    )
    check(
        "removing positive charge is allowed on a PPC-flagged clone",
        scan.developability_direction("K", "A", {"ppc_cdr_vicinity": "amber"}, kd),
        None,
    )
    check(
        "adding negative charge is refused on a PNC-flagged clone",
        scan.developability_direction("A", "D", {"pnc_cdr_vicinity": "red"}, kd) is not None,
        True,
    )
    check(
        "adding hydrophobicity is refused on a PSH-red clone",
        scan.developability_direction("S", "I", {"psh_cdr_vicinity": "red"}, kd) is not None,
        True,
    )
    check(
        "the same substitution is allowed when PSH is only amber",
        scan.developability_direction("S", "I", {"psh_cdr_vicinity": "amber"}, kd),
        None,
    )
    check(
        "an unflagged clone is unconstrained",
        scan.developability_direction("A", "K", {}, kd),
        None,
    )

    residues = [
        {"region": "FR1"}, {"region": "CDR1"}, {"region": "FR2"},
        {"region": "CDR3"}, {"region": "FR4"},
    ]
    check("scan_region=cdr selects the CDRs",
          scan.scannable_positions("VH", residues, "cdr"), [1, 3])
    check("scan_region=cdrh3 selects only heavy CDR3",
          scan.scannable_positions("VH", residues, "cdrh3"), [3])
    check("scan_region=cdrh3 selects nothing on the light chain",
          scan.scannable_positions("VL", residues, "cdrh3"), [])
    check("scan_region=all selects everything",
          scan.scannable_positions("VH", residues, "all"), [0, 1, 2, 3, 4])


# --------------------------------------------------------------------------- node 06


def test_rank(rank):
    print("\n== node 06: ranking order ==")

    rows = [
        {"clone_id": "amber_only", "status": "ok", "n_red_flags": "0",
         "n_exposed_high_liabilities": "0", "n_amber_flags": "2",
         "mean_peptide_content": "0.80"},
        {"clone_id": "clean_less_human", "status": "ok", "n_red_flags": "0",
         "n_exposed_high_liabilities": "0", "n_amber_flags": "0",
         "mean_peptide_content": "0.50"},
        {"clone_id": "clean_more_human", "status": "ok", "n_red_flags": "0",
         "n_exposed_high_liabilities": "0", "n_amber_flags": "0",
         "mean_peptide_content": "0.90"},
        {"clone_id": "red_flag", "status": "ok", "n_red_flags": "1",
         "n_exposed_high_liabilities": "0", "n_amber_flags": "0",
         "mean_peptide_content": "0.99"},
        {"clone_id": "exposed_liability", "status": "ok", "n_red_flags": "0",
         "n_exposed_high_liabilities": "2", "n_amber_flags": "0",
         "mean_peptide_content": "0.95"},
        {"clone_id": "failed_clone", "status": "failed", "n_red_flags": "",
         "n_exposed_high_liabilities": "", "n_amber_flags": "",
         "mean_peptide_content": ""},
    ]
    ordered = [row["clone_id"] for row in sorted(rows, key=rank.main_sort_key)]
    check(
        "red flags then exposed liabilities then amber; failures last overall",
        ordered,
        [
            "clean_more_human",
            "clean_less_human",
            "amber_only",
            "exposed_liability",
            "red_flag",
            "failed_clone",
        ],
    )
    check("to_int tolerates blanks", rank.to_int("", 99), 99)
    check("to_float tolerates blanks", rank.to_float(""), None)
    check("to_float parses numbers", rank.to_float("1.5"), 1.5)


def main():
    install_fake_biopython()
    annotate = load("01-annotate", "wf039_annotate")
    dev = load("03-developability", "wf039_developability")
    scan = load("05-mutate-scan", "wf039_scan")
    rank = load("06-rank", "wf039_rank")

    test_annotate(annotate)
    test_developability(dev)
    test_mutate_scan(scan)
    test_rank(rank)

    print(f"\n{CHECKS[0]} check(s) run")
    if FAILURES:
        print(f"{len(FAILURES)} FAILURE(S):")
        for failure in FAILURES:
            print(f"  - {failure}")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
