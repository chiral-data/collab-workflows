#!/usr/bin/env python3
"""
Node 6: Comparison with wet-lab structure.

This script is intentionally tolerant to non-standard DiffDock-style PDB lines,
especially malformed coordinate columns in predicted poses.
"""

import argparse
import json
import logging
import os
import re
import shutil
import sys
from datetime import datetime

import numpy as np

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


def parse_chain_list(chain_csv):
    return [c.strip() for c in chain_csv.split(",") if c.strip()]


def _parse_float_triplet_tolerant(line):
    """
    Parse xyz from a possibly malformed PDB ATOM line.

    Tries strict PDB columns first, then falls back to regex extraction in the
    coordinate/occupancy region.
    """
    if len(line) >= 54:
        xs = line[30:38].strip()
        ys = line[38:46].strip()
        zs = line[46:54].strip()
        try:
            return np.array([float(xs), float(ys), float(zs)], dtype=float)
        except ValueError:
            pass

    # Fallback for lines like: "... 1594.864 575.259-558.879  1.00 ..."
    tail = line[26:] if len(line) > 26 else line
    nums = re.findall(r"[-+]?\d+\.\d+|[-+]?\d+", tail)
    if len(nums) < 3:
        return None
    try:
        return np.array([float(nums[0]), float(nums[1]), float(nums[2])], dtype=float)
    except ValueError:
        return None


def _parse_ca_line_tolerant(line, chain_index_counters):
    if not (line.startswith("ATOM") or line.startswith("HETATM")):
        return None

    atom_name = line[12:16].strip() if len(line) >= 16 else ""
    if not atom_name:
        toks = line.split()
        if len(toks) >= 3:
            atom_name = toks[2]

    if atom_name != "CA":
        return None

    # Try strict fields, then token fallback.
    resname = line[17:20].strip() if len(line) >= 20 else ""
    chain_id = line[21].strip() if len(line) >= 22 else ""
    resseq = None

    if len(line) >= 26:
        rs = line[22:26].strip()
        if rs:
            try:
                resseq = int(rs)
            except ValueError:
                resseq = None

    if not resname or not chain_id:
        toks = line.split()
        # Common fallback tokenization: ATOM serial CA RES CHAIN ...
        if len(toks) >= 5:
            if not resname:
                resname = toks[3]
            if not chain_id:
                chain_id = toks[4][:1]

    coord = _parse_float_triplet_tolerant(line)
    if coord is None or not chain_id:
        return None

    idx = chain_index_counters.get(chain_id, 0) + 1
    chain_index_counters[chain_id] = idx

    key_res = str(resseq) if resseq is not None else f"idx{idx}"
    key = f"{chain_id}:{key_res}"

    return {
        "chain": chain_id,
        "resname": resname or "UNK",
        "resseq": resseq,
        "index": idx,
        "key": key,
        "coord": coord,
    }


def parse_ca_atoms_tolerant(pdb_path, chain_ids=None):
    chain_ids = set(chain_ids or [])
    by_chain = {}
    all_chains = set()
    counters = {}

    with open(pdb_path, "r", errors="ignore") as f:
        for line in f:
            entry = _parse_ca_line_tolerant(line, counters)
            if entry is None:
                continue

            all_chains.add(entry["chain"])
            if chain_ids and entry["chain"] not in chain_ids:
                continue

            by_chain.setdefault(entry["chain"], []).append(entry)

    return by_chain, all_chains


def pair_chain_entries(ref_by_chain, pred_by_chain, chain_ids):
    """
    Pair residues per chain.

    Preferred: match by residue number when both chains are sufficiently annotated.
    Fallback: pair by sequence order and trim to min length.
    """
    pairs = []
    stats = {
        "chains_used": [],
        "matched": 0,
        "trimmed_by_order": 0,
        "matched_by_resseq": 0,
        "missing_chains_ref": [],
        "missing_chains_pred": [],
    }

    for chain in chain_ids:
        ref_list = ref_by_chain.get(chain, [])
        pred_list = pred_by_chain.get(chain, [])

        if not ref_list:
            stats["missing_chains_ref"].append(chain)
            continue
        if not pred_list:
            stats["missing_chains_pred"].append(chain)
            continue

        stats["chains_used"].append(chain)

        ref_annot = sum(1 for r in ref_list if r["resseq"] is not None)
        pred_annot = sum(1 for r in pred_list if r["resseq"] is not None)
        use_resseq = (
            ref_annot / len(ref_list) >= 0.8 and pred_annot / len(pred_list) >= 0.8
        )

        if use_resseq:
            ref_map = {r["resseq"]: r for r in ref_list if r["resseq"] is not None}
            pred_map = {r["resseq"]: r for r in pred_list if r["resseq"] is not None}
            common = sorted(set(ref_map.keys()) & set(pred_map.keys()))
            for k in common:
                pairs.append((ref_map[k], pred_map[k]))
            stats["matched_by_resseq"] += len(common)
        else:
            n = min(len(ref_list), len(pred_list))
            for i in range(n):
                pairs.append((ref_list[i], pred_list[i]))
            stats["trimmed_by_order"] += abs(len(ref_list) - len(pred_list))

        stats["matched"] = len(pairs)

    return pairs, stats


def kabsch_align(mobile_coords, target_coords):
    """Return R, t, aligned_mobile where aligned = mobile @ R + t."""
    if mobile_coords.shape != target_coords.shape:
        raise ValueError("Coordinate arrays must have same shape")
    if mobile_coords.shape[0] < 3:
        raise ValueError("Need at least 3 atom pairs for rigid alignment")

    cm = np.mean(mobile_coords, axis=0)
    ct = np.mean(target_coords, axis=0)

    mobile_centered = mobile_coords - cm
    target_centered = target_coords - ct

    h = mobile_centered.T @ target_centered
    u, _, vt = np.linalg.svd(h)
    r = vt.T @ u.T

    # Ensure proper rotation.
    if np.linalg.det(r) < 0:
        vt[-1, :] *= -1
        r = vt.T @ u.T

    t = ct - (cm @ r)
    aligned = mobile_coords @ r + t
    return r, t, aligned


def rmsd(a, b):
    if a.shape != b.shape or a.shape[0] == 0:
        return None
    d2 = np.sum((a - b) ** 2, axis=1)
    return float(np.sqrt(np.mean(d2)))


def interface_keys(ref_ab_coords, ref_ag_entries, cutoff=5.0):
    """Interface antigen residue keys in reference structure."""
    if len(ref_ab_coords) == 0:
        return set()

    ab = np.asarray(ref_ab_coords)
    keys = set()
    for e in ref_ag_entries:
        # Vectorized distance to all antibody CA atoms.
        d = np.linalg.norm(ab - e["coord"], axis=1)
        if np.any(d < cutoff):
            keys.add(e["key"])
    return keys


def interpret_rmsd(rmsd_value):
    if rmsd_value is None:
        return "unknown", "RMSD could not be computed"
    if rmsd_value < 2.0:
        return "excellent", "Excellent match - predicted pose is very close to experimental structure"
    if rmsd_value < 5.0:
        return "good", "Good match - reasonable agreement with experimental structure"
    if rmsd_value < 10.0:
        return "moderate", "Moderate match - significant conformational differences"
    return "poor", "Poor match - predicted pose differs substantially from experimental structure"


def generate_pymol_script(output_dir, ref_pdb, pred_pdb, ab_chains, ag_chains):
    script_path = os.path.join(output_dir, "align_structures.pml")
    script = f"""# PyMOL alignment script
load {os.path.basename(ref_pdb)}, ref_complex
load {os.path.basename(pred_pdb)}, pred_pose

select antibody_ref, (ref_complex and chain {'+'.join(ab_chains)})
select antibody_pred, (pred_pose and chain {'+'.join(ab_chains)})
select antigen_ref, (ref_complex and chain {'+'.join(ag_chains)})
select antigen_pred, (pred_pose and chain {'+'.join(ag_chains)})

align antibody_pred, antibody_ref
color red, ref_complex
color cyan, pred_pose
show cartoon, antibody_ref
show cartoon, antibody_pred
show surface, antigen_ref
show surface, antigen_pred
set transparency, 0.3
zoom
"""
    with open(script_path, "w") as f:
        f.write(script)
    return script_path


def calc_rmsd_tolerant(ref_complex, pred_pose, ab_chains, ag_chains):
    ref_ab_by_chain, ref_all = parse_ca_atoms_tolerant(ref_complex, ab_chains)
    pred_ab_by_chain, pred_all = parse_ca_atoms_tolerant(pred_pose, ab_chains)

    ab_pairs, ab_stats = pair_chain_entries(ref_ab_by_chain, pred_ab_by_chain, ab_chains)

    fallback_note = None
    superposition_mode = "antibody_superposition"
    if len(ab_pairs) >= 3:
        ref_ab_coords = np.asarray([p[0]["coord"] for p in ab_pairs], dtype=float)
        pred_ab_coords = np.asarray([p[1]["coord"] for p in ab_pairs], dtype=float)
        r, t, _ = kabsch_align(pred_ab_coords, ref_ab_coords)
    else:
        # Some DiffDock outputs contain ligand-only poses; in that case we cannot
        # apply the intended antibody-based superposition.
        superposition_mode = "none_missing_antibody"
        fallback_note = (
            "Predicted pose does not contain enough antibody CA anchors for superposition; "
            "antigen RMSD was computed in the original coordinate frame (no rigid alignment)."
        )
        logger.warning(fallback_note)
        r = np.eye(3, dtype=float)
        t = np.zeros(3, dtype=float)

    ref_ag_by_chain, _ = parse_ca_atoms_tolerant(ref_complex, ag_chains)
    pred_ag_by_chain, _ = parse_ca_atoms_tolerant(pred_pose, ag_chains)
    ag_pairs, ag_stats = pair_chain_entries(ref_ag_by_chain, pred_ag_by_chain, ag_chains)

    if len(ag_pairs) == 0:
        raise RuntimeError("No matched antigen CA atoms found")

    ref_ag_coords = np.asarray([p[0]["coord"] for p in ag_pairs], dtype=float)
    pred_ag_coords = np.asarray([p[1]["coord"] for p in ag_pairs], dtype=float)
    pred_ag_aligned = pred_ag_coords @ r + t

    full_rmsd = rmsd(ref_ag_coords, pred_ag_aligned)

    ref_ag_entries_all = []
    for c in ag_chains:
        ref_ag_entries_all.extend(ref_ag_by_chain.get(c, []))
    ref_ab_entries_all = []
    for c in ab_chains:
        ref_ab_entries_all.extend(ref_ab_by_chain.get(c, []))

    iface_keys = interface_keys(
        [e["coord"] for e in ref_ab_entries_all],
        ref_ag_entries_all,
        cutoff=5.0,
    )

    iface_idx = [i for i, p in enumerate(ag_pairs) if p[0]["key"] in iface_keys]
    interface_rmsd = None
    if iface_idx:
        ref_iface = ref_ag_coords[iface_idx]
        pred_iface = pred_ag_aligned[iface_idx]
        interface_rmsd = rmsd(ref_iface, pred_iface)

    return {
        "full_rmsd": full_rmsd,
        "interface_rmsd": interface_rmsd,
        "num_ref_atoms": int(ref_ag_coords.shape[0]),
        "num_pred_atoms": int(pred_ag_coords.shape[0]),
        "num_interface_atoms": int(len(iface_idx)),
        "rotation_matrix": r.tolist(),
        "translation": t.tolist(),
        "ab_pairing_stats": ab_stats,
        "ag_pairing_stats": ag_stats,
        "ref_detected_chains": sorted(ref_all),
        "pred_detected_chains": sorted(pred_all),
        "superposition_mode": superposition_mode,
        "fallback_note": fallback_note,
    }


def main():
    parser = argparse.ArgumentParser(description="Node 6 comparison with tolerant PDB parsing")
    parser.add_argument("--original_complex", required=True)
    parser.add_argument("--pred_pose", required=True)
    parser.add_argument("--ab_chains", required=True, help="Comma-separated antibody chains")
    parser.add_argument("--ag_chains", required=True, help="Comma-separated antigen chains")
    parser.add_argument("--output_dir", default="outputs")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    ab_chains = parse_chain_list(args.ab_chains)
    ag_chains = parse_chain_list(args.ag_chains)

    logger.info("Starting tolerant Node 6 comparison")
    logger.info("Reference complex: %s", args.original_complex)
    logger.info("Predicted pose: %s", args.pred_pose)
    logger.info("Antibody chains: %s", ab_chains)
    logger.info("Antigen chains: %s", ag_chains)

    if not os.path.exists(args.original_complex):
        logger.error("Reference complex not found: %s", args.original_complex)
        return 1
    if not os.path.exists(args.pred_pose):
        logger.error("Predicted pose not found: %s", args.pred_pose)
        return 1

    try:
        result = calc_rmsd_tolerant(args.original_complex, args.pred_pose, ab_chains, ag_chains)
        quality, interpretation = interpret_rmsd(result["full_rmsd"])

        pymol_script = generate_pymol_script(
            args.output_dir, args.original_complex, args.pred_pose, ab_chains, ag_chains
        )

        data = {
            "comparison_performed": True,
            "timestamp": datetime.now().isoformat(),
            "full_rmsd": result["full_rmsd"],
            "interface_rmsd": result["interface_rmsd"],
            "quality": quality,
            "interpretation": interpretation,
            "num_atoms_compared": result["num_pred_atoms"],
            "num_interface_atoms": result["num_interface_atoms"],
            "reference_complex": os.path.basename(args.original_complex),
            "predicted_pose": os.path.basename(args.pred_pose),
            "antibody_chains": ab_chains,
            "antigen_chains": ag_chains,
            "pymol_script": os.path.basename(pymol_script),
            "rotation_matrix": result["rotation_matrix"],
            "translation": result["translation"],
            "ab_pairing_stats": result["ab_pairing_stats"],
            "ag_pairing_stats": result["ag_pairing_stats"],
            "ref_detected_chains": result["ref_detected_chains"],
            "pred_detected_chains": result["pred_detected_chains"],
            "superposition_mode": result["superposition_mode"],
            "fallback_note": result["fallback_note"],
        }
    except Exception as e:
        logger.exception("Comparison failed")
        data = {
            "comparison_performed": False,
            "error": str(e),
            "timestamp": datetime.now().isoformat(),
        }

    with open(os.path.join(args.output_dir, "data.json"), "w") as f:
        json.dump(data, f, indent=2)
    with open(os.path.join(args.output_dir, "rmsd_analysis.json"), "w") as f:
        json.dump(data, f, indent=2)

    shutil.copy(args.pred_pose, os.path.join(args.output_dir, "best_match.pdb"))

    report_path = os.path.join(args.output_dir, "final_comparison_report.txt")
    with open(report_path, "w") as f:
        f.write("=" * 70 + "\n")
        f.write("ANTIBODY-ANTIGEN DOCKING COMPARISON REPORT\n")
        f.write("=" * 70 + "\n\n")
        if data.get("comparison_performed"):
            f.write(f"Timestamp: {data['timestamp']}\n")
            f.write("Status: SUCCESS\n\n")
            f.write("RMSD METRICS\n")
            f.write("-" * 70 + "\n")
            f.write(f"Full antigen RMSD:     {data['full_rmsd']:.3f} A\n")
            if data.get("interface_rmsd") is not None:
                f.write(f"Interface RMSD:        {data['interface_rmsd']:.3f} A\n")
            f.write(f"Atoms compared:        {data['num_atoms_compared']}\n")
            f.write(f"Interface atoms:       {data['num_interface_atoms']}\n\n")
            f.write("QUALITY ASSESSMENT\n")
            f.write("-" * 70 + "\n")
            f.write(f"Quality:      {data['quality'].upper()}\n")
            f.write(f"Interpretation: {data['interpretation']}\n\n")
            f.write("ALIGNMENT DETAILS\n")
            f.write("-" * 70 + "\n")
            f.write(f"Superposition mode: {data.get('superposition_mode', 'unknown')}\n")
            if data.get("fallback_note"):
                f.write(f"Note: {data['fallback_note']}\n")
            f.write("\n")
        else:
            f.write("Status: COMPARISON FAILED\n")
            f.write(f"Error: {data.get('error', 'Unknown error')}\n")

    print(json.dumps(data))
    return 0


if __name__ == "__main__":
    sys.exit(main())
