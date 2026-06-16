#!/usr/bin/env python3
"""
Node 03: Target Redocking QC

1. Extract the native ligand's crystal pose from receptor.pdb → crystal_ligand.pdbqt
2. Parse pocket_config.txt for grid box parameters
3. Run GNINA: self-dock native_ligand.pdbqt with --crystal_pose for RMSD
4. Extract top pose's CNNscore and RMSD from GNINA output
5. Pass/fail QC; exit 1 on failure (stops the Silva workflow)
"""

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

EXCLUDE_RESIDUES = {
    "HOH", "WAT", "H2O",
    "SO4", "PO4", "NO3",
    "CL", "NA", "MG", "ZN", "CA", "K", "FE", "MN", "CU", "CO",
    "GOL", "EDO", "PEG", "MPD", "IPA", "DMS", "ACT", "ACE", "MSE",
    "UNX", "UNL",
}


# ── Pocket config ─────────────────────────────────────────────────────────────

def parse_pocket_config(path):
    """Return dict with center_x/y/z and size_x/y/z from pocket_config.txt."""
    cfg = {}
    for line in Path(path).read_text().splitlines():
        m = re.match(r'(\w+)\s*=\s*([-\d.]+)', line.strip())
        if m:
            cfg[m.group(1)] = float(m.group(2))
    required = {"center_x", "center_y", "center_z", "size_x", "size_y", "size_z"}
    missing = required - cfg.keys()
    if missing:
        print(f"ERROR: pocket_config.txt missing fields: {missing}", flush=True)
        sys.exit(1)
    return cfg


# ── Crystal reference extraction ──────────────────────────────────────────────

def extract_crystal_reference_pdb(receptor_pdb, out_pdb):
    """Write the native (non-solvent) HETATM residue from receptor.pdb to out_pdb."""
    # Count atoms per HETATM residue to identify the ligand
    counts = {}
    for line in Path(receptor_pdb).read_text().splitlines():
        if not line.startswith("HETATM"):
            continue
        res = line[17:20].strip()
        if res in EXCLUDE_RESIDUES:
            continue
        key = (res, line[21].strip(), line[22:26].strip())
        counts[key] = counts.get(key, 0) + 1

    if not counts:
        print("ERROR: no non-solvent HETATM residue found in receptor.pdb", flush=True)
        sys.exit(1)

    best = max(counts, key=lambda k: counts[k])
    lig_res, lig_chain, lig_seq = best
    print(f"  Crystal ligand: {lig_res} chain {lig_chain} res {lig_seq} "
          f"({counts[best]} atoms)", flush=True)

    # Write those atoms as a minimal PDB
    lines_out = []
    for line in Path(receptor_pdb).read_text().splitlines():
        if not line.startswith("HETATM"):
            continue
        if (line[17:20].strip() == lig_res and
                line[21].strip() == lig_chain and
                line[22:26].strip() == lig_seq):
            lines_out.append(line)
    lines_out.append("END")
    Path(out_pdb).write_text("\n".join(lines_out) + "\n")
    print(f"  Wrote crystal reference: {out_pdb} ({len(lines_out)-1} atoms)", flush=True)
    return lig_res


def pdb_to_pdbqt(pdb_path, pdbqt_path):
    """Convert PDB to PDBQT with OpenBabel."""
    res = subprocess.run(
        ["obabel", str(pdb_path), "-O", str(pdbqt_path), "-xh"],
        capture_output=True, text=True
    )
    if res.returncode != 0 or not Path(pdbqt_path).exists():
        print(f"  WARNING: obabel failed to convert {pdb_path} to PDBQT", flush=True)
        print(res.stderr, flush=True)
        return False
    return True


# ── GNINA docking ─────────────────────────────────────────────────────────────

def run_gnina(receptor_pdbqt, ligand_pdbqt, crystal_pdbqt, cfg, out_sdf, log_file,
              cnn_scoring="all"):
    """Run GNINA redocking; return True on success."""
    cmd = [
        "gnina",
        "-r", str(receptor_pdbqt),
        "-l", str(ligand_pdbqt),
        "--center_x", str(cfg["center_x"]),
        "--center_y", str(cfg["center_y"]),
        "--center_z", str(cfg["center_z"]),
        "--size_x",   str(cfg["size_x"]),
        "--size_y",   str(cfg["size_y"]),
        "--size_z",   str(cfg["size_z"]),
        "--num_modes",      "15",
        "--exhaustiveness", "32",
        "--cnn_scoring",    cnn_scoring,
        "--no_gpu",
        "-o", str(out_sdf),
        "--log", str(log_file),
    ]
    print(f"  Running: {' '.join(cmd)}", flush=True)
    res = subprocess.run(cmd, capture_output=True, text=True)
    # GNINA exits 0 on success; pipe its output live
    print(res.stdout, flush=True)
    if res.returncode != 0:
        print("GNINA stderr:", res.stderr, flush=True)
        return False
    return True


# ── Parse GNINA output ────────────────────────────────────────────────────────

def parse_gnina_scores(sdf_path, log_path):
    """
    Return (cnn_score, rmsd) for the top-ranked pose.

    CNNscore comes from SDF properties or the log score table (column 3).
    RMSD is NOT parsed from the log: without --crystal_pose the 4th column
    is CNN affinity, not RMSD.  The caller must compute RMSD via RDKit.
    """
    cnn_score, rmsd = None, None

    # Try SDF properties first
    try:
        from rdkit import Chem
        suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
        for mol in suppl:
            if mol is None:
                continue
            if mol.HasProp("CNNscore"):
                cnn_score = float(mol.GetProp("CNNscore"))
            if mol.HasProp("RMSD"):
                rmsd = float(mol.GetProp("RMSD"))
            break  # top pose only
    except Exception as e:
        print(f"  WARNING: could not parse SDF properties: {e}", flush=True)

    # Fall back to log score table for CNNscore only (column 3 = CNN pose score)
    if cnn_score is None:
        log_text = Path(log_path).read_text() if Path(log_path).exists() else ""
        for line in log_text.splitlines():
            # "   1   -3.55   0.87   0.3807   2.923"
            # cols: mode  affinity  intramol  CNNscore  CNNaffinity
            m = re.match(r'\s*1\s+[-\d.]+\s+[\d.]+\s+([\d.]+)', line)
            if m:
                cnn_score = float(m.group(1))
                break

    return cnn_score, rmsd


def compute_rmsd_rdkit(docked_sdf, crystal_pdbqt):
    """
    Heavy-atom RMSD between the top docked pose and the crystal pose.

    Reads reference coords from crystal_ligand.pdbqt (not the PDB) because
    OpenBabel reorders atoms during PDB→PDBQT conversion and GNINA's SDF output
    preserves the PDBQT atom ordering.  Direct positional RMSD is correct here
    since both structures share that ordering.
    """
    try:
        import numpy as np
        from rdkit import Chem

        # Docked top pose (heavy atoms only)
        suppl = Chem.SDMolSupplier(str(docked_sdf), removeHs=True, sanitize=False)
        docked = next((m for m in suppl if m is not None), None)
        if docked is None:
            return None
        docked_pos = np.array(docked.GetConformer().GetPositions())

        # Crystal coords from PDBQT (same atom ordering as the docked SDF)
        crystal_coords = []
        for line in Path(crystal_pdbqt).read_text().splitlines():
            if line.startswith(("ATOM", "HETATM")):
                try:
                    crystal_coords.append([
                        float(line[30:38]), float(line[38:46]), float(line[46:54])
                    ])
                except ValueError:
                    pass
        if not crystal_coords:
            return None
        ref_pos = np.array(crystal_coords)

        n = min(len(docked_pos), len(ref_pos))
        rmsd = float(np.sqrt(np.mean(np.sum((docked_pos[:n] - ref_pos[:n])**2, axis=1))))
        print(f"  Fallback RMSD (heavy-atom, PDBQT-ordered): {rmsd:.3f} A "
              f"({len(docked_pos)} docked, {len(ref_pos)} crystal atoms)", flush=True)
        return rmsd
    except Exception as e:
        print(f"  WARNING: fallback RMSD failed: {e}", flush=True)
        return None


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Node 03: Target Redocking QC")
    parser.add_argument("--receptor-pdbqt",      default="receptor.pdbqt")
    parser.add_argument("--native-ligand-pdbqt", default="native_ligand.pdbqt")
    parser.add_argument("--receptor-pdb",        default="receptor.pdb")
    parser.add_argument("--pocket-config",       default="pocket_config.txt")
    parser.add_argument("--cnn-score-threshold", type=float, default=0.90)
    parser.add_argument("--rmsd-threshold",      type=float, default=2.0)
    parser.add_argument("--cnn-scoring",         default="rescore",
                        choices=["all", "rescore", "none"])
    args = parser.parse_args()

    cnn_thr  = args.cnn_score_threshold
    rmsd_thr = args.rmsd_threshold
    cnn_scoring = args.cnn_scoring

    print(f"CNN score threshold : >= {cnn_thr}", flush=True)
    print(f"RMSD threshold      : <= {rmsd_thr} A", flush=True)
    print(f"CNN scoring mode    : {cnn_scoring}", flush=True)

    # 1. Parse pocket box
    print("\n--- Step 1: Parse pocket grid box ---", flush=True)
    cfg = parse_pocket_config(args.pocket_config)
    print(f"  Center: ({cfg['center_x']:.2f}, {cfg['center_y']:.2f}, {cfg['center_z']:.2f})", flush=True)
    print(f"  Size:   ({cfg['size_x']:.2f}, {cfg['size_y']:.2f}, {cfg['size_z']:.2f})", flush=True)

    # 2. Extract crystal reference
    print("\n--- Step 2: Extract crystal reference pose ---", flush=True)
    crystal_pdb   = Path("crystal_ligand.pdb")
    crystal_pdbqt = Path("crystal_ligand.pdbqt")
    lig_id = extract_crystal_reference_pdb(args.receptor_pdb, crystal_pdb)
    has_crystal_pdbqt = pdb_to_pdbqt(crystal_pdb, crystal_pdbqt)

    # 3. Run GNINA
    print("\n--- Step 3: Run GNINA redocking ---", flush=True)
    out_sdf  = Path("redocked_native.sdf")
    log_file = Path("gnina_redock_log.txt")

    gnina_ok = False
    # Prefer crystal_pdbqt as ligand: it has the correct bound-state geometry
    # extracted from the PDB, vs native_ligand_pdbqt which may come from a 2D SDF.
    ligand_for_docking = str(crystal_pdbqt) if has_crystal_pdbqt else args.native_ligand_pdbqt
    if has_crystal_pdbqt:
        gnina_ok = run_gnina(
            args.receptor_pdbqt, ligand_for_docking,
            crystal_pdbqt, cfg, out_sdf, log_file,
            cnn_scoring=cnn_scoring
        )
    else:
        cmd_no_crystal = [
            "gnina",
            "-r", args.receptor_pdbqt,
            "-l", args.native_ligand_pdbqt,
            "--center_x", str(cfg["center_x"]),
            "--center_y", str(cfg["center_y"]),
            "--center_z", str(cfg["center_z"]),
            "--size_x",   str(cfg["size_x"]),
            "--size_y",   str(cfg["size_y"]),
            "--size_z",   str(cfg["size_z"]),
            "--num_modes",      "15",
            "--exhaustiveness", "32",
            "--cnn_scoring",    cnn_scoring,
            "--no_gpu",
            "-o", str(out_sdf),
            "--log", str(log_file),
        ]
        res = subprocess.run(cmd_no_crystal, capture_output=True, text=True)
        print(res.stdout, flush=True)
        gnina_ok = (res.returncode == 0)

    if not gnina_ok or not out_sdf.exists():
        print("ERROR: GNINA docking failed", flush=True)
        sys.exit(1)

    # 4. Extract scores
    print("\n--- Step 4: Extract scores ---", flush=True)
    cnn_score, rmsd = parse_gnina_scores(out_sdf, log_file)

    # Fallback RMSD if GNINA didn't report it
    if rmsd is None:
        print("  GNINA did not report RMSD; computing via RDKit fallback ...", flush=True)
        rmsd = compute_rmsd_rdkit(out_sdf, crystal_pdbqt)

    print(f"  CNNscore (top pose): {cnn_score}", flush=True)
    print(f"  RMSD vs crystal:     {rmsd} A", flush=True)

    # 5. QC decision
    print("\n--- Step 5: QC decision ---", flush=True)
    cnn_pass  = (cnn_score  is not None) and (cnn_score  >= cnn_thr)
    rmsd_pass = (rmsd       is not None) and (rmsd       <= rmsd_thr)
    qc_pass   = cnn_pass and rmsd_pass

    cnn_str  = f"{cnn_score:.3f}" if cnn_score is not None else "None"
    rmsd_str = f"{rmsd:.3f}"     if rmsd      is not None else "None"
    print(f"  CNN score  {cnn_str} >= {cnn_thr}: {'PASS' if cnn_pass  else 'FAIL'}", flush=True)
    print(f"  RMSD       {rmsd_str} A <= {rmsd_thr} A: {'PASS' if rmsd_pass else 'FAIL'}", flush=True)
    print(f"  Overall QC: {'PASS' if qc_pass else 'FAIL'}", flush=True)

    result = {
        "native_ligand_id":   lig_id,
        "top_pose_cnn_score": cnn_score,
        "top_pose_rmsd_A":    rmsd,
        "cnn_score_threshold":  cnn_thr,
        "rmsd_threshold_A":     rmsd_thr,
        "cnn_pass":  cnn_pass,
        "rmsd_pass": rmsd_pass,
        "qc_pass":   qc_pass,
    }
    with open("qc_validation_results.json", "w") as f:
        json.dump(result, f, indent=2)
    print("\nWrote qc_validation_results.json", flush=True)
    print(json.dumps(result, indent=2), flush=True)

    if not qc_pass:
        print("\nQC FAILED — stopping workflow. "
              "Check pocket_config.txt and native ligand preparation.", flush=True)
        sys.exit(1)

    print("\nNode 03 completed", flush=True)


if __name__ == "__main__":
    main()
