#!/usr/bin/env python3
"""
Node 04b: GNINA Screening

1. Parse pocket_config.txt for grid box
2. Split combined PDBQT into individual ligand files
3. Dock each ligand with GNINA CLI (CNN rescoring)
4. Collect all top poses → gnina_screening_poses.sdf (with CNN score properties)
5. Write scores → gnina_screening_scores.csv + gnina_docking_report.json
"""

import argparse
import csv
import json
import re
import subprocess
import sys
import time
from pathlib import Path


# ── Config ────────────────────────────────────────────────────────────────────

def parse_pocket_config(path):
    cfg = {}
    for line in Path(path).read_text().splitlines():
        m = re.match(r'(\w+)\s*=\s*([-\d.]+)', line.strip())
        if m:
            cfg[m.group(1)] = float(m.group(2))
    required = {"center_x", "center_y", "center_z", "size_x", "size_y", "size_z"}
    if missing := required - cfg.keys():
        print(f"ERROR: pocket_config.txt missing: {missing}", flush=True)
        sys.exit(1)
    return cfg


# ── PDBQT splitter (shared logic with 04a) ────────────────────────────────────

def split_pdbqt(combined_path, output_dir):
    """Split a combined multi-molecule PDBQT into individual files.

    Returns list of (name, Path) tuples.
    """
    content = Path(combined_path).read_text()
    blocks = []

    if re.search(r'^MODEL\s+\d+', content, re.MULTILINE):
        current = []
        for line in content.splitlines(keepends=True):
            if line.startswith("MODEL") and current:
                blocks.append("".join(current))
                current = []
            current.append(line)
        if current:
            blocks.append("".join(current))
    else:
        current = []
        for line in content.splitlines(keepends=True):
            current.append(line)
            if line.strip().startswith("TORSDOF"):
                blocks.append("".join(current))
                current = []
        if current and any(l.startswith(("ATOM", "HETATM")) for l in current):
            blocks.append("".join(current))

    files = []
    for i, block in enumerate(blocks):
        if not block.strip():
            continue
        m = re.search(r'REMARK\s+Name\s*=\s*(.+)', block)
        raw_name = m.group(1).strip() if m else f"mol_{i+1}"
        safe_name = re.sub(r'[^\w.-]', '_', raw_name)
        out_path = output_dir / f"lig_{i+1:04d}_{safe_name}.pdbqt"
        out_path.write_text(block)
        files.append((raw_name, out_path))

    return files


# ── GNINA docking ─────────────────────────────────────────────────────────────

def dock_with_gnina(receptor, ligand_path, cfg, out_sdf, exhaustiveness, num_modes):
    """Run GNINA CLI for one ligand. Returns (success, log_text)."""
    cmd = [
        "gnina",
        "-r",                  str(receptor),
        "-l",                  str(ligand_path),
        "--center_x",          str(cfg["center_x"]),
        "--center_y",          str(cfg["center_y"]),
        "--center_z",          str(cfg["center_z"]),
        "--size_x",            str(cfg["size_x"]),
        "--size_y",            str(cfg["size_y"]),
        "--size_z",            str(cfg["size_z"]),
        "--exhaustiveness",    str(exhaustiveness),
        "--num_modes",         str(num_modes),
        "--cnn_scoring",       "rescore",
        "--no_gpu",
        "-o",                  str(out_sdf),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    return res.returncode == 0, res.stdout + res.stderr


def parse_gnina_top_scores(sdf_path):
    """Return (cnn_score, cnn_affinity, vina_affinity) for the top pose, or (None,None,None)."""
    try:
        from rdkit import Chem
        suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False, sanitize=False)
        for mol in suppl:
            if mol is None:
                continue
            cnn   = float(mol.GetProp("CNNscore"))    if mol.HasProp("CNNscore")          else None
            cnn_a = float(mol.GetProp("CNNaffinity")) if mol.HasProp("CNNaffinity")        else None
            vina  = float(mol.GetProp("minimizedAffinity")) if mol.HasProp("minimizedAffinity") else None
            return cnn, cnn_a, vina
    except Exception:
        pass
    return None, None, None


def count_sdf_mols(sdf_path):
    try:
        return Path(sdf_path).read_text().count("$$$$")
    except Exception:
        return 0


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Node 04b: GNINA Screening")
    parser.add_argument("--receptor",       default="receptor.pdbqt")
    parser.add_argument("--library",        default="optimized_screening_library.pdbqt")
    parser.add_argument("--pocket-config",  default="pocket_config.txt")
    parser.add_argument("--exhaustiveness", type=int, default=8)
    parser.add_argument("--num-modes",      type=int, default=9)
    args = parser.parse_args()

    for p in (args.receptor, args.library, args.pocket_config):
        if not Path(p).exists():
            print(f"ERROR: input not found: {p}", flush=True)
            sys.exit(1)

    cfg = parse_pocket_config(args.pocket_config)
    print(f"Grid center: ({cfg['center_x']:.2f}, {cfg['center_y']:.2f}, {cfg['center_z']:.2f})", flush=True)
    print(f"Box size:    ({cfg['size_x']:.2f}, {cfg['size_y']:.2f}, {cfg['size_z']:.2f})", flush=True)

    # Split library
    print("\n--- Splitting screening library ---", flush=True)
    lig_dir = Path("ligands_split")
    lig_dir.mkdir(exist_ok=True)
    ligands = split_pdbqt(args.library, lig_dir)
    print(f"  {len(ligands)} ligands to dock", flush=True)
    if not ligands:
        print("ERROR: could not split PDBQT library", flush=True)
        sys.exit(1)

    # Dock each ligand
    print("\n--- Docking ---", flush=True)
    out_dir = Path("gnina_poses")
    out_dir.mkdir(exist_ok=True)
    log_lines = []
    results = []

    t_dock_start = time.monotonic()
    for name, lig_path in ligands:
        out_sdf = out_dir / f"{lig_path.stem}_out.sdf"
        ok, log = dock_with_gnina(
            args.receptor, lig_path, cfg, out_sdf,
            args.exhaustiveness, args.num_modes
        )
        log_lines.append(f"=== {name} ===\n{log}")

        if ok and out_sdf.exists():
            cnn, cnn_a, vina_a = parse_gnina_top_scores(out_sdf)
            n_poses = count_sdf_mols(out_sdf)
            cnn_str  = f"{cnn:.3f}"    if cnn    is not None else "N/A"
            cnna_str = f"{cnn_a:.2f}"  if cnn_a  is not None else "N/A"
            vina_str = f"{vina_a:.2f}" if vina_a is not None else "N/A"
            print(f"  {name}: CNNscore={cnn_str}  CNN_affinity={cnna_str}  Vina_affinity={vina_str}  ({n_poses} poses)",
                  flush=True)
            results.append({
                "name":          name,
                "status":        "success",
                "cnn_score":     cnn,
                "cnn_affinity":  cnn_a,
                "vina_affinity": vina_a,
                "num_poses":     n_poses,
                "sdf_file":      str(out_sdf),
            })
        else:
            print(f"  {name}: FAILED", flush=True)
            results.append({"name": name, "status": "failed",
                             "cnn_score": None, "cnn_affinity": None, "vina_affinity": None})
    t_dock_end = time.monotonic()
    wall_clock_s = round(t_dock_end - t_dock_start, 1)
    print(f"\n  Docking wall-clock: {wall_clock_s} s", flush=True)

    # Write combined log
    Path("gnina_runtime_log.txt").write_text("\n\n".join(log_lines))

    # Rank by CNN score (higher = better for GNINA)
    successful = [r for r in results if r["status"] == "success" and r["cnn_score"] is not None]
    successful.sort(key=lambda r: r["cnn_score"], reverse=True)
    for rank, r in enumerate(successful, 1):
        r["rank"] = rank

    # Collect top poses → combined SDF
    print("\n--- Collecting top poses ---", flush=True)
    combined_sdf_blocks = []
    for r in results:
        if r["status"] != "success":
            continue
        sdf_text = Path(r["sdf_file"]).read_text()
        # Take only the first molecule (top pose) from the SDF
        first_mol_end = sdf_text.find("$$$$")
        if first_mol_end >= 0:
            block = sdf_text[: first_mol_end + 4]
        else:
            block = sdf_text
        combined_sdf_blocks.append(block)

    Path("gnina_screening_poses.sdf").write_text("\n".join(combined_sdf_blocks))
    print(f"  Wrote gnina_screening_poses.sdf ({len(combined_sdf_blocks)} top poses)", flush=True)

    # Scores CSV
    with open("gnina_screening_scores.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["rank", "name", "cnn_score",
                                                "cnn_affinity_kcal_mol",
                                                "vina_affinity_kcal_mol",
                                                "num_poses", "status"])
        writer.writeheader()
        for r in results:
            writer.writerow({
                "rank":                    r.get("rank", ""),
                "name":                    r["name"],
                "cnn_score":               r.get("cnn_score", ""),
                "cnn_affinity_kcal_mol":   r.get("cnn_affinity", ""),
                "vina_affinity_kcal_mol":  r.get("vina_affinity", ""),
                "num_poses":               r.get("num_poses", ""),
                "status":                  r["status"],
            })
    print("  Wrote gnina_screening_scores.csv", flush=True)

    # Report JSON
    n_total = len(results)
    report = {
        "tool":           "GNINA",
        "scoring":        "CNN rescore",
        "total_ligands":  n_total,
        "successful":     len(successful),
        "failed":         n_total - len(successful),
        "top_hits": [
            {"rank": r["rank"], "name": r["name"],
             "cnn_score": r["cnn_score"],
             "cnn_affinity_kcal_mol": r["cnn_affinity"]}
            for r in successful[:5]
        ],
        "exhaustiveness":               args.exhaustiveness,
        "num_modes":                    args.num_modes,
        "wall_clock_seconds":           wall_clock_s,
        "wall_clock_per_ligand_seconds": round(wall_clock_s / n_total, 2) if n_total else None,
    }
    with open("gnina_docking_report.json", "w") as f:
        json.dump(report, f, indent=2)

    print("\nTop 5 by GNINA CNNscore:", flush=True)
    for r in successful[:5]:
        print(f"  #{r['rank']:2d}  {r['name']:<40s}  CNNscore={r['cnn_score']:.3f}  "
              f"CNN_aff={r['cnn_affinity']:.2f} kcal/mol", flush=True)

    print("\nNode 04b completed", flush=True)


if __name__ == "__main__":
    main()
