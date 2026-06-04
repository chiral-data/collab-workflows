#!/usr/bin/env python3
"""
Node 04a: AutoDock Vina Screening

1. Parse pocket_config.txt for grid box
2. Split combined PDBQT into individual ligand files
3. Dock each ligand with AutoDock Vina CLI
4. Collect all top poses → vina_screening_poses.pdbqt
5. Write scores → vina_screening_scores.csv + vina_docking_report.json
"""

import argparse
import csv
import json
import re
import subprocess
import sys
import tempfile
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


# ── PDBQT splitter ────────────────────────────────────────────────────────────

def split_pdbqt(combined_path, output_dir):
    """Split a combined multi-molecule PDBQT into individual files.

    Handles both Meeko format (ROOT…TORSDOF) and obabel format (MODEL…ENDMDL).
    Returns list of (name, Path) tuples.
    """
    content = Path(combined_path).read_text()
    blocks = []

    if re.search(r'^MODEL\s+\d+', content, re.MULTILINE):
        # obabel MODEL/ENDMDL format
        current = []
        for line in content.splitlines(keepends=True):
            if line.startswith("MODEL") and current:
                blocks.append("".join(current))
                current = []
            current.append(line)
        if current:
            blocks.append("".join(current))
    else:
        # Meeko: each molecule ends with TORSDOF
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


# ── Vina docking ──────────────────────────────────────────────────────────────

def dock_with_vina(receptor, ligand_path, cfg, out_pdbqt, exhaustiveness, num_modes):
    """Run Vina CLI for one ligand. Returns True on success."""
    cmd = [
        "vina",
        "--receptor",      str(receptor),
        "--ligand",        str(ligand_path),
        "--center_x",      str(cfg["center_x"]),
        "--center_y",      str(cfg["center_y"]),
        "--center_z",      str(cfg["center_z"]),
        "--size_x",        str(cfg["size_x"]),
        "--size_y",        str(cfg["size_y"]),
        "--size_z",        str(cfg["size_z"]),
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes",      str(num_modes),
        "--out",           str(out_pdbqt),
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    return res.returncode == 0, res.stdout + res.stderr


def parse_vina_scores(pdbqt_path):
    """Return list of (pose_idx, affinity) from a Vina output PDBQT."""
    scores = []
    idx = 0
    for line in Path(pdbqt_path).read_text().splitlines():
        if "REMARK VINA RESULT:" in line:
            parts = line.split()
            try:
                scores.append((idx, float(parts[3])))
            except (IndexError, ValueError):
                pass
            idx += 1
    return scores


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Node 04a: AutoDock Vina Screening")
    parser.add_argument("--receptor",        default="receptor.pdbqt")
    parser.add_argument("--library",         default="optimized_screening_library.pdbqt")
    parser.add_argument("--pocket-config",   default="pocket_config.txt")
    parser.add_argument("--exhaustiveness",  type=int, default=8)
    parser.add_argument("--num-modes",       type=int, default=9)
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
    out_dir = Path("vina_poses")
    out_dir.mkdir(exist_ok=True)
    log_lines = []
    results = []

    for name, lig_path in ligands:
        out_pdbqt = out_dir / f"{lig_path.stem}_out.pdbqt"
        ok, log = dock_with_vina(
            args.receptor, lig_path, cfg, out_pdbqt,
            args.exhaustiveness, args.num_modes
        )
        log_lines.append(f"=== {name} ===\n{log}")

        if ok and out_pdbqt.exists():
            scores = parse_vina_scores(out_pdbqt)
            best = scores[0][1] if scores else None
            print(f"  {name}: best ΔG = {best:.2f} kcal/mol  ({len(scores)} poses)", flush=True)
            results.append({
                "name":          name,
                "status":        "success",
                "best_affinity": best,
                "num_poses":     len(scores),
                "all_affinities": [s[1] for s in scores],
                "pdbqt_file":    str(out_pdbqt),
            })
        else:
            print(f"  {name}: FAILED", flush=True)
            results.append({"name": name, "status": "failed", "best_affinity": None})

    # Write combined log
    Path("vina_runtime_log.txt").write_text("\n\n".join(log_lines))

    # Rank by best affinity (most negative = best)
    successful = [r for r in results if r["status"] == "success" and r["best_affinity"] is not None]
    successful.sort(key=lambda r: r["best_affinity"])
    for rank, r in enumerate(successful, 1):
        r["rank"] = rank

    # Collect top poses → combined PDBQT
    print("\n--- Collecting top poses ---", flush=True)
    combined_blocks = []
    for r in results:
        if r["status"] != "success":
            continue
        pdbqt_text = Path(r["pdbqt_file"]).read_text()
        # Extract first MODEL block (top pose) or the whole file if no MODEL markers
        if "MODEL" in pdbqt_text:
            m = re.search(r'(MODEL\s+1.*?ENDMDL)', pdbqt_text, re.DOTALL)
            block = m.group(1) if m else pdbqt_text
        else:
            block = pdbqt_text
        combined_blocks.append(f"REMARK Name = {r['name']}\n{block}")
    Path("vina_screening_poses.pdbqt").write_text("\n".join(combined_blocks))
    print(f"  Wrote vina_screening_poses.pdbqt ({len(combined_blocks)} top poses)", flush=True)

    # Scores CSV
    with open("vina_screening_scores.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["rank", "name", "best_affinity_kcal_mol",
                                                "num_poses", "status"])
        writer.writeheader()
        for r in results:
            writer.writerow({
                "rank":                  r.get("rank", ""),
                "name":                  r["name"],
                "best_affinity_kcal_mol": r.get("best_affinity", ""),
                "num_poses":             r.get("num_poses", ""),
                "status":                r["status"],
            })
    print("  Wrote vina_screening_scores.csv", flush=True)

    # Report JSON
    report = {
        "tool":             "AutoDock Vina",
        "total_ligands":    len(results),
        "successful":       len(successful),
        "failed":           len(results) - len(successful),
        "top_hits": [
            {"rank": r["rank"], "name": r["name"],
             "best_affinity_kcal_mol": r["best_affinity"]}
            for r in successful[:5]
        ],
        "exhaustiveness": args.exhaustiveness,
        "num_modes":      args.num_modes,
    }
    with open("vina_docking_report.json", "w") as f:
        json.dump(report, f, indent=2)

    print("\nTop 5 by Vina ΔG:", flush=True)
    for r in successful[:5]:
        print(f"  #{r['rank']:2d}  {r['name']:<40s}  {r['best_affinity']:.2f} kcal/mol", flush=True)

    print("\nNode 04a completed", flush=True)


if __name__ == "__main__":
    main()
