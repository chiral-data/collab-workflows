#!/usr/bin/env python3
"""Run LightDock docking simulation, post-process results, and collect top poses."""

import json
import os
import re
import shutil
import subprocess
import sys
import tarfile


def run_cmd(cmd, check=True):
    print(f"+ {' '.join(str(a) for a in cmd)}", flush=True)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout, flush=True)
    if result.stderr:
        print(result.stderr, file=sys.stderr, flush=True)
    if check and result.returncode != 0:
        print(f"ERROR: command exited with code {result.returncode}", flush=True)
        sys.exit(result.returncode)
    return result


def parse_gso_scores(gso_path):
    """Parse gso_N.out and return list of (glowworm_id, score) sorted by score descending."""
    poses = []
    with open(gso_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                glowworm_id = int(parts[0])
                score = float(parts[-1])
                poses.append((glowworm_id, score))
            except (ValueError, IndexError):
                pass
    poses.sort(key=lambda x: x[1], reverse=True)
    return poses


def main():
    steps = int(os.environ.get("PARAM_STEPS", "50"))
    scoring_function = os.environ.get("PARAM_SCORING_FUNCTION", "fastdfire")
    num_conformations = int(os.environ.get("PARAM_NUM_CONFORMATIONS", "200"))
    top_n = int(os.environ.get("PARAM_TOP_N", "10"))

    for fname in ["receptor.pdb", "ligand.pdb"]:
        src = os.path.join("inputs", fname)
        if not os.path.exists(src):
            print(f"ERROR: {src} not found", flush=True)
            sys.exit(1)
        shutil.copy(src, fname)

    workspace_tar = os.path.join("inputs", "lightdock_workspace.tar.gz")
    if not os.path.exists(workspace_tar):
        print(f"ERROR: {workspace_tar} not found", flush=True)
        sys.exit(1)
    print(f"Extracting {workspace_tar}...", flush=True)
    with tarfile.open(workspace_tar, "r:gz") as tar:
        tar.extractall(".")

    if not os.path.exists("setup.json"):
        print("ERROR: setup.json not found after extraction", flush=True)
        sys.exit(1)

    swarm_dirs = sorted(d for d in os.listdir(".") if d.startswith("swarm_") and os.path.isdir(d))
    num_swarms = len(swarm_dirs)
    print(f"Found {num_swarms} swarm directories", flush=True)

    print(f"\nRunning LightDock: {steps} steps, scoring={scoring_function}...", flush=True)
    run_cmd(["lightdock3.py", "setup.json", str(steps), "-s", scoring_function])

    gso_file = f"gso_{steps}.out"

    print("\nPost-processing swarms...", flush=True)
    for swarm_dir in swarm_dirs:
        swarm_gso = os.path.join(swarm_dir, gso_file)
        if not os.path.exists(swarm_gso):
            print(f"  WARNING: {swarm_gso} not found, skipping", flush=True)
            continue

        orig_dir = os.getcwd()
        os.chdir(swarm_dir)
        try:
            run_cmd(
                ["lgd_generate_conformations.py", "../receptor.pdb", "../ligand.pdb",
                 gso_file, str(num_conformations)],
                check=True,
            )
            run_cmd(["lgd_cluster_bsas.py", gso_file], check=True)
        finally:
            os.chdir(orig_dir)

    print("\nRanking poses...", flush=True)
    # lgd_rank.py takes num_swarms and steps (not conformations)
    run_cmd(["lgd_rank.py", str(num_swarms), str(steps)], check=True)

    top_poses = []

    # lightdock produces rank_by_scoring.list (sorted best→worst by scoring function)
    rank_file = "rank_by_scoring.list"
    if os.path.exists(rank_file):
        print(f"Parsing {rank_file}...", flush=True)
        with open(rank_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("Swarm"):
                    continue
                parts = line.split()
                try:
                    swarm_id = int(parts[0])
                    score = float(parts[-1])
                    # PDB name is like "lightdock_N.pdb"; build path relative to CWD
                    pdb_match = re.search(r'lightdock_\d+\.pdb', line)
                    pdb_file = f"swarm_{swarm_id}/{pdb_match.group()}" if pdb_match else None
                    top_poses.append({"score": score, "pdb_file": pdb_file, "rank": len(top_poses) + 1})
                except (ValueError, IndexError, AttributeError):
                    pass
        print(f"Parsed {len(top_poses)} ranked poses from {rank_file}", flush=True)

    if not top_poses:
        print("Falling back to per-swarm gso score extraction...", flush=True)
        swarm_bests = []
        for swarm_dir in swarm_dirs:
            swarm_gso = os.path.join(swarm_dir, gso_file)
            if not os.path.exists(swarm_gso):
                continue
            scored = parse_gso_scores(swarm_gso)
            if scored:
                best_id, best_score = scored[0]
                pdb_file = os.path.join(swarm_dir, f"lightdock_{best_id}.pdb")
                swarm_bests.append({"score": best_score, "pdb_file": pdb_file, "swarm": swarm_dir})
        swarm_bests.sort(key=lambda x: x["score"], reverse=True)
        for i, p in enumerate(swarm_bests):
            p["rank"] = i + 1
        top_poses = swarm_bests

    top_poses = top_poses[:top_n]

    print(f"\nCollecting top {len(top_poses)} poses into top_predictions.pdb...", flush=True)
    with open("top_predictions.pdb", "w") as out_f:
        for pose in top_poses:
            out_f.write(f"MODEL     {pose['rank']:4d}\n")
            out_f.write(f"REMARK  Rank {pose['rank']}  Score {pose['score']:.4f}\n")
            pdb_path = pose.get("pdb_file")
            if pdb_path and os.path.exists(pdb_path):
                with open(pdb_path) as pdb_f:
                    for line in pdb_f:
                        if not line.startswith(("END", "MODEL", "ENDMDL")):
                            out_f.write(line)
            out_f.write("ENDMDL\n")

    rank_data = {
        "num_swarms": num_swarms,
        "steps": steps,
        "scoring_function": scoring_function,
        "num_conformations": num_conformations,
        "top_n": top_n,
        "top_poses": [
            {"rank": p["rank"], "score": p["score"]}
            for p in top_poses
        ],
    }
    with open("rank_results.json", "w") as f:
        json.dump(rank_data, f, indent=2)

    print("Saved top_predictions.pdb and rank_results.json", flush=True)
    if top_poses:
        print(f"Best score: {top_poses[0]['score']:.4f}", flush=True)
    print("Docking complete.", flush=True)


if __name__ == "__main__":
    main()
