#!/usr/bin/env python3
"""Run lightdock3_setup.py to generate swarm positions and initial configurations."""

import os
import shutil
import subprocess
import sys
import tarfile


def run_cmd(cmd):
    print(f"+ {' '.join(cmd)}", flush=True)
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.stdout:
        print(result.stdout, flush=True)
    if result.stderr:
        print(result.stderr, file=sys.stderr, flush=True)
    if result.returncode != 0:
        print(f"ERROR: command exited with code {result.returncode}", flush=True)
        sys.exit(result.returncode)


def main():
    num_swarms = int(os.environ.get("PARAM_NUM_SWARMS", "25"))
    num_glowworms = int(os.environ.get("PARAM_NUM_GLOWWORMS", "200"))

    for fname in ["receptor.pdb", "ligand.pdb"]:
        src = os.path.join("inputs", fname)
        if not os.path.exists(src):
            print(f"ERROR: {src} not found", flush=True)
            sys.exit(1)
        shutil.copy(src, fname)
        print(f"Copied {src} -> {fname}", flush=True)

    run_cmd([
        "lightdock3_setup.py",
        "receptor.pdb", "ligand.pdb",
        "--swarms", str(num_swarms),
        "--glowworms", str(num_glowworms),
        "--noxt", "--noh", "--now",
    ])

    if not os.path.exists("setup.json"):
        print("ERROR: setup.json not created by lightdock3_setup.py", flush=True)
        sys.exit(1)

    swarm_dirs = sorted(d for d in os.listdir(".") if d.startswith("swarm_") and os.path.isdir(d))
    print(f"Setup complete: {len(swarm_dirs)} swarm directories created.", flush=True)

    archive_name = "lightdock_workspace.tar.gz"
    print(f"Archiving workspace to {archive_name}...", flush=True)
    with tarfile.open(archive_name, "w:gz") as tar:
        tar.add("setup.json")
        if os.path.isdir("init"):
            tar.add("init")
        for swarm_dir in swarm_dirs:
            tar.add(swarm_dir)
        # lightdock3.py reads these processed files from the working directory
        for fname in os.listdir("."):
            if fname.startswith("lightdock_") and not fname.endswith(".tar.gz"):
                tar.add(fname)

    archived_ld = [f for f in os.listdir(".") if f.startswith("lightdock_") and not f.endswith(".tar.gz")]
    print(f"Archived setup.json, init/, {len(swarm_dirs)} swarm dirs, and {len(archived_ld)} lightdock files.", flush=True)


if __name__ == "__main__":
    main()
