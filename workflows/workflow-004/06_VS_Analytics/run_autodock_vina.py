#!/usr/bin/env python3
"""
Node 6: Virtual Screening
Performs molecular docking using AutoDock Vina
"""

import json
import subprocess
import shutil
from pathlib import Path
from datetime import datetime
import os
import sys

def run_vina():
    workdir = Path("./")
    outputs_dir = workdir / "outputs"
    outputs_dir.mkdir(exist_ok=True)
    
    print("[VS Analytics] ========== Identifying Components ==========")
    print(f"[VS Analytics] Working directory: {workdir.resolve()}")
    
    # Flatten inputs
    for ext in ["*.pdbqt", "*.json"]:
        for file in workdir.glob(f"**/{ext}"):
            if file.parent != workdir:
                target = workdir / file.name
                if not target.exists():
                    shutil.copy(file, target)
                    print(f"[VS Analytics] Flattened input: {file} → {target}")

    # Identify Receptor and Ligands
    all_pdbqts = list(workdir.glob("*.pdbqt"))
    if not all_pdbqts:
        print(f"[VS Analytics] ERROR: No PDBQT files found")
        sys.exit(1)
        
    # Heuristic: Largest file is receptor
    all_pdbqts.sort(key=lambda x: x.stat().st_size, reverse=True)
    receptor = all_pdbqts[0]
    ligands = all_pdbqts[1:]
    
    print(f"[VS Analytics] Receptor: {receptor.name}")
    print(f"[VS Analytics] Ligands found: {len(ligands)}")

    if not ligands:
        print("[VS Analytics] ERROR: No ligands found")
        sys.exit(1)

    # Load Grid Configuration
    print("[VS Analytics] ========== Loading Grid Configuration ==========")
    grid_config = None
    possible_configs = [
        workdir / "pocket_discovery_metadata.json",
        workdir / "outputs" / "pocket_discovery_metadata.json",
        workdir / "grid_config.json",
        workdir / "outputs" / "grid_config.json"
    ]
    
    for p in possible_configs:
        if p.exists():
            try:
                with open(p) as f:
                    data = json.load(f)
                    grid_config = data.get("grid_config", data)
                    if "center_x" in grid_config:
                        print(f"[VS Analytics] Loaded grid config from {p}")
                        break
            except:
                continue

    if not grid_config or "center_x" not in grid_config:
        print("[VS Analytics] ERROR: Valid grid configuration not found")
        sys.exit(1)

    print(f"[VS Analytics] Grid Center: ({grid_config['center_x']:.2f}, {grid_config['center_y']:.2f}, {grid_config['center_z']:.2f})")

    # Setup Vina Parameters
    params = {
        "exhaustiveness": int(os.environ.get("PARAM_EXHAUSTIVENESS", "32")),
        "num_modes": int(os.environ.get("PARAM_NUM_MODES", "10")),
        "energy_range": float(os.environ.get("PARAM_ENERGY_RANGE", "5.0"))
    }

    # GPU Check
    use_gpu = False
    try:
        check_gpu = subprocess.run(["vina", "--help"], capture_output=True, text=True)
        if "--gpu" in check_gpu.stdout:
            use_gpu = True
            print("[VS Analytics] GPU acceleration enabled")
    except:
        pass

    # Docking Loop
    print("[VS Analytics] ========== Starting Docking ==========")
    docking_results = []
    
    def split_poses(pdbqt_file, ligand_name):
        """Split multi-model PDBQT into individual pose files"""
        poses_data = []
        if not pdbqt_file.exists():
            return poses_data
            
        with open(pdbqt_file, "r") as f:
            content = f.read()
        
        models = content.split("MODEL")[1:]
        
        for i, model_str in enumerate(models):
            pose_num = i + 1
            pose_text = "MODEL" + model_str
            
            affinity = None
            for line in pose_text.splitlines():
                if "REMARK VINA RESULT:" in line:
                    parts = line.split()
                    if len(parts) >= 4:
                        affinity = float(parts[3])
                    break
            
            pose_filename = f"{ligand_name}_pose{pose_num}.pdbqt"
            pose_path = outputs_dir / pose_filename
            with open(pose_path, "w") as f:
                f.write(pose_text)
            
            poses_data.append({
                "pose": pose_num,
                "affinity": affinity,
                "file": pose_filename
            })
            
        return poses_data

    for ligand in ligands:
        lig_name = ligand.stem
        out_pdbqt = outputs_dir / f"{lig_name}_docked.pdbqt"
        log_file = outputs_dir / f"{lig_name}_docking.log"
        
        cmd = [
            "vina",
            "--receptor", str(receptor.resolve()),
            "--ligand", str(ligand.resolve()),
            "--center_x", str(grid_config["center_x"]),
            "--center_y", str(grid_config["center_y"]),
            "--center_z", str(grid_config["center_z"]),
            "--size_x", str(grid_config.get("size_x", 20.0)),
            "--size_y", str(grid_config.get("size_y", 20.0)),
            "--size_z", str(grid_config.get("size_z", 20.0)),
            "--exhaustiveness", str(params["exhaustiveness"]),
            "--num_modes", str(params["num_modes"]),
            "--energy_range", str(params["energy_range"]),
            "--out", str(out_pdbqt.resolve())
        ]
        
        if use_gpu:
            cmd.append("--gpu")

        print(f" → Docking {lig_name}...", end=" ", flush=True)
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            with open(log_file, "w") as f:
                f.write(result.stdout)
            
            poses = split_poses(out_pdbqt, lig_name)
            best_affinity = poses[0]["affinity"] if poses else 0.0
            
            docking_results.append({
                "ligand": lig_name,
                "best_affinity": best_affinity,
                "poses": poses,
                "output_pdbqt": out_pdbqt.name,
                "status": "success"
            })
            print(f"Done. Best: {best_affinity:.2f} kcal/mol")
            
        except subprocess.CalledProcessError as e:
            print(f"Failed")
            with open(log_file, "w") as f:
                f.write(e.stdout + "\n" + e.stderr)
            docking_results.append({
                "ligand": lig_name, "status": "failed", "error": str(e)
            })

    print("[VS Analytics] ========== Ranking Results ==========")
    
    # Ranking
    successful = [r for r in docking_results if r.get("status") == "success"]
    successful.sort(key=lambda x: x.get("best_affinity", 0))
    
    for i, res in enumerate(successful):
        res["rank"] = i + 1
        print(f" → Rank #{i+1}: {res['ligand']} ({res['best_affinity']:.2f} kcal/mol)")

    # Finalize Assets
    shutil.copy(receptor, outputs_dir / "receptor.pdbqt")
    
    final_data = {
        "summary": {
            "total_compounds": len(docking_results),
            "successful_dockings": len(successful),
            "failed_dockings": len([r for r in docking_results if r.get("status") == "failed"]),
            "top_hit": successful[0]["ligand"] if successful else None,
            "best_affinity": successful[0]["best_affinity"] if successful else None
        },
        "receptor_file": "receptor.pdbqt",
        "docking_results": docking_results,
        "ligand_results": successful,
        "parameters": params,
        "grid_config": grid_config,
        "timestamp": datetime.now().isoformat(),
        "software": "AutoDock Vina"
    }
    
    data_path = outputs_dir / "virtual_screening_metadata.json"
    with open(data_path, "w") as f:
        json.dump(final_data, f, indent=2)
        
    print(f"\n[VS Analytics] ✅ Docking complete")
    print(f"[VS Analytics] ✓ Results saved to {data_path}")

if __name__ == "__main__":
    run_vina()