#!/usr/bin/env python3
import argparse
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import yaml


def find_diffdock_repo(user_path: str = "") -> Tuple[Optional[str], List[str]]:
    candidates = []
    if user_path:
        candidates.append(user_path)

    env_path = os.environ.get("DIFFDOCK_PP_PATH", "").strip()
    if env_path:
        candidates.append(env_path)

    conda_prefix = os.environ.get("CONDA_PREFIX", "").strip()
    if conda_prefix:
        candidates.append(os.path.join(conda_prefix, "DiffDock-PP"))

    candidates.extend(
        [
            "/opt/DiffDock-PP",
            "/workspace/DiffDock-PP",
        ]
    )

    candidates.extend(
        [
            "./DiffDock-PP",
            "../DiffDock-PP",
            str(Path.home() / "DiffDock-PP"),
            str(Path(__file__).resolve().parent / "DiffDock-PP"),
            str(Path(__file__).resolve().parent.parent / "DiffDock-PP"),
        ]
    )

    normalized = []
    seen = set()
    for candidate in candidates:
        repo = os.path.abspath(candidate)
        if repo in seen:
            continue
        seen.add(repo)
        normalized.append(repo)

    for repo in normalized:
        main_inf = os.path.join(repo, "src", "main_inf.py")
        if os.path.isdir(repo) and os.path.isfile(main_inf):
            return repo, normalized
    return None, normalized


def setup_db5_layout(receptor_pdb: str, ligand_pdb: str, output_dir: str) -> Dict[str, str]:
    temp_root = os.path.join(output_dir, "temp_inference_data")
    structures = os.path.join(temp_root, "structures")
    os.makedirs(structures, exist_ok=True)

    complex_name = "complex"
    rec_dst = os.path.join(structures, f"{complex_name}_r_b.pdb")
    lig_dst = os.path.join(structures, f"{complex_name}_l_b.pdb")

    shutil.copy2(receptor_pdb, rec_dst)
    shutil.copy2(ligand_pdb, lig_dst)

    split_csv = os.path.join(temp_root, "splits_test.csv")
    with open(split_csv, "w", encoding="utf-8") as f:
        f.write("path,split\n")
        f.write(f"{complex_name},test\n")

    return {
        "temp_root": temp_root,
        "structures": structures,
        "split_csv": split_csv,
        "complex_name": complex_name,
        "receptor_copy": rec_dst,
        "ligand_copy": lig_dst,
    }


def create_inference_yaml(
    output_dir: str,
    split_csv: str,
    data_path: str,
    num_samples: int,
    inference_steps: int,
) -> str:
    config = {
        "data": {
            "dataset": "db5",
            "data_file": split_csv,
            "data_path": data_path,
            "resolution": "residue",
            "no_graph_cache": True,
            "knn_size": 20,
            "use_orientation_features": False,
            "multiplicity": 1,
            "use_unbound": False,
        },
        "model": {
            "model_type": "e3nn",
            "no_torsion": True,
            "no_batch_norm": True,
            "lm_embed_dim": 1280,
            "dropout": 0.0,
            "dynamic_max_cross": True,
            "cross_cutoff_weight": 3,
            "cross_cutoff_bias": 40,
            "cross_max_dist": 80,
            "num_conv_layers": 4,
            "ns": 16,
            "nv": 4,
            "dist_embed_dim": 32,
            "cross_dist_embed_dim": 32,
            "sigma_embed_dim": 32,
            "max_radius": 5.0,
        },
        "diffusion": {
            "tr_s_min": 0.01,
            "tr_s_max": 30.0,
            "rot_s_min": 0.01,
            "rot_s_max": 1.65,
            "sample_train": True,
            "num_inference_complexes_train_data": max(100, num_samples),
        },
        "inference": {
            "mirror_ligand": False,
            "run_inference_without_confidence_model": False,
            "wandb_sweep": False,
            "no_final_noise": True,
            "temp_sampling": 2.439,
            "temp_psi": 0.216,
            "temp_sigma_data_tr": 0.593,
            "temp_sigma_data_rot": 0.228,
        },
        "train": {
            "num_steps": inference_steps,
        },
    }

    yaml_path = os.path.join(output_dir, "inference_config.yaml")
    with open(yaml_path, "w", encoding="utf-8") as f:
        yaml.safe_dump(config, f, sort_keys=False)
    return yaml_path


def find_model_path(base_dir: str, model_name: str) -> str:
    fold_path = os.path.join(base_dir, "checkpoints", model_name, "fold_0")
    if os.path.isdir(fold_path):
        return fold_path

    fallback = os.path.join(base_dir, "checkpoints", model_name)
    if os.path.isdir(fallback):
        return fallback

    raise FileNotFoundError(
        f"Missing DiffDock-PP checkpoint for {model_name}. Expected {fold_path} or {fallback}."
    )


def run_diffdock(
    diffdock_repo: str,
    config_yaml: str,
    output_dir: str,
    num_samples: int,
    inference_steps: int,
    batch_size: int,
    seed: int,
    use_gpu: bool,
) -> Tuple[str, str]:
    output_dir = os.path.abspath(output_dir)
    config_yaml = os.path.abspath(config_yaml)
    visualization_dir = os.path.join(output_dir, "poses_raw")
    storage_dir = os.path.join(output_dir, "diffdock_storage")
    prediction_storage = os.path.join(storage_dir, "predictions.pkl")
    os.makedirs(visualization_dir, exist_ok=True)
    os.makedirs(storage_dir, exist_ok=True)

    score_model = find_model_path(diffdock_repo, "large_model_dips")
    confidence_model = find_model_path(diffdock_repo, "confidence_model_dips")

    cmd = [
        sys.executable,
        os.path.join(diffdock_repo, "src", "main_inf.py"),
        "--mode",
        "test",
        "--config_file",
        config_yaml,
        "--run_name",
        "inference_run",
        "--save_path",
        storage_dir,
        "--batch_size",
        str(batch_size),
        "--seed",
        str(seed),
        "--score_model_path",
        score_model,
        "--filtering_model_path",
        confidence_model,
        "--num_samples",
        str(num_samples),
        "--actual_steps",
        str(inference_steps),
        "--visualization_path",
        visualization_dir,
        "--visualize_n_val_graphs",
        str(num_samples),
        "--prediction_storage",
        prediction_storage,
    ]

    if use_gpu:
        cmd.extend(["--num_gpu", "1", "--gpu", "0"])
    else:
        cmd.extend(["--num_gpu", "0"])

    env = os.environ.copy()
    if not use_gpu:
        env["CUDA_VISIBLE_DEVICES"] = ""

    # Run from a writable directory. Some DiffDock dependencies write temporary
    # cache files to the current working directory at import time.
    result = subprocess.run(
        cmd,
        cwd=output_dir,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )

    log_path = os.path.join(output_dir, "inference_log.txt")
    with open(log_path, "w", encoding="utf-8") as f:
        f.write("COMMAND:\n")
        f.write(" ".join(cmd) + "\n\n")
        f.write(f"RETURN_CODE: {result.returncode}\n\n")
        f.write("STDOUT:\n")
        f.write(result.stdout or "")
        f.write("\n\nSTDERR:\n")
        f.write(result.stderr or "")

    if result.returncode != 0:
        raise RuntimeError(
            f"DiffDock-PP inference failed with code {result.returncode}. See {log_path} for details."
        )

    return visualization_dir, log_path


def parse_confidence_from_name(filename: str) -> Tuple[float, Optional[int], Optional[float]]:
    stem = Path(filename).stem

    # Expected: complex_1_0.0_0.95.pdb
    m = re.match(r"^.+_(\d+)_(-?\d+(?:\.\d+)?)_(-?\d+(?:\.\d+)?)$", stem)
    if m:
        return float(m.group(3)), int(m.group(1)), float(m.group(2))

    # Fallback: find the last floating-point number in filename.
    all_floats = re.findall(r"-?\d+(?:\.\d+)?", stem)
    if all_floats:
        try:
            return float(all_floats[-1]), None, None
        except ValueError:
            pass

    return 0.0, None, None


def gather_pose_files(visualization_dir: str) -> List[str]:
    paths = []
    if not os.path.isdir(visualization_dir):
        return paths
    for root, _, files in os.walk(visualization_dir):
        for name in files:
            if name.lower().endswith(".pdb"):
                paths.append(os.path.join(root, name))
    return sorted(paths)


def is_predicted_pose_file(path: str) -> bool:
    name = os.path.basename(path).lower()
    if "receptor" in name:
        return False
    if "ligand-gt" in name:
        return False
    if "ligand-0" in name:
        return False
    return "ligand" in name or "rank" in name


def collect_real_pose_files(visualization_dir: str, output_dir: str, diffdock_repo: str) -> List[str]:
    candidates = [
        os.path.abspath(visualization_dir),
        os.path.join(os.path.abspath(output_dir), "visualization"),
        os.path.join(os.path.abspath(diffdock_repo), "visualization"),
    ]

    all_paths: List[str] = []
    seen = set()
    for directory in candidates:
        for path in gather_pose_files(directory):
            full = os.path.abspath(path)
            if full in seen:
                continue
            seen.add(full)
            all_paths.append(full)

    preferred = [p for p in all_paths if is_predicted_pose_file(p)]
    return preferred if preferred else all_paths


def rank_and_export_poses(pose_files: List[str], output_dir: str, mode: str) -> Tuple[List[Dict[str, object]], List[float]]:
    parsed = []
    for pose in pose_files:
        conf, source_rank, rmsd = parse_confidence_from_name(os.path.basename(pose))
        parsed.append(
            {
                "source_file": pose,
                "source_name": os.path.basename(pose),
                "confidence": conf,
                "source_rank": source_rank,
                "source_rmsd": rmsd,
            }
        )

    parsed.sort(key=lambda x: x["confidence"], reverse=True)

    pose_details: List[Dict[str, object]] = []
    confidence_scores: List[float] = []

    for i, item in enumerate(parsed, start=1):
        out_name = f"rank{i}.pdb"
        out_path = os.path.join(output_dir, out_name)
        shutil.copy2(item["source_file"], out_path)

        detail = {
            "rank": i,
            "pose_file": os.path.abspath(out_path),
            "pose_filename": out_name,
            "confidence": float(item["confidence"]),
            "source_name": item["source_name"],
            "source_rank": item["source_rank"],
            "source_rmsd": item["source_rmsd"],
            "mode": mode,
        }
        pose_details.append(detail)
        confidence_scores.append(float(item["confidence"]))

    return pose_details, confidence_scores


def create_mock_pose(receptor_pdb: str, ligand_pdb: str, out_pdb: str) -> None:
    with open(receptor_pdb, "r", encoding="utf-8", errors="ignore") as fr:
        rec_lines = [line for line in fr if not line.startswith("END")]
    with open(ligand_pdb, "r", encoding="utf-8", errors="ignore") as fl:
        lig_lines = [line for line in fl if not line.startswith("END")]

    with open(out_pdb, "w", encoding="utf-8") as f:
        f.writelines(rec_lines)
        f.writelines(lig_lines)
        f.write("END\n")


def generate_mock_outputs(receptor_pdb: str, ligand_pdb: str, output_dir: str, num_samples: int) -> Tuple[List[Dict[str, object]], List[float], str]:
    poses_raw = os.path.join(output_dir, "poses_raw")
    os.makedirs(poses_raw, exist_ok=True)

    mock_count = max(1, min(num_samples, 10))
    for idx in range(1, mock_count + 1):
        conf = max(0.0, 1.0 - (idx - 1) * 0.08)
        fname = f"complex_{idx}_0.0_{conf:.2f}.pdb"
        create_mock_pose(receptor_pdb, ligand_pdb, os.path.join(poses_raw, fname))

    log_path = os.path.join(output_dir, "inference_log.txt")
    log_mode = "a" if os.path.isfile(log_path) else "w"
    with open(log_path, log_mode, encoding="utf-8") as f:
        if log_mode == "a":
            f.write("\n")
        f.write("Mock inference mode enabled.\n")
        f.write("No real DiffDock-PP model execution was performed.\n")

    pose_files = gather_pose_files(poses_raw)
    pose_details, confidence_scores = rank_and_export_poses(pose_files, output_dir, mode="mock")
    return pose_details, confidence_scores, log_path


def write_outputs(
    output_dir: str,
    receptor_pdb: str,
    ligand_pdb: str,
    receptor_features: str,
    ligand_features: str,
    mode: str,
    diffdock_repo: str,
    diffdock_candidates: List[str],
    num_samples: int,
    inference_steps: int,
    batch_size: int,
    pose_details: List[Dict[str, object]],
    confidence_scores: List[float],
    log_path: str,
) -> Dict[str, object]:
    confidence_json = os.path.join(output_dir, "confidence_scores.json")
    with open(confidence_json, "w", encoding="utf-8") as f:
        json.dump(
            {
                "num_poses": len(pose_details),
                "confidence_scores": confidence_scores,
                "pose_details": pose_details,
                "mode": mode,
            },
            f,
            indent=2,
        )

    data = {
        "status": "success",
        "node": "diffdock_inference",
        "mode": mode,
        "receptor": {
            "pdb": os.path.abspath(receptor_pdb),
            "role": "receptor",
            "esm2_features": os.path.abspath(receptor_features) if receptor_features else None,
        },
        "ligand": {
            "pdb": os.path.abspath(ligand_pdb),
            "role": "ligand",
            "esm2_features": os.path.abspath(ligand_features) if ligand_features else None,
        },
        "diffdock": {
            "repo_path": os.path.abspath(diffdock_repo) if diffdock_repo else None,
            "auto_detect_candidates": diffdock_candidates,
            "num_samples": int(num_samples),
            "inference_steps": int(inference_steps),
            "batch_size": int(batch_size),
        },
        "num_poses": len(pose_details),
        "confidence_scores": confidence_scores,
        "poses": [d["pose_file"] for d in pose_details],
        "pose_details": pose_details,
        "outputs": {
            "confidence_scores_json": os.path.abspath(confidence_json),
            "data_json": os.path.abspath(os.path.join(output_dir, "data.json")),
            "inference_log": os.path.abspath(log_path),
        },
    }

    data_json = os.path.join(output_dir, "data.json")
    with open(data_json, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2)

    return data


def main() -> int:
    parser = argparse.ArgumentParser(description="Run DiffDock-PP inference and rank poses by confidence.")
    parser.add_argument("--receptor_pdb", required=True)
    parser.add_argument("--ligand_pdb", required=True)
    parser.add_argument("--receptor_features", default="")
    parser.add_argument("--ligand_features", default="")
    parser.add_argument("--diffdock_path", default="")
    parser.add_argument("--num_samples", type=int, default=10)
    parser.add_argument("--inference_steps", type=int, default=20)
    parser.add_argument("--batch_size", type=int, default=1)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--use_gpu", action="store_true")
    parser.add_argument("--allow_mock_fallback", action="store_true")
    parser.add_argument("--output_dir", default="4_diffdock_inference_outputs")
    args = parser.parse_args()

    if not os.path.isfile(args.receptor_pdb):
        raise FileNotFoundError(f"Receptor PDB not found: {args.receptor_pdb}")
    if not os.path.isfile(args.ligand_pdb):
        raise FileNotFoundError(f"Ligand PDB not found: {args.ligand_pdb}")

    args.output_dir = os.path.abspath(args.output_dir)
    os.makedirs(args.output_dir, exist_ok=True)

    db5 = setup_db5_layout(args.receptor_pdb, args.ligand_pdb, args.output_dir)
    config_yaml = create_inference_yaml(
        output_dir=args.output_dir,
        split_csv=db5["split_csv"],
        data_path=db5["temp_root"],
        num_samples=args.num_samples,
        inference_steps=args.inference_steps,
    )

    diffdock_repo, diffdock_candidates = find_diffdock_repo(args.diffdock_path)

    mode = "real"
    pose_details: List[Dict[str, object]]
    confidence_scores: List[float]
    log_path: str

    try:
        if diffdock_repo is None:
            raise FileNotFoundError(
                "DiffDock-PP repository not found. Set --diffdock_path or DIFFDOCK_PP_PATH to a valid clone."
            )

        visualization_dir, log_path = run_diffdock(
            diffdock_repo=diffdock_repo,
            config_yaml=config_yaml,
            output_dir=args.output_dir,
            num_samples=args.num_samples,
            inference_steps=args.inference_steps,
            batch_size=args.batch_size,
            seed=args.seed,
            use_gpu=args.use_gpu,
        )

        pose_files = collect_real_pose_files(
            visualization_dir=visualization_dir,
            output_dir=args.output_dir,
            diffdock_repo=diffdock_repo,
        )
        if not pose_files:
            raise RuntimeError("DiffDock-PP finished but no PDB pose files were generated.")

        pose_details, confidence_scores = rank_and_export_poses(
            pose_files=pose_files,
            output_dir=args.output_dir,
            mode="real",
        )

    except Exception as e:
        if not args.allow_mock_fallback:
            raise

        mode = "mock"
        diffdock_repo = diffdock_repo or ""
        pose_details, confidence_scores, log_path = generate_mock_outputs(
            receptor_pdb=args.receptor_pdb,
            ligand_pdb=args.ligand_pdb,
            output_dir=args.output_dir,
            num_samples=args.num_samples,
        )

        with open(log_path, "a", encoding="utf-8") as f:
            f.write("\nFallback reason:\n")
            f.write(str(e) + "\n")

    data = write_outputs(
        output_dir=args.output_dir,
        receptor_pdb=args.receptor_pdb,
        ligand_pdb=args.ligand_pdb,
        receptor_features=args.receptor_features,
        ligand_features=args.ligand_features,
        mode=mode,
        diffdock_repo=diffdock_repo or "",
        diffdock_candidates=diffdock_candidates,
        num_samples=args.num_samples,
        inference_steps=args.inference_steps,
        batch_size=args.batch_size,
        pose_details=pose_details,
        confidence_scores=confidence_scores,
        log_path=log_path,
    )

    print(json.dumps(data, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
