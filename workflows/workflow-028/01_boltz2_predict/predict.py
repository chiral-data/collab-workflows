#!/usr/bin/env python3
"""
Validate Boltz-2 holo YAML input and run structure prediction.

Sample target: EGFR kinase domain + Lapatinib (PDB 1XKK, 2.40 Å).

Validation checks:
  - Valid YAML with required 'version' and 'sequences' fields
  - At least one protein entity with a non-empty sequence
  - At least one ligand entity with a SMILES string (holo mode required)
  - A 'properties.affinity' block naming the ligand chain (for affinity output)

Prediction:
  - Runs 'boltz predict' with --output_format mmcif
    (--output_format pdb is broken for protein-ligand complexes: boltz#298)
  - Collects *.cif, confidence_*.json, affinity_*.json, pae/pde/plddt *.npz
    to ./outputs/
"""

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

try:
    import yaml
except ImportError:
    yaml = None


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

def _parse_yaml(filepath):
    if yaml is not None:
        with open(filepath) as f:
            return yaml.safe_load(f)
    # Minimal fallback when PyYAML is unavailable
    content = Path(filepath).read_text()
    for required in ("version", "sequences", "sequence", "smiles"):
        if required not in content:
            raise ValueError(f"Missing required keyword: '{required}'")
    return None  # can't return structured data without PyYAML


def validate_holo_yaml(filepath):
    """
    Validate a Boltz-2 holo prediction YAML.
    Returns a summary dict; raises ValueError / FileNotFoundError on failure.
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {filepath}")
    if path.suffix not in (".yaml", ".yml"):
        raise ValueError(f"Expected .yaml file, got: {path.suffix}")

    data = _parse_yaml(filepath)

    summary = {"filename": path.name, "entities": [], "has_affinity_block": False}

    if data is None:
        # PyYAML unavailable — passed minimal keyword check
        print("WARNING: PyYAML not available; only basic format check performed")
        return summary

    if not isinstance(data, dict):
        raise ValueError("YAML root must be a mapping")
    if "version" not in data:
        raise ValueError("Missing required field: 'version'")

    sequences = data.get("sequences", [])
    if not sequences:
        raise ValueError("'sequences' list is empty")

    valid_types = {"protein", "rna", "dna", "ligand", "ccd"}
    has_protein = False
    has_ligand = False

    for i, entry in enumerate(sequences):
        if not isinstance(entry, dict):
            raise ValueError(f"Sequence entry {i} must be a mapping")
        for entity_type, entity_data in entry.items():
            if entity_type not in valid_types:
                raise ValueError(f"Unknown entity type '{entity_type}' in entry {i}")
            info = {"type": entity_type, "id": entity_data.get("id", f"chain_{i}")}
            if entity_type == "protein":
                seq = entity_data.get("sequence", "")
                if not seq:
                    raise ValueError(f"Empty protein sequence for chain {info['id']}")
                info["length"] = len(seq)
                has_protein = True
            elif entity_type == "ligand":
                smiles = entity_data.get("smiles", "")
                if not smiles:
                    raise ValueError(f"Ligand entity {info['id']} has no 'smiles' field")
                info["smiles"] = smiles
                has_ligand = True
            summary["entities"].append(info)

    if not has_protein:
        raise ValueError("Input YAML must contain at least one 'protein' entity")
    if not has_ligand:
        raise ValueError(
            "Input YAML must contain a 'ligand' entity with a SMILES string. "
            "Holo prediction (protein + ligand together) is required for this pipeline — "
            "it significantly improves binding site geometry over apo prediction."
        )

    # Check affinity block
    properties = data.get("properties", []) or []
    for prop in properties:
        if isinstance(prop, dict) and "affinity" in prop:
            summary["has_affinity_block"] = True
            break
    if not summary["has_affinity_block"]:
        print(
            "WARNING: No 'properties.affinity' block found. "
            "Add one to get Boltz-2 binding affinity predictions (affinity_*.json)."
        )

    return summary


# ---------------------------------------------------------------------------
# Prediction
# ---------------------------------------------------------------------------

def run_boltz_predict(input_file, diffusion_samples, recycling_steps, use_msa_server, accelerator="gpu"):
    """Run 'boltz predict' and collect all outputs to ./outputs/."""
    os.makedirs("./outputs", exist_ok=True)

    cmd = [
        "boltz", "predict", input_file,
        "--output_format", "mmcif",   # pdb flag broken for protein-ligand (boltz#298)
        "--diffusion_samples", str(diffusion_samples),
        "--recycling_steps", str(recycling_steps),
        "--devices", "1",
        "--accelerator", accelerator,
        "--no_kernels",               # cuequivariance_ops_torch missing from base image; use reference impl
        "--num_workers", "0",         # Docker /dev/shm too small for multiprocessing workers
        "--max_parallel_samples", "1",  # prevents GPU OOM on consumer GPUs (8GB VRAM)
        "--override",                 # always reprocess; avoids stale cache between runs
    ]
    if str(use_msa_server).lower() == "true":
        cmd.append("--use_msa_server")

    print(f"Running: {' '.join(cmd)}", flush=True)
    result = subprocess.run(cmd, capture_output=False, text=True)
    if result.returncode != 0:
        print(f"ERROR: boltz predict exited with code {result.returncode}")
        sys.exit(1)

    # Collect outputs from Boltz results directory
    # Boltz writes to: boltz_results_{stem}/predictions/{stem}/
    stem = Path(input_file).stem
    search_roots = [
        f"boltz_results_{stem}/predictions/{stem}",
        f"boltz_results_{stem}",
    ]

    collected = []
    for root in search_roots:
        for fpath in glob.glob(f"{root}/*"):
            if os.path.isfile(fpath):
                dest = os.path.join("./outputs", os.path.basename(fpath))
                shutil.copy2(fpath, dest)
                collected.append(os.path.basename(fpath))
                print(f"  Collected: {os.path.basename(fpath)}", flush=True)

    if not collected:
        # Fallback: walk the boltz_results tree
        for dirpath, _, filenames in os.walk("."):
            if "predictions" in dirpath:
                for fname in filenames:
                    src = os.path.join(dirpath, fname)
                    dest = os.path.join("./outputs", fname)
                    shutil.copy2(src, dest)
                    collected.append(fname)
                    print(f"  Collected (fallback): {fname}", flush=True)

    print(f"Total files collected: {len(collected)}", flush=True)

    # Verify key outputs
    cif_files = glob.glob("./outputs/*_model_*.cif")
    confidence_files = glob.glob("./outputs/confidence_*.json")
    plddt_files = glob.glob("./outputs/plddt_*.npz")
    affinity_files = glob.glob("./outputs/affinity_*.json")

    print(f"  Structure files (.cif): {len(cif_files)}", flush=True)
    print(f"  Confidence JSON:        {len(confidence_files)}", flush=True)
    print(f"  pLDDT npz:              {len(plddt_files)}", flush=True)
    print(f"  Affinity JSON:          {len(affinity_files)}", flush=True)

    if not cif_files:
        print("WARNING: No .cif structure files found — check boltz output directory")
    if not plddt_files:
        print("WARNING: No plddt_*.npz files found — downstream QC node will fail")

    # Clean up intermediate directories
    for d in glob.glob("boltz_results_*"):
        if os.path.isdir(d):
            shutil.rmtree(d)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Validate and run Boltz-2 holo prediction")
    parser.add_argument("--input", required=True)
    parser.add_argument("--diffusion-samples", type=int, default=2)
    parser.add_argument("--recycling-steps", type=int, default=3)
    parser.add_argument("--use-msa-server", default="true")
    parser.add_argument("--accelerator", default="gpu")
    args = parser.parse_args()

    print(f"Validating: {args.input}", flush=True)
    summary = validate_holo_yaml(args.input)

    for entity in summary["entities"]:
        etype = entity["type"]
        eid = entity.get("id", "?")
        extra = f", length={entity['length']}" if "length" in entity else ""
        extra += f", smiles={entity.get('smiles', '')}" if "smiles" in entity else ""
        print(f"  Entity {eid}: {etype}{extra}", flush=True)

    print("Validation passed", flush=True)

    # Write summary alongside outputs for downstream nodes
    os.makedirs("./outputs", exist_ok=True)
    with open("./outputs/input_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    run_boltz_predict(
        input_file=args.input,
        diffusion_samples=args.diffusion_samples,
        recycling_steps=args.recycling_steps,
        use_msa_server=args.use_msa_server,
        accelerator=args.accelerator,
    )


if __name__ == "__main__":
    main()
