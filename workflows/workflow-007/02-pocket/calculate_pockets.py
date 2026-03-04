# Part 2: Calculate pockets and save to JSON
# github repo: https://github.com/cch1999/pocketeer
# doc: https://pocketeer.readthedocs.io/en/latest/

import json
import os
import shutil

import pocketeer as pt

# =============================================================================
# CONFIGURATION
# =============================================================================

# Input/Output directories (silva 0.4.0+)
input_dir = "inputs"
output_dir = "outputs"

# Read PDB IDs from config.json (produced by 01-download)
config_path = os.path.join(input_dir, "config.json")
with open(config_path, "r") as f:
    config = json.load(f)
pdb_ids = config["pdb_ids"]

# Pocketeer find_pockets parameters (from job parameters)
r_min = float(os.environ.get("PARAM_R_MIN", "3.0"))
r_max = float(os.environ.get("PARAM_R_MAX", "6.0"))
polar_probe_radius = float(os.environ.get("PARAM_POLAR_PROBE_RADIUS", "1.8"))
merge_distance = float(os.environ.get("PARAM_MERGE_DISTANCE", "1.2"))
min_spheres = int(os.environ.get("PARAM_MIN_SPHERES", "35"))
ignore_hydrogens = os.environ.get("PARAM_IGNORE_HYDROGENS", "true").lower() == "true"
ignore_water = os.environ.get("PARAM_IGNORE_WATER", "true").lower() == "true"
ignore_hetero = os.environ.get("PARAM_IGNORE_HETERO", "true").lower() == "true"

# =============================================================================

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

print(f"Processing {len(pdb_ids)} protein(s): {', '.join(pdb_ids)}", flush=True)
print(f"Detection parameters:", flush=True)
print(f"  r_min={r_min}, r_max={r_max}", flush=True)
print(f"  polar_probe_radius={polar_probe_radius}, merge_distance={merge_distance}", flush=True)
print(f"  min_spheres={min_spheres}", flush=True)
print(f"  ignore_hydrogens={ignore_hydrogens}, ignore_water={ignore_water}, ignore_hetero={ignore_hetero}", flush=True)

for pdb_id in pdb_ids:
    print(f"\n{'='*50}", flush=True)
    print(f"Processing {pdb_id}", flush=True)
    print(f"{'='*50}", flush=True)

    pdb_filename = f"{pdb_id}.pdb"
    pdb_path = os.path.join(input_dir, pdb_filename)

    # Load structure
    print(f"Loading structure from {pdb_path}...", flush=True)
    atomarray = pt.load_structure(pdb_path)

    # Detect pockets
    pockets = pt.find_pockets(
        atomarray,
        r_min=r_min,
        r_max=r_max,
        polar_probe_radius=polar_probe_radius,
        merge_distance=merge_distance,
        min_spheres=min_spheres,
        ignore_hydrogens=ignore_hydrogens,
        ignore_water=ignore_water,
        ignore_hetero=ignore_hetero,
    )

    # Display results
    print(f"Found {len(pockets)} pockets", flush=True)
    if len(pockets) == 0:
        print("  WARNING: No pockets detected. Try adjusting parameters (lower min_spheres, wider r_min/r_max range).", flush=True)
    for pocket in pockets[:5]:
        print(
            f"  Pocket {pocket.pocket_id}: score={pocket.score:.2f}, "
            f"volume={pocket.volume:.1f} Å³, "
            f"spheres={pocket.n_spheres}"
        )

    # Save pockets to per-protein JSON
    output_json = os.path.join(output_dir, f"{pdb_id}_pockets.json")
    print(f"Saving pockets to {output_json}...", flush=True)
    pt.write_pockets_json(output_json, pockets)

    # Copy PDB file to outputs for next step
    output_pdb = os.path.join(output_dir, pdb_filename)
    shutil.copy(pdb_path, output_pdb)
    print(f"Copied PDB to {output_pdb}", flush=True)

# Forward config.json to outputs for downstream jobs
output_config = os.path.join(output_dir, "config.json")
shutil.copy(config_path, output_config)
print(f"\nConfig forwarded to {output_config}", flush=True)
print(f"All {len(pdb_ids)} protein(s) processed.", flush=True)
