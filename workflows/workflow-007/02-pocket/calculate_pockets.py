# Part 2: Calculate pockets and save to JSON
# github repo: https://github.com/cch1999/pocketeer
# doc: https://pocketeer.readthedocs.io/en/latest/

import glob
import json
import os
import re
import shutil
import sys

import pocketeer as pt

# =============================================================================
# CONFIGURATION
# =============================================================================

# Input/Output directories (silva 0.4.0+)
input_dir = "inputs"
output_dir = "outputs"

# Derive PDB IDs from *.pdb filenames in inputs/
pdb_files = sorted(glob.glob(os.path.join(input_dir, "*.pdb")))
if not pdb_files:
    print("ERROR: No .pdb files found in inputs/", flush=True)
    sys.exit(1)

# Validate filenames — reject unsafe characters that would break HTML/JS in 03-visualize
VALID_STEM = re.compile(r'^[A-Za-z0-9._-]+$')
pdb_ids = []
for pdb_file in pdb_files:
    stem = os.path.splitext(os.path.basename(pdb_file))[0]
    if not VALID_STEM.match(stem):
        print(f'ERROR: Invalid filename "{os.path.basename(pdb_file)}" — filenames must only contain alphanumeric characters, hyphens, underscores, and dots.', flush=True)
        sys.exit(1)
    pdb_ids.append(stem)

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

# Generate config.json for downstream jobs (03-visualize)
output_config = os.path.join(output_dir, "config.json")
with open(output_config, "w") as f:
    json.dump({"pdb_ids": pdb_ids}, f)
print(f"\nConfig saved to {output_config}", flush=True)
print(f"All {len(pdb_ids)} protein(s) processed.", flush=True)
