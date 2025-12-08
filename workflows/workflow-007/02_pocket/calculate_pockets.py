# Part 2: Calculate pockets and save to JSON
# github repo: https://github.com/cch1999/pocketeer
# doc: https://pocketeer.readthedocs.io/en/latest/

import os

import pocketeer as pt

# =============================================================================
# CONFIGURATION FROM ENVIRONMENT VARIABLES
# =============================================================================

# PDB ID (from global workflow parameter)
pdb_id = os.environ.get("PARAM_PDB_ID", "4TOS")
pdb_filename = f"{pdb_id.upper()}.pdb"

# Pocketeer find_pockets parameters (from job parameters)
r_min = float(os.environ.get("PARAM_R_MIN", "3.0"))
r_max = float(os.environ.get("PARAM_R_MAX", "6.0"))
polar_probe_radius = float(os.environ.get("PARAM_POLAR_PROBE_RADIUS", "1.8"))
merge_distance = float(os.environ.get("PARAM_MERGE_DISTANCE", "1.2"))
min_spheres = int(os.environ.get("PARAM_MIN_SPHERES", "35"))
ignore_hydrogens = os.environ.get("PARAM_IGNORE_HYDROGENS", "true").lower() == "true"
ignore_water = os.environ.get("PARAM_IGNORE_WATER", "true").lower() == "true"
ignore_hetero = os.environ.get("PARAM_IGNORE_HETERO", "true").lower() == "true"

# Output JSON file for pocket data
output_json = "pockets.json"

# =============================================================================

# Load structure
print("Loading structure...", flush=True)
atomarray = pt.load_structure(pdb_filename)

# Detect pockets with configurable parameters
print(f"Detecting pockets with parameters:", flush=True)
print(f"  r_min={r_min}, r_max={r_max}", flush=True)
print(f"  polar_probe_radius={polar_probe_radius}, merge_distance={merge_distance}", flush=True)
print(f"  min_spheres={min_spheres}", flush=True)
print(f"  ignore_hydrogens={ignore_hydrogens}, ignore_water={ignore_water}, ignore_hetero={ignore_hetero}", flush=True)

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
print(f"\nFound {len(pockets)} pockets:")
for pocket in pockets[:5]:  # Show top 5
    print(
        f"  Pocket {pocket.pocket_id}: score={pocket.score:.2f}, "
        f"volume={pocket.volume:.1f} Å³, "
        f"spheres={pocket.n_spheres}"
    )

# Save pockets to JSON
print(f"\nSaving pockets to {output_json}...", flush=True)
pt.write_pockets_json(output_json, pockets)
print(f"Pockets saved to {output_json}", flush=True)
