# Part 2: Calculate pockets and save to JSON
# github repo: https://github.com/cch1999/pocketeer
# doc: https://pocketeer.readthedocs.io/en/latest/

import pocketeer as pt

# =============================================================================
# CONFIGURATION OPTIONS
# =============================================================================

# PDB file to analyze (output from part 1)
pdb_filename = "4TOS.pdb"

# Output JSON file for pocket data
output_json = "pockets.json"

# =============================================================================

# Load structure
print("Loading structure...", flush=True)
atomarray = pt.load_structure(pdb_filename)

# Detect pockets
print("Detecting pockets (this may take a moment)...", flush=True)
pockets = pt.find_pockets(atomarray)

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
