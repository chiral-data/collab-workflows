# github repo: https://github.com/cch1999/pocketeer
# doc: https://pocketeer.readthedocs.io/en/latest/ doc: https://pocketeer.readthedocs.io/en/latest/

import urllib.request

import pocketeer as pt

pdb_code = "4tos"

# Download the pdb file for demonstration
pdb_filename = f"{pdb_code.upper()}.pdb"
url = f"https://files.rcsb.org/download/{pdb_code.upper()}.pdb"
urllib.request.urlretrieve(url, pdb_filename)

# Load structure
atomarray = pt.load_structure(pdb_filename)

# Detect pockets
pockets = pt.find_pockets(atomarray)

# Display results
print(f"\nFound {len(pockets)} pockets:")
for pocket in pockets[:5]:  # Show top 5
    print(
        f"  Pocket {pocket.pocket_id}: score={pocket.score:.2f}, "
        f"volume={pocket.volume:.1f} Å³, "
        f"spheres={pocket.n_spheres}"
    )

# Create visualization
viewer = pt.view_pockets(atomarray, pockets)
viewer.show()
