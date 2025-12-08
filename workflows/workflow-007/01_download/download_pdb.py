# Part 1: Download PDB file
# github repo: https://github.com/cch1999/pocketeer
# doc: https://pocketeer.readthedocs.io/en/latest/

import urllib.request

# =============================================================================
# CONFIGURATION OPTIONS
# =============================================================================

# PDB code to analyze
pdb_code = "4tos"

# =============================================================================

# Download the pdb file
print(f"Downloading PDB file: {pdb_code}...", flush=True)
pdb_filename = f"{pdb_code.upper()}.pdb"
url = f"https://files.rcsb.org/download/{pdb_code.upper()}.pdb"
urllib.request.urlretrieve(url, pdb_filename)
print(f"Downloaded {pdb_filename}", flush=True)
