import json
import os
import sys
import urllib.request

uniprot_id = os.environ.get("PARAM_UNIPROT_ID", "P69905")

# ── FASTA sequence ─────────────────────────────────────────────────────────────
fasta_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
fasta_file = f"{uniprot_id}.fasta"
print(f"Downloading sequence for {uniprot_id} from UniProt...", flush=True)
urllib.request.urlretrieve(fasta_url, fasta_file)
print(f"Saved {fasta_file}", flush=True)

# ── Experimental PDB structure ─────────────────────────────────────────────────
print(f"\nLooking up experimental PDB structures for {uniprot_id}...", flush=True)
info_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
with urllib.request.urlopen(info_url) as resp:
    data = json.load(resp)

pdb_entries = []
for ref in data.get("uniProtKBCrossReferences", []):
    if ref.get("database") != "PDB":
        continue
    pdb_id = ref["id"]
    props = {p["key"]: p["value"] for p in ref.get("properties", [])}
    method = props.get("Method", "")
    try:
        resolution = float(props.get("Resolution", "999").replace(" A", "").strip())
    except ValueError:
        resolution = 999.0
    pdb_entries.append((pdb_id, method, resolution))

if not pdb_entries:
    print(f"Warning: no PDB structures found for {uniprot_id}, skipping ground truth.", flush=True)
    sys.exit(0)

# Prefer X-ray or EM; among those pick highest resolution (lowest Å value)
preferred = [e for e in pdb_entries if "X-ray" in e[1] or "EM" in e[1]]
best_pdb_id, best_method, best_res = sorted(preferred or pdb_entries, key=lambda e: e[2])[0]
print(f"Selected: {best_pdb_id} ({best_method}, {best_res:.2f} Å) from {len(pdb_entries)} available", flush=True)

pdb_url = f"https://files.rcsb.org/download/{best_pdb_id}.pdb"
pdb_file = f"{uniprot_id}_reference.pdb"
print(f"Downloading {best_pdb_id} from RCSB...", flush=True)
try:
    urllib.request.urlretrieve(pdb_url, pdb_file)
    print(f"Saved {pdb_file}", flush=True)
except Exception as e:
    print(f"Warning: could not download {best_pdb_id}: {e} — skipping ground truth.", flush=True)
