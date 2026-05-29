import json
import os
import re
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
    # "Chains" format: "A/C=2-142" or "A=1-141" — first letter is the chain to use
    chains_str = props.get("Chains", "")
    m = re.match(r'([A-Za-z0-9])', chains_str)
    chain_letter = m.group(1) if m else 'A'
    pdb_entries.append((pdb_id, method, resolution, chain_letter))

if not pdb_entries:
    print(f"Warning: no PDB structures found for {uniprot_id}, skipping ground truth.", flush=True)
    sys.exit(0)

# Prefer X-ray or EM; among those pick highest resolution (lowest Å value)
preferred = [e for e in pdb_entries if "X-ray" in e[1] or "EM" in e[1]]
best_pdb_id, best_method, best_res, best_chain = sorted(preferred or pdb_entries, key=lambda e: e[2])[0]
print(f"Selected: {best_pdb_id} chain {best_chain} ({best_method}, {best_res:.2f} Å) from {len(pdb_entries)} available", flush=True)

pdb_url = f"https://files.rcsb.org/download/{best_pdb_id}.pdb"
pdb_file = f"{uniprot_id}_reference.pdb"
print(f"Downloading {best_pdb_id} from RCSB...", flush=True)
try:
    tmp_file = pdb_file + ".tmp"
    urllib.request.urlretrieve(pdb_url, tmp_file)

    # Extract only the chain that maps to the UniProt sequence
    with open(tmp_file) as f:
        raw_lines = f.readlines()

    kept = []
    for line in raw_lines:
        record = line[:6].strip()
        if record in ('ATOM', 'ANISOU', 'TER') and len(line) > 21:
            if line[21] == best_chain:
                kept.append(line)
        elif record in ('HEADER', 'TITLE', 'REMARK', 'SEQRES', 'CRYST1'):
            kept.append(line)
    kept.append('END\n')

    os.remove(tmp_file)
    with open(pdb_file, 'w') as f:
        f.writelines(kept)

    atom_count = sum(1 for l in kept if l[:6].strip() == 'ATOM')
    print(f"Saved {pdb_file} (chain {best_chain}, {atom_count} ATOM records)", flush=True)
except Exception as e:
    print(f"Warning: could not download {best_pdb_id}: {e} — skipping ground truth.", flush=True)
    if os.path.exists(pdb_file + ".tmp"):
        os.remove(pdb_file + ".tmp")
