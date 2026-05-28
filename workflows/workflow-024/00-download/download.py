import os
import urllib.error
import urllib.request
import yaml

uniprot_id = os.environ.get("PARAM_UNIPROT_ID", "P69905")
fasta_url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"

print(f"Downloading {uniprot_id} from UniProt...", flush=True)
try:
    response = urllib.request.urlopen(fasta_url)
    fasta_text = response.read().decode("utf-8").strip()
except urllib.error.HTTPError as e:
    raise SystemExit(f"Error: could not fetch UniProt ID '{uniprot_id}' (HTTP {e.code}). Check the accession ID is valid.")

lines = fasta_text.splitlines()
sequence = "".join(l.strip() for l in lines if not l.startswith(">")).upper()

if not sequence:
    raise ValueError(f"No sequence found for UniProt ID: {uniprot_id}")

os.makedirs("./outputs", exist_ok=True)
output_path = f"./outputs/{uniprot_id}.yaml"
with open(output_path, "w") as f:
    yaml.dump({"sequences": [{"id": "A", "type": "protein", "sequence": sequence}]}, f, default_flow_style=False)

print(f"Saved {output_path} ({len(sequence)} residues)", flush=True)
