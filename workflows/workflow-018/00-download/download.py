import os
import urllib.request

uniprot_id = os.environ.get("PARAM_UNIPROT_ID", "P69905")
url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
output_file = f"{uniprot_id}.fasta"

print(f"Downloading {uniprot_id} from UniProt...", flush=True)
urllib.request.urlretrieve(url, output_file)
print(f"Saved {output_file}", flush=True)
