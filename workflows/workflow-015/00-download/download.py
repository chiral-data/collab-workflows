import os
import urllib.request

os.makedirs("outputs", exist_ok=True)

url = "https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-015/input_files/sample_proteins.fasta"
dst = os.path.join("outputs", "sample_proteins.fasta")
print(f"Downloading sample_proteins.fasta ...", flush=True)
urllib.request.urlretrieve(url, dst)
print(f"Saved {dst}", flush=True)
