import urllib.request

url = "https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-015/input_files/sample_proteins.fasta"
dst = "sample_proteins.fasta"
print(f"Downloading {dst} ...", flush=True)
urllib.request.urlretrieve(url, dst)
print(f"Saved {dst}", flush=True)
