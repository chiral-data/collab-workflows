import json
import os
import urllib.request

record_id = os.environ.get("PARAM_ZENODO_RECORD_ID", "18301020")
gtf_url = os.environ.get(
    "PARAM_GTF_URL",
    "https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-020/input_files/Saccharomyces_cerevisiae.gtf",
)

print(f"Fetching file list from Zenodo record {record_id}...")
with urllib.request.urlopen(f"https://zenodo.org/api/records/{record_id}") as r:
    record = json.load(r)

files = record["files"]
print(f"Found {len(files)} file(s)")

for f in files:
    filename = f["key"]
    url = f["links"]["self"]
    size_mb = f["size"] / 1024 / 1024
    print(f"Downloading {filename} ({size_mb:.1f} MB)...")
    urllib.request.urlretrieve(url, filename)
    print(f"  -> {filename} done")

gtf_filename = gtf_url.split("/")[-1]
print(f"Downloading {gtf_filename} from GitHub...")
urllib.request.urlretrieve(gtf_url, gtf_filename)
print(f"  -> {gtf_filename} done")

print("All files downloaded.")
