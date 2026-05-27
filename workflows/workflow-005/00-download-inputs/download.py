import os
import urllib.request

csv_url = os.environ.get(
    "PARAM_CSV_URL",
    "https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows/workflow-005/input_files/SpikeRBD_DD.csv",
)

filename = csv_url.split("/")[-1]
print(f"Downloading {filename}...")
urllib.request.urlretrieve(csv_url, filename)
print(f"  -> {filename} done")
