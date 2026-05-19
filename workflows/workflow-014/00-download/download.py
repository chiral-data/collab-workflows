import os
import urllib.request

os.makedirs("outputs", exist_ok=True)

url = "https://raw.githubusercontent.com/swansonk14/admet_ai/main/admet_ai/resources/data/drugbank_approved.csv"
dst = os.path.join("outputs", "drugbank_approved.csv")
print(f"Downloading drugbank_approved.csv ...", flush=True)
urllib.request.urlretrieve(url, dst)
print(f"Saved {dst}", flush=True)
