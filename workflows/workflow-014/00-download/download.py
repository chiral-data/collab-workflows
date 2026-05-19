import urllib.request

url = "https://raw.githubusercontent.com/swansonk14/admet_ai/main/admet_ai/resources/data/drugbank_approved.csv"
dst = "drugbank_approved.csv"
print("Downloading drugbank_approved.csv ...", flush=True)
urllib.request.urlretrieve(url, dst)
print(f"Saved {dst}", flush=True)
