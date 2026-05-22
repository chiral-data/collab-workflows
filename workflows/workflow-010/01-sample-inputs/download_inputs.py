# Download sample input files from GitHub repo
# Repo: https://github.com/chiral-data/collab-workflows

import os
import urllib.request

REPO_BASE_URL = "https://raw.githubusercontent.com/chiral-data/collab-workflows/main/workflows"

SOURCE_MAP = {
    "FEGrow": {
        "path": "workflow-009/input_files",
        "files": ["protein.pdb", "ligand.sdf"],
    },
    "Chiral Challenge #1": {
        "path": "workflow-010/input_files",
        "files": ["2GQG_A_prepared_for_dock.pdb"],
    },
}


def main():
    source = os.environ.get("PARAM_SOURCE", "Chiral Challenge #1")

    if source not in SOURCE_MAP:
        raise ValueError(f"Unknown source: {source}. Available: {list(SOURCE_MAP.keys())}")

    config = SOURCE_MAP[source]
    print(f"Downloading input files for source: {source}", flush=True)

    for filename in config["files"]:
        url = f"{REPO_BASE_URL}/{config['path']}/{filename}"
        print(f"  Downloading {filename} from {url}...", flush=True)
        urllib.request.urlretrieve(url, filename)
        print(f"  Saved {filename}", flush=True)

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
