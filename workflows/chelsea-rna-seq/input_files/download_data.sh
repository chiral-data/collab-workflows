#!/bin/bash

RECORD_ID="18301020"

mkdir -p data
cd data

echo "Downloading files from Zenodo record: $RECORD_ID"

curl -s "https://zenodo.org/api/records/$RECORD_ID" \
  | jq -r '.files[].links.self' \
  | while read -r url; do
        # URL format: .../files/{filename}/content - extract filename from parent dir
        filename=$(basename "$(dirname "$url")")
        echo "Downloading $filename"
        curl -L -o "$filename" "$url"
    done

echo "Download complete. Files saved in ./data/"

