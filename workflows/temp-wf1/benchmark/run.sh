#!/bin/bash
set -e

URL="https://www.mpinat.mpg.de/3925290/hif2a_eq.zip"
FILENAME="hif2a_eq.zip"

echo "Attempting to download $FILENAME from $URL..."
# Download the zip file using wget
# -q: quiet (suppress output)
# -O: output document to filename
if ! wget -q -O "$FILENAME" "$URL"; then
    echo "Error: Failed to download $FILENAME from $URL"
    exit 1
fi
echo "$FILENAME downloaded successfully."

echo "Attempting to unzip $FILENAME..."
# Unzip the downloaded file
# -o: overwrite files without prompting
# -d: extract files into the specified directory (optional, current directory by default)
if ! unzip -o "$FILENAME"; then
    echo "Error: Failed to unzip $FILENAME"
    exit 1
fi
echo "$FILENAME unzipped successfully."
echo "You can now find the extracted files in the current directory."
