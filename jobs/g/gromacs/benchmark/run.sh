#!/bin/bash
set -e

SAMPLE="hif2a_eq"
URL="https://www.mpinat.mpg.de/3925290/${SAMPLE}.zip"
FILENAME="${SAMPLE}.zip"

echo "Attempting to download $FILENAME from $URL..."
# Download the zip file using wget
if ! wget -q -O "$FILENAME" "$URL"; then
    echo "Error: Failed to download $FILENAME from $URL"
    exit 1
fi
echo "$FILENAME downloaded successfully."

echo "Attempting to unzip $FILENAME..."
# Unzip the downloaded file
if ! unzip -o "$FILENAME"; then
    echo "Error: Failed to unzip $FILENAME"
    exit 1
fi
echo "$FILENAME unzipped successfully."
echo "You can now find the extracted files in the current directory."

# run gromacs benchmark
gmx mdrun -s ${SAMPLE}.tpr -nsteps 10000
