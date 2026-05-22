#!/bin/bash
set -e

# Copy input files from inputs/ to current directory
cp inputs/*.bam inputs/*.gtf . 2>/dev/null || true

# Create outputs directory
mkdir -p outputs

Rscript featurecounts.R

# Move outputs
mv *.csv outputs/ 2>/dev/null || true
mv *.txt outputs/ 2>/dev/null || true

# Clean up copied input files
rm -f *.bam *.gtf
