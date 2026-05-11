#!/bin/bash
set -e

mkdir -p outputs

Rscript plots.R

mv *.png merged_with_regulation.csv outputs/ 2>/dev/null || true
