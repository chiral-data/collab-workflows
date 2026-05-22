#!/bin/bash
set -e

mkdir -p outputs

Rscript deseq2.R

mv merged.csv dds.rds res.rds res_alpha.rds pca_data.csv outputs/ 2>/dev/null || true
mv deseq2_summary*.csv outputs/ 2>/dev/null || true
