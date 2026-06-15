#!/bin/bash
set -e
docker build -t qpcr-pipeline:latest .
echo "Built qpcr-pipeline:latest"
