#!/bin/bash
set -e
docker build -t lightdock-pipeline:latest .
echo "Built lightdock-pipeline:latest"
