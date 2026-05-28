#!/bin/bash
set -e
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
docker build -t qpcr-pipeline:latest "$SCRIPT_DIR"
echo "Built qpcr-pipeline:latest"
