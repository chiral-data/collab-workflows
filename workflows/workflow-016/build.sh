#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

IMAGE_NAME="diffdock_abag"
VERSION="v1"

echo "Building Docker image: ${IMAGE_NAME}:${VERSION}"
docker build --file "${SCRIPT_DIR}/Dockerfile" --tag "${IMAGE_NAME}:${VERSION}" --tag "${IMAGE_NAME}:latest" "${REPO_ROOT}"
echo "Build complete!"