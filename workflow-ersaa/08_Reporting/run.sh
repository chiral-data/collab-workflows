#!/bin/bash
set -e

LOG_DIR="/workspace/input"

python3 /workspace/rank_vina.py "$LOG_DIR"
