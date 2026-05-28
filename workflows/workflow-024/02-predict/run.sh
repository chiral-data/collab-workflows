#!/bin/bash
set -e
export LOGNAME="${LOGNAME:-user}"
export NUMBA_CACHE_DIR="/tmp/numba_cache"
python3 predict.py
