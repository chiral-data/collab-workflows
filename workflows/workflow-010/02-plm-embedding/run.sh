#!/bin/bash
set -e
python3 generate_embeddings.py --inputs inputs --outputs outputs --device "${PARAM_DEVICE:-cuda}" --use-plm "${PARAM_USE_PLM:-0}"
