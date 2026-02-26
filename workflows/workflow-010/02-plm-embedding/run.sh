#!/bin/bash
set -e
python3 node2.py --inputs inputs --outputs outputs --device "${PARAM_DEVICE:-cuda}"
