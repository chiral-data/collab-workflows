#!/bin/bash
set -e
python3 node4.py --inputs inputs --outputs outputs --title "${PARAM_REPORT_TITLE:-ABB3 Structure Predictions}"
