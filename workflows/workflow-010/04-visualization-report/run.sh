#!/bin/bash
set -e
python3 generate_report.py --inputs inputs --outputs outputs --title "${PARAM_REPORT_TITLE:-ABB3 Structure Predictions}"
