#!/bin/bash
set -e
echo "Starting Node 01: Data Ingestion and Preprocessing"
python load_data.py
echo "Node 01 completed"
