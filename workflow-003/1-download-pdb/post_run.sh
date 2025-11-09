#!/bin/bash
set -e

mkdir outputs/
cp pdb_id.txt outputs/
while IFS= read -r pdb_id; do
    cp "./${pdb_id}.pdb" "./outputs/"
done <pdb_id.txt
