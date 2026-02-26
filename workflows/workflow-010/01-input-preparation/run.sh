#!/bin/bash
set -e
python3 node1.py --heavy_fasta inputs/heavy.fasta --light_fasta inputs/light.fasta --out_dir outputs
