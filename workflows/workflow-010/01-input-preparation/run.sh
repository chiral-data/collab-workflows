#!/bin/bash
set -e
python3 prepare_input.py --heavy_fasta inputs/heavy.fasta --light_fasta inputs/light.fasta --out_dir outputs
