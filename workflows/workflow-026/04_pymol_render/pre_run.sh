#!/bin/bash
set -e
mamba install -y -c conda-forge pymol-open-source=3.1.0
pip install py3Dmol matplotlib numpy
